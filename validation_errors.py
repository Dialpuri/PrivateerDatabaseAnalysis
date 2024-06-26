import os
import sys
import numpy as np
import json
from datetime import datetime
import gzip
import shutil 
import Bio
import Bio.PDB
import Bio.SeqRecord
from pprint import pprint
from tqdm import tqdm
import requests
import enum

def get_pdb_path(pdb_code: str, unzip: bool = True) -> str:
    """Get PDB path, if unzip is True, unzip and return the new path, else, return the gzipped path

    Args:
        pdb_code (str): 4 letter PDB code

    Returns:
        str: path to the PDB file

    JD
    """
    gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/pdb/pdb{pdb_code}.ent.gz"
    
    if not unzip:
        return gzipped_pdb
    
    unzipped_pdb = f"/vault/tmp_extracted_pdbs/pdb{pdb_code}.ent"

    # If the PDB doesn't exist, then try to find the mmCIF
    if not os.path.exists(gzipped_pdb):
        gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/mmCIF/{pdb_code}.cif.gz"
        unzipped_pdb = f"/vault/tmp_extracted_cifs/{pdb_code}.cif"

        if not os.path.exists(gzipped_pdb):
            return None
 
    if os.path.exists(unzipped_pdb):
        return unzipped_pdb

    with gzip.open(gzipped_pdb, 'rb') as f_in:
        with open(unzipped_pdb, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return unzipped_pdb

def get_file_paths(root_directory):
    """
    Function to find all the files in the database 
    LH
    """
    filepaths = []
    for root, dirs, files in os.walk(root_directory):
        for f in files:
            filepaths.append(os.path.join(root,f))
    return filepaths

def get_year(json_path):
    """
    Function to find the year the structure corresponding to a particular database entry was deposited
    LH
    """
    json_filename = os.path.basename(json_path)
    pdb_code = json_filename.rpartition('.')[0]
    pdb_path = get_pdb_path(pdb_code, unzip=True)
    if not pdb_path:
        return None

    try:
        if "cif" in pdb_path:
            raise RuntimeError()

        pdbheader = Bio.PDB.parse_pdb_header(pdb_path)
        datestring = pdbheader['deposition_date']
        dt = datetime.strptime(datestring, '%Y-%m-%d')
        return int(dt.year)
    except:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_code.lower()}"
        response = requests.get(url)
        json = response.json()
        deposition_date = json[pdb_code.lower()][0]["deposition_date"][:4]
        return int(deposition_date)

def get_resolution(json_path):
    """
    Function to find resolution of structure corresponding to particular database entry
    LH/JD
    """

    json_filename = os.path.basename(json_path)
    pdb_code = json_filename.rpartition('.')[0]
    pdb_path = get_pdb_path(pdb_code, unzip=True)

    if not pdb_path: 
        return None

    try:
        pdb_header = Bio.PDB.parse_pdb_header(pdb_path)
        return pdb_header['resolution']
    except:
        print(pdb_code.lower())
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pdb_code.lower()}"
        response = requests.get(url)
        json = response.json()
        resolution = json[pdb_code.lower()]["resolution"] 
        return resolution 

def validation_errors_per_year(database_root):

    files = get_file_paths(database_root)
    
    year_data = {}

    for file in tqdm(files):
        error_count = 0

        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            for k, v in glycan_data.items(): 
                if k != "N-glycan":
                    continue
                if not v: continue
                for glycan in v:
                    sugars = glycan["sugars"]
                    for sugar in sugars:
                        if sugar["diagnostic"] == "no":
                            error_count += 1

        year = get_year(file)
        if not year: continue

        year_data.setdefault(year, {"totalGlycans": 0, "totalErrors": 0})
        year_data[year]["totalGlycans"] += 1
        if error_count:
            year_data[year]["totalErrors"] += 1

    output_file = "validation_errors_per_year.json"
    sorted_keys = sorted(list(year_data))
    data_range = np.arange(float(sorted_keys[0]), float(sorted_keys[-1]), dtype=np.intc)
    for value in data_range: 
        year_data.setdefault(int(value), {"totalGlycans": 0, "totalErrors": 0})

    return dict(sorted(year_data.items()))
    
def validation_errors_per_resolution(database_root):

    files = get_file_paths(database_root)
    
    resolution_data = {}

    for file in tqdm(files):
        error_count = 0

        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            for k, v in glycan_data.items():
                if k != "N-glycan":
                    continue 
                if not v: continue
                for glycan in v:
                    sugars = glycan["sugars"]
                    for sugar in sugars:
                        if sugar["diagnostic"] == "no":
                            error_count += 1

        resolution = get_resolution(file)
        if not resolution: continue

        round_resolution = round(resolution, 1)

        resolution_data.setdefault(round_resolution, {"totalGlycans": 0, "totalErrors": 0})
        resolution_data[round_resolution]["totalGlycans"] += 1
        if error_count:
            resolution_data[round_resolution]["totalErrors"] += 1
    
    return dict(sorted(resolution_data.items()))

def conformational_errors_per_year(database_root):

    files = get_file_paths(database_root)
    
    year_data = {}

    for file in tqdm(files):
        check_count = 0
        all_count = 0 

        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            for k, v in glycan_data.items(): 
                if k != "N-glycan":
                    continue

                if not v: continue

                for glycan in v:
                    sugars = glycan["sugars"]
                    for sugar in sugars:
                        if sugar["diagnostic"] == "check":
                            check_count += 1
                        all_count += 1

        year = get_year(file)
        if not year: continue

        year_data.setdefault(year, {"totalCheck": 0, "total": 0})
        year_data[year]["totalCheck"] += check_count
        year_data[year]["total"] += all_count

    sorted_keys = sorted(list(year_data))
    data_range = np.arange(float(sorted_keys[0]), float(sorted_keys[-1]), dtype=np.intc)
    for value in data_range: 
        year_data.setdefault(int(value), {"totalCheck": 0, "total": 0})

    return dict(sorted(year_data.items()))

def conformational_errors_per_resolution(database_root):

    files = get_file_paths(database_root)
    
    resolution_data = {}

    for file in tqdm(files):
        check_count = 0
        all_count = 0 
        
        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            for k, v in glycan_data.items(): 

                if k != "N-glycan":
                    continue
                    
                if not v: continue

                for glycan in v:
                    sugars = glycan["sugars"]
                    for sugar in sugars:
                        if sugar["diagnostic"] == "check":
                            check_count += 1
                        all_count += 1


        resolution = get_resolution(file)
        if not resolution: continue

        round_resolution = round(resolution, 1)

        resolution_data.setdefault(round_resolution, {"totalCheck": 0, "total": 0})
        resolution_data[round_resolution]["totalCheck"] += check_count
        resolution_data[round_resolution]["total"] += all_count
    
    return dict(sorted(resolution_data.items()))


def calculate_all_validation(database_root): 
    files = get_file_paths(database_root)
    
    resolution_data_val = {}
    resolution_data_conf = {}

    year_data_val = {}
    year_data_conf = {}

    for file in tqdm(files):
        check_count = 0
        no_count = 0         
        total_count = 0

        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            for k, v in glycan_data.items(): 

                if not v: continue
                if k != "n-glycan": continue

                for glycan in v:
                    sugars = glycan["sugars"]
                    for sugar in sugars:
                        if sugar["diagnostic"] == "check":
                            check_count += 1

                        if sugar["diagnostic"] == "no":
                            no_count += 1

                        total_count += 1

        year = get_year(file)
        if year: 
            year_data_val.setdefault(year, {"totalSugars": 0, "totalCheck": 0, "totalNo": 0})
            year_data_val[year]["totalSugars"] += total_count
            year_data_val[year]["totalCheck"] += check_count
            year_data_val[year]["totalNo"] += no_count

        resolution = get_resolution(file)
        if resolution: 
            round_resolution = round(resolution, 1)
            resolution_data_conf.setdefault(round_resolution, {"totalSugars": 0, "totalCheck": 0, "totalNo": 0})
            resolution_data_conf[round_resolution]["totalSugars"] += total_count
            resolution_data_conf[round_resolution]["totalCheck"] += check_count
            resolution_data_conf[round_resolution]["totalNo"] += no_count
    
    sorted_keys_val = sorted(list(year_data_val))
    data_range = np.arange(float(sorted_keys_val[0]), float(sorted_keys_val[-1]), dtype=np.intc)
    for value in data_range: 
        year_data_val.setdefault(int(value), {"totalSugars": 0, "totalCheck": 0, "totalNo": 0})
        
    return dict(sorted(year_data_val.items())), dict(sorted(resolution_data_conf.items()))

    # output_file = "validation_errors_per_year.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump(dict(sorted(year_data_val.items())), output_file, indent=4, sort_keys=True)

    # output_file = "conformational_errors_per_year.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump(dict(sorted(year_data_conf.items())), output_file, indent=4, sort_keys=True)

    # output_file = "validation_errors_per_resolution.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump(dict(sorted(resolution_data_val.items())), output_file, indent=4, sort_keys=True)     

    # output_file = "conformational_errors_per_resolution.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump(dict(sorted(resolution_data_conf.items())), output_file, indent=4, sort_keys=True)


class StatType(enum.Enum): 
    validation_errors_per_year = 1
    conformation_errors_per_year = 2
    validation_errors_per_resolution = 3
    conformation_errors_per_resolution = 4


if __name__ == "__main__":

    pdb_database_root = "/vault/privateer_database/pdb" 
    pdbredo_database_root = "/vault/privateer_database/pdbredo" 

    pdb_year_data_val, pdb_res_data_conf = calculate_all_validation(pdb_database_root)
    pdbredo_year_data_val, pdbredo_res_data_conf = calculate_all_validation(pdbredo_database_root)

    output_file = "validation_errors_per_year.json"
    with open(output_file, "w") as output_file: 
        json.dump({"pdb": pdb_year_data_val, "pdbredo": pdbredo_year_data_val}, output_file, indent=4, sort_keys=True)

    # output_file = "conformational_errors_per_year.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump({"pdb": pdb_year_data_conf, "pdbredo": pdbredo_year_data_conf}, output_file, indent=4, sort_keys=True)

    output_file = "validation_errors_per_resolution.json"
    with open(output_file, "w") as output_file: 
        json.dump({"pdb": pdb_res_data_conf, "pdbredo": pdbredo_res_data_conf}, output_file, indent=4, sort_keys=True)

    # output_file = "conformational_errors_per_resolution.json"
    # with open(output_file, "w") as output_file: 
    #     json.dump({"pdb": pdb_res_data_conf, "pdbredo": pdbredo_res_data_conf}, output_file, indent=4, sort_keys=True)

    # stats = [StatType.conformation_errors_per_resolution]

    # for stat in StatType:
    #     if stat == StatType.validation_errors_per_year:
    #         year_data = {}
    #         year_data["pdbredo"] = validation_errors_per_year(pdbredo_database_root)
    #         year_data["pdb"] = validation_errors_per_year(pdb_database_root)

    #         output_file = "validation_errors_per_year.json"
    #         with open(output_file, "w") as output_file: 
    #             json.dump(year_data, output_file, indent=4, sort_keys=True)
                
    #     elif stat == StatType.conformation_errors_per_year:
    #         year_data = {}              
    #         year_data["pdbredo"] = conformational_errors_per_year(pdbredo_database_root)
    #         year_data["pdb"] = conformational_errors_per_year(pdb_database_root)

    #         output_file = "conformational_errors_per_year.json"
    #         with open(output_file, "w") as output_file: 
    #             json.dump(year_data, output_file, indent=4, sort_keys=True)

    #     elif stat == StatType.validation_errors_per_resolution: 
    #         resolution_data = {}
    #         resolution_data["pdb"] = validation_errors_per_resolution(pdb_database_root)
    #         resolution_data["pdbredo"] = validation_errors_per_resolution(pdbredo_database_root)

    #         output_file = "validation_errors_per_resolution.json"
    #         with open(output_file, "w") as output_file: 
    #             json.dump(resolution_data, output_file, indent=4, sort_keys=True)
        
    #     elif stat == StatType.conformation_errors_per_resolution: 
    #         resolution_data = {}
    #         resolution_data["pdb"] = conformational_errors_per_resolution(pdb_database_root)
    #         resolution_data["pdbredo"] = conformational_errors_per_resolution(pdbredo_database_root)

    #         output_file = "conformational_errors_per_resolution.json"

    #         with open(output_file, "w") as output_file: 
    #             json.dump(resolution_data, output_file, indent=4, sort_keys=True)

    with open("last_updated.json", "w") as date_file: 
        x = datetime.now()

        data = {"date": f"{x.year}-{x.month}-{x.day}"}

        json.dump(data, date_file, indent=4, sort_keys=True)