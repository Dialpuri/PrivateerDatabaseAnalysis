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
from collections import Counter

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
            filepaths.append((os.path.join(root,f), f))
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
        if not pdb_header["resolution"]: 
            raise RuntimeError()
        return pdb_header['resolution']
    except:
        
        if os.path.exists("cached_resolutions.json"):
            with open("cached_resolutions.json") as json_file:
                data = json.load(json_file)
        else:
            data = {}
        
        if pdb_code.lower() not in data.keys(): 
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pdb_code.lower()}"
            response = requests.get(url)
            json_repsonse = response.json()
            try:
                resolution = json_repsonse[pdb_code.lower()][0]["resolution"] 
            except:
                return None
            
            if resolution != None:
                data[pdb_code.lower()] = resolution
                with open("cached_resolutions.json", "w") as json_file:
                    json.dump(data, json_file)

            return resolution 
        else:
            return data[pdb_code.lower()]

def get_statistics(database_root): 
    files = get_file_paths(database_root)
    
    stats = {}

    for file_data in tqdm(files):
        file, name = file_data

        name = name[:4]
        with open(file, 'r') as f:
            data = json.load(f)
            glycan_data = data['glycans']
            n_glycans = glycan_data["n-glycan"]
            o_glycans = glycan_data["o-glycan"]
            s_glycans = glycan_data["s-glycan"]
            c_glycans = glycan_data["c-glycan"]
            ligands = glycan_data["ligand"]

            try:
                n_glycan_linkages = dict(Counter([x for g in n_glycans for x in set(list(g["linkages"].keys()))]))
                o_glycan_linkages = dict(Counter([x for g in o_glycans for x in set(list(g["linkages"].keys()))]))
                s_glycan_linkages = dict(Counter([x for g in s_glycans for x in set(list(g["linkages"].keys()))]))
                c_glycan_linkages = dict(Counter([x for g in c_glycans for x in set(list(g["linkages"].keys()))]))
                ligand_linkages = dict(Counter([x for g in ligands for x in set(list(g["linkages"].keys()))]))
            except KeyError as e: 
                print(name, e)
                continue

            stats[name] = {
                "n-glycans": len(n_glycans), 
                "o-glycans": len(o_glycans), 
                "s-glycans": len(s_glycans), 
                "c-glycans": len(c_glycans), 
                "ligands": len(ligands), 
                "path": file, 
                "n_glycan_linkages": n_glycan_linkages,
                "o_glycan_linkages": o_glycan_linkages,
                "s_glycan_linkages": s_glycan_linkages,
                "c_glycan_linkages": c_glycan_linkages,
                "ligand_linkages": ligand_linkages
            }
            
        year = get_year(file)
        stats[name]["year"] = year if year else "N/A"
        resolution = get_resolution(file)
        stats[name]["resolution"] = resolution if resolution else "N/A"

    return stats
        

class StatType(enum.Enum): 
    validation_errors_per_year = 1
    conformation_errors_per_year = 2
    validation_errors_per_resolution = 3
    conformation_errors_per_resolution = 4


if __name__ == "__main__":

    pdb_database_root = "/vault/privateer_database_backup2/pdb" 

    stats = get_statistics(pdb_database_root)

    n_glycans = {}
    o_glycans = {}
    s_glycans = {}
    c_glycans = {}
    ligands = {}

    output_dirs = [
        "linkages/n-glycan",
        "linkages/o-glycan",
        "linkages/s-glycan",
        "linkages/c-glycan",
        "linkages/ligand",
                  ]

    for o in output_dirs:
        shutil.rmtree(o)
        os.makedirs(o, exist_ok=True)
    

    for k, v in stats.items(): 
        n_glycan_linkages = v["n_glycan_linkages"]
        for k1, v1 in n_glycan_linkages.items(): 
            n_glycans.setdefault(k1, []).append({"pdb": k, "count": v1, "resolution": v["resolution"], "year": v["year"]})
        
        o_glycan_linkages = v["o_glycan_linkages"]
        for k1, v1 in o_glycan_linkages.items(): 
            o_glycans.setdefault(k1, []).append({"pdb": k, "count": v1, "resolution": v["resolution"], "year": v["year"]})
        
        s_glycan_linkages = v["s_glycan_linkages"]
        for k1, v1 in s_glycan_linkages.items(): 
            s_glycans.setdefault(k1, []).append({"pdb": k, "count": v1, "resolution": v["resolution"], "year": v["year"]})
        
        c_glycan_linkages = v["c_glycan_linkages"]
        for k1, v1 in c_glycan_linkages.items(): 
            c_glycans.setdefault(k1, []).append({"pdb": k, "count": v1, "resolution": v["resolution"], "year": v["year"]})
        
        ligand_linkages = v["ligand_linkages"]
        for k1, v1 in ligand_linkages.items(): 
            ligands.setdefault(k1, []).append({"pdb": k, "count": v1, "resolution": v["resolution"], "year": v["year"]})
    
    
    NA = ["C", "U", "A", "G", "DG", "DT", "DA", "DC"]

    for k, v in n_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            output_path: f"linkages/other/{k}.json"
        else:
            output_path = f"linkages/n-glycan/{k}.json"
            
        with open(output_path, "w") as output_file: 
            json.dump(v, output_file, indent=4, sort_keys=True)

    for k, v in o_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            output_path: f"linkages/other/{k}.json"
        else:
            output_path = f"linkages/o-glycan/{k}.json"
        with open(output_path, "w") as output_file: 
            json.dump(v, output_file, indent=4, sort_keys=True)
    
    for k, v in s_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            output_path: f"linkages/other/{k}.json"
        else:
            output_path = f"linkages/s-glycan/{k}.json"
        with open(output_path, "w") as output_file: 
            json.dump(v, output_file, indent=4, sort_keys=True)

    for k, v1 in c_glycans.items(): 
        print(k)
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            # output_path: f"linkages/other/{k}.json"
            # print(output_path, k, v1)

            # with open(output_path, "w") as output_file: 
            #     json.dump(v1, output_file, indent=4, sort_keys=True)
            continue
        else:
            output_path = f"linkages/c-glycan/{k}.json"
            print(output_path, k, v1)
            with open(output_path, "w") as output_file: 
                json.dump(v1, output_file, indent=4, sort_keys=True)
            continue
            

    for k, v in ligands.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            output_path: f"linkages/other/{k}.json"
        else:
            output_path = f"linkages/ligand/{k}.json"
        with open(output_path, "w") as output_file: 
            json.dump(v, output_file, indent=4, sort_keys=True)

    n_glycans_filtered = {}
    o_glycans_filtered = {}
    c_glycans_filtered = {}
    s_glycans_filtered = {}

    for k, v in n_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            continue
        n_glycans_filtered[k] = v

    for k, v in o_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            continue
        o_glycans_filtered[k] = v
        

    for k, v in c_glycans.items():
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            continue
        c_glycans_filtered[k] = v
        
    for k, v in s_glycans.items(): 
        s = k.split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA: 
            continue
        s_glycans_filtered[k] = v

    output_path = "linkages/n-glycan/any.json"
    with open(output_path, "w") as output_file: 
        json.dump(n_glycans_filtered, output_file, indent=4, sort_keys=True)

    output_path = "linkages/o-glycan/any.json"
    with open(output_path, "w") as output_file: 
        json.dump(o_glycans_filtered, output_file, indent=4, sort_keys=True)

    output_path = "linkages/s-glycan/any.json"
    with open(output_path, "w") as output_file: 
        json.dump(s_glycans_filtered, output_file, indent=4, sort_keys=True)

    output_path = "linkages/c-glycan/any.json"
    with open(output_path, "w") as output_file: 
        json.dump(c_glycans_filtered, output_file, indent=4, sort_keys=True)



    # output_path = "individual/individual_counts.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(stats, output_file, indent=4, sort_keys=True)

    # all_n_glycans = {}
    # all_o_glycans = {}
    # all_s_glycans = {}
    # all_c_glycans = {}
    # all_ligands = {}

    # for pdb, stats in stats.items():
    #     n_glycans = stats["n-glycans"]
    #     o_glycans = stats["o-glycans"]
    #     s_glycans = stats["s-glycans"]
    #     c_glycans = stats["c-glycans"]
    #     ligands = stats["ligands"]

    #     if n_glycans:   all_n_glycans[pdb] = {"count": stats["n-glycans"], "year": stats["year"], "resolution": stats["resolution"]}
    #     if o_glycans:   all_o_glycans[pdb] ={"count": stats["o-glycans"], "year": stats["year"], "resolution": stats["resolution"]}
    #     if s_glycans:   all_s_glycans[pdb] ={"count": stats["s-glycans"], "year": stats["year"], "resolution": stats["resolution"]}
    #     if c_glycans:   all_c_glycans[pdb] = {"count": stats["c-glycans"], "year": stats["year"], "resolution": stats["resolution"]}
    #     if ligands:   all_ligands[pdb] = {"count": stats["ligands"], "year": stats["year"], "resolution": stats["resolution"]}

    # output_path = "individual/n-glycans.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(all_n_glycans, output_file, indent=4, sort_keys=True)

    # output_path = "individual/o-glycans.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(all_o_glycans, output_file, indent=4, sort_keys=True)

    # output_path = "individual/s-glycans.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(all_s_glycans, output_file, indent=4, sort_keys=True)

    # output_path = "individual/c-glycans.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(all_c_glycans, output_file, indent=4, sort_keys=True)

    # output_path = "individual/ligands.json"
    # with open(output_path, "w") as output_file: 
    #     json.dump(all_ligands, output_file, indent=4, sort_keys=True)
