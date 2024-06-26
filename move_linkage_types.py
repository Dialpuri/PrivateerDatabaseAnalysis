import os, shutil

other_folder = "linkages/other"
os.makedirs(other_folder, exist_ok=True)
NA = ["C", "U", "A", "G", "DG", "DT", "DA", "DC"]

for path in os.scandir("linkages"): 
    if "other" in path.name or os.path.isfile(path.path): 
        continue
    for subpath in os.scandir(path.path):
        
        s = subpath.name.rstrip(".json").split("-")
        s1 = s[0]
        s2 = s[-1]
        if s1 in NA or s2 in NA:
            shutil.move(subpath.path, other_folder)