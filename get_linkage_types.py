import os, shutil, json

other_folder = "linkages/other"



for path in os.scandir("linkages"): 
    if "other" in path.name or os.path.isfile(path.path): 
        continue
    
    dirs = []
    for subpath in os.scandir(path.path):
        with open(subpath.path) as f:
            d = json.load(f)
        
        dirs.append((subpath.name.rstrip(".json"), len(d)))

    print(f'"{path.name}": {[x[0] for x in sorted(dirs, key=lambda x: x[1], reverse=True )]},')
        