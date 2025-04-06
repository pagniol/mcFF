import os
import re
import json
import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# Liste des paires de bases acceptées
ACCEPTED_PAIRINGS = [
    "Ww/Ww pairing antiparallel cis",
    "Ww/Ww pairing antiparallel trans",
    "Ww/Ws pairing antiparallel cis",
    "Ww/Ws pairing antiparallel trans",
    "Ww/Bs pairing parallel trans",
    "Hh/Ww pairing antiparallel trans",
    "Ws/Hh pairing antiparallel trans",
    "Hh/Bs pairing parallel cis",
    "Ss/Hh pairing antiparallel trans",
    "Ss/Ww pairing antiparallel cis",
    "Ws/Hh pairing antiparallel trans",
    "Wh/Ss pairing antiparallel cis",
    "Hh/Ss pairing antiparallel trans",
    "Hh/Ss pairing antiparallel cis"
]

pairing_regex = "|".join(map(re.escape, ACCEPTED_PAIRINGS))
pairing_pattern = re.compile(rf"([AUGC])-([AUGC]).*?({pairing_regex})(?:\\s|$)", re.IGNORECASE)

def extract_base_pairs(file_path):
    base_pairs = set()
    try:
        with open(file_path, "r") as f:
            inside_bp_section = False
            for line in f:
                line = line.strip()
                if line.startswith("Base-pairs"):
                    inside_bp_section = True
                    continue
                elif line.startswith("Residue conformations"):
                    inside_bp_section = False
                if inside_bp_section:
                    match = pairing_pattern.search(line)
                    if match:
                        base1, base2, _ = match.groups()
                        base_pairs.add(f"{base1.upper()}-{base2.upper()}")
    except Exception as e:
        print(f"Erreur lors de la lecture de {file_path}: {e}")
    return base_pairs

def process_hinge_pair(args):
    dir1, dir2, root_directory, file_map = args
    hinge_counts = defaultdict(int)
    hinge_key = f"{dir1}-{dir2}"
    common_files = file_map[dir1] & file_map[dir2]
    total_bps = 0
    for file_name in common_files:
        path1 = os.path.join(root_directory, dir1, file_name)
        path2 = os.path.join(root_directory, dir2, file_name)
        bp_set1 = extract_base_pairs(path1)
        bp_set2 = extract_base_pairs(path2)
        common_pairs = bp_set1 & bp_set2
        total_bps += len(common_pairs)
        for pair in common_pairs:
            hinge_counts[pair] += 1
    return hinge_key, hinge_counts, total_bps

def count_hinges(root_directory, num_processes, use_multiprocessing=True):
    hinge_counts = defaultdict(lambda: defaultdict(int))
    directories = sorted([d for d in os.listdir(root_directory) if os.path.isdir(os.path.join(root_directory, d))])
    file_map = {d: set(os.listdir(os.path.join(root_directory, d))) for d in directories}
    tasks = [(directories[i], directories[j], root_directory, file_map)
             for i in range(len(directories)) for j in range(i + 1, len(directories))]
    
    total_bps = 0
    
    if use_multiprocessing:
        num_processes = min(num_processes, len(tasks), cpu_count())
        with Pool(processes=num_processes) as pool:
            for hinge_key, counts, bps in tqdm(pool.imap_unordered(process_hinge_pair, tasks),
                                               total=len(tasks), desc="Traitement des hinges"):
                total_bps += bps
                for pair, count in counts.items():
                    hinge_counts[hinge_key][pair] += count
    else:
        for args in tqdm(tasks, desc="Traitement des hinges"):
            hinge_key, counts, bps = process_hinge_pair(args)
            total_bps += bps
            for pair, count in counts.items():
                hinge_counts[hinge_key][pair] += count
    
    print(f"Nombre total de paires de bases trouvées : {total_bps}")
    return hinge_counts

def main():
    parser = argparse.ArgumentParser(description="Analyse les hinges entre motifs NCMs à partir de fichiers .mc-annotate.")
    parser.add_argument("root_directory", help="Répertoire contenant les sous-répertoires de motifs (ex: 2_2, 4_3, etc.).")
    parser.add_argument("-o", "--output", required=True, help="Fichier de sortie pour enregistrer les résultats en JSON.")
    parser.add_argument("-p", "--processes", type=int, default=max(1, cpu_count() // 2),
                        help="Nombre maximal de processus à utiliser (par défaut : moitié des cœurs CPU).")
    parser.add_argument("--multiprocessing", action="store_true",
                        help="Activer le multiprocessing pour accélérer le traitement.")
    args = parser.parse_args()
    
    hinge_counts = count_hinges(args.root_directory, args.processes, args.multiprocessing)

    with open(args.output, "w") as f:
        json.dump(hinge_counts, f, indent=4)
    print(f"Résultat enregistré dans {args.output}")

if __name__ == "__main__":
    main()

