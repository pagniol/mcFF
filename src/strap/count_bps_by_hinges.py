import os
import re
import glob
import json
import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# Liste des conformations acceptées
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

# Construction de la regex à partir de la liste ACCEPTED_PAIRINGS.
# On s'assure que l'annotation se termine par un espace ou la fin de la ligne.
pairing_regex = "|".join(map(re.escape, ACCEPTED_PAIRINGS))
pairing_pattern = re.compile(rf"([AUGC])-([AUGC]).*?({pairing_regex})(?:\s|$)", re.IGNORECASE)

def extract_base_pairs(file_path):
    """
    Extrait les paires de bases d'un fichier .mc-annotate.
    Seules les lignes dont l'annotation correspond à un des motifs ACCEPTED_PAIRINGS sont retenues.
    Renvoie un ensemble de tuples : (base_pair, pairing_type).
    """
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
                        base1, base2, pairing_type = match.groups()
                        base_pair = f"{base1.upper()}-{base2.upper()}"
                        pairing_type_cleaned = pairing_type.replace("pairing", "").strip()
                        base_pairs.add((base_pair, pairing_type_cleaned))
    except Exception as e:
        print(f"Erreur lors de la lecture de {file_path}: {e}")
    return base_pairs

def process_hinge_pair(args):
    """
    Traite une paire de répertoires pour identifier les hinges.
    Pour chaque fichier commun aux deux répertoires, on extrait les paires de bases,
    puis on considère comme hinge la paire qui apparaît dans les deux fichiers.
    Retourne (hinge_key, { (pair, pairing_type): count }).
    """
    dir1, dir2, root_directory, file_map = args
    hinge_counts = defaultdict(int)
    hinge_key = f"{dir1}-{dir2}"
    common_files = file_map[dir1] & file_map[dir2]
    for file_name in common_files:
        path1 = os.path.join(root_directory, dir1, file_name)
        path2 = os.path.join(root_directory, dir2, file_name)
        bp_set1 = extract_base_pairs(path1)
        bp_set2 = extract_base_pairs(path2)
        common_pairs = bp_set1 & bp_set2
        for pair in common_pairs:
            hinge_counts[pair] += 1
    return hinge_key, hinge_counts

def transform_pairing_type(pt):
    """
    Transforme l'annotation d'une paire pour la simplifier.
    Remplace les chaînes selon les règles suivantes :
      - "Ww/Ww" devient "W/W"
      - "Ww/Ws" devient "W/S"
      - "Ww/Bs" devient "W/B"
      - "Hh/Ww" devient "H/W"
      - "Ws/Hh" devient "W/H"
      - "Hh/Bs" devient "H/B"
      - "Ss/Hh" devient "S/H"
      - "Ss/Ww" devient "S/W"
      - "Wh/Ss" devient "W/S"
      - "Hh/Ss" devient "H/S"
    Puis remplace "antiparallel" par "anti" et "parallel" par "para".  
    """
    replacements = {
        "Ww/Ww": "W/W",
        "Ww/Ws": "W/S",
        "Ww/Bs": "W/B",
        "Hh/Ww": "H/W",
        "Ws/Hh": "W/H",
        "Hh/Bs": "H/B",
        "Ss/Hh": "S/H",
        "Ss/Ww": "S/W",
        "Wh/Ss": "W/S",
        "Hh/Ss": "H/S"
    }
    for old, new in replacements.items():
        pt = pt.replace(old, new)
    pt = pt.replace("antiparallel", "anti").replace("parallel", "para")
    return pt.strip()

def convert_keys_to_str(hinge_counts):
    """
    Convertit les clés du dictionnaire hinge_counts (initialement des tuples) en chaînes.
    La nouvelle clé sera du type 'A-U | W/S anti cis' par exemple.
    """
    new_hinge_counts = {}
    for hinge_key, counts in hinge_counts.items():
        new_counts = {}
        for (pair, pairing_type), count in counts.items():
            new_key = f"{pair} | {transform_pairing_type(pairing_type)}"
            new_counts[new_key] = count
        new_hinge_counts[hinge_key] = new_counts
    return new_hinge_counts

def count_hinges(root_directory, num_processes):
    """
    Parcourt les sous-répertoires de root_directory et, pour chaque paire de répertoires,
    compte les hinges (paires de bases communes).
    Utilise multiprocessing pour accélérer le traitement.
    Retourne un dictionnaire sous la forme :
      { 'dir1-dir2': { (pair, pairing_type): count, ... }, ... }
    """
    hinge_counts = defaultdict(lambda: defaultdict(int))
    directories = sorted([d for d in os.listdir(root_directory) if os.path.isdir(os.path.join(root_directory, d))])
    file_map = {d: set(os.listdir(os.path.join(root_directory, d))) for d in directories}
    tasks = [(directories[i], directories[j], root_directory, file_map)
             for i in range(len(directories)) for j in range(i + 1, len(directories))]
    num_processes = min(num_processes, len(tasks), cpu_count())
    from multiprocessing import Pool  # Assure-toi de l'importer ici si ce n'est pas déjà fait
    with Pool(processes=num_processes) as pool:
        for hinge_key, counts in tqdm(pool.imap_unordered(process_hinge_pair, tasks),
                                      total=len(tasks), desc="Traitement des hinges"):
            for pair, count in counts.items():
                hinge_counts[hinge_key][pair] += count
    return hinge_counts

def main():
    parser = argparse.ArgumentParser(description="Analyse les hinges entre motifs NCMs à partir de fichiers .mc-annotate.")
    parser.add_argument("root_directory", help="Répertoire contenant les sous-répertoires de motifs (ex: 2_2, 4_3, etc.).")
    parser.add_argument("-o", "--output", required=True, help="Fichier de sortie pour enregistrer les résultats en JSON.")
    parser.add_argument("-p", "--processes", type=int, default=max(1, cpu_count() // 2),
                        help="Nombre maximal de processus à utiliser (par défaut : moitié des cœurs CPU).")
    args = parser.parse_args()
    
    hinge_counts = count_hinges(args.root_directory, args.processes)
    # Conversion des clés (tuples) en chaînes simplifiées
    hinge_counts_str = convert_keys_to_str(hinge_counts)
    
    with open(args.output, "w") as f:
        json.dump(hinge_counts_str, f, indent=4)
    print(f"Résultat enregistré dans {args.output}")

if __name__ == "__main__":
    main()
