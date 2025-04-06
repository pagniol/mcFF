import os
import re
import glob
import csv
import argparse
from collections import defaultdict

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
    Ensuite, remplace "antiparallel" par "anti" et "parallel" par "para".
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
    pt = pt.replace("antiparallel", "anti")
    pt = pt.replace("parallel", "para")
    pt = pt.replace("pairing", '')
    return pt.strip()

def extract_base_pairs(file_path):
    """
    Extrait les paires de bases et leur type depuis un fichier .mc-annotate.
    Seules les lignes dont l'annotation correspond à un des motifs ACCEPTED_PAIRINGS sont retenues.
    Renvoie une liste de tuples : (base_pair, pairing_type).
    """
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
    pairing_regex = "|".join(map(re.escape, ACCEPTED_PAIRINGS))
    base_pair_pattern = re.compile(rf"([AUGC])-([AUGC]).*?({pairing_regex})", re.IGNORECASE)
    
    base_pairs = []
    with open(file_path, 'r') as file:
        for line in file:
            match = base_pair_pattern.search(line)
            if match:
                pair = f"{match.group(1)}-{match.group(2)}"
                pairing_type = match.group(3)
                base_pairs.append((pair, pairing_type))
    return base_pairs

def process_files(input_path):
    """
    Analyse un fichier ou un répertoire contenant des fichiers .mc-annotate
    et compte les paires de bases.
    """
    base_pair_counts = defaultdict(int)
    if os.path.isfile(input_path):
        files = [input_path]
    elif os.path.isdir(input_path):
        files = [os.path.join(root, f) 
                 for root, _, filenames in os.walk(input_path) 
                 for f in filenames if f.endswith(".mc-annotate")]
    else:
        raise ValueError(f"Erreur : {input_path} n'est ni un fichier ni un répertoire valide.")
    
    for file in files:
        base_pairs = extract_base_pairs(file)
        for pair, pairing_type in base_pairs:
            base_pair_counts[(pair, pairing_type)] += 1
    return base_pair_counts

def write_csv(output_file, base_pair_counts):
    """
    Écrit les résultats dans un fichier CSV en affichant le type de liaison simplifié.
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Pairs", "Number of occurrences", "Base Pairing Type"])
        for (pair, pairing_type), count in base_pair_counts.items():
            simplified = transform_pairing_type(pairing_type)
            writer.writerow([pair, count, simplified])

def main():
    parser = argparse.ArgumentParser(description="Compter les paires de bases et simplifier l'annotation des liaisons.")
    parser.add_argument("input", help="Fichier .mc-annotate ou répertoire contenant ces fichiers")
    parser.add_argument("-o", "--output", required=True, help="Fichier CSV de sortie")
    args = parser.parse_args()
    
    try:
        base_pair_counts = process_files(args.input)
        write_csv(args.output, base_pair_counts)
        print(f"Résultats enregistrés dans {args.output}")
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()
