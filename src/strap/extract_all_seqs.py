import os
import csv
import argparse
from collections import Counter
from Bio.PDB import PDBParser

def extract_sequences_from_pdb(pdb_file):
    """
    Extrait les séquences de tous les modèles dans un fichier PDB.
    Retourne une liste des séquences ARN.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    sequences = []

    for model in structure:
        sequence = ""
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name in ["A", "C", "G", "U"]:  # Filtrage pour ARN uniquement
                    sequence += res_name
        if sequence:
            sequences.append(sequence)  # Ajouter la séquence extraite

    return sequences

def process_pdb_files(input_dir):
    """
    Parcourt tous les fichiers PDB de la base et extrait les séquences.
    Retourne un dictionnaire {sequence: nombre d'occurrences}.
    """
    sequence_counter = Counter()

    for root, _, files in os.walk(input_dir):
        for pdb_file in files:
            if pdb_file.endswith(".pdb"):
                pdb_path = os.path.join(root, pdb_file)
                try:
                    sequences = extract_sequences_from_pdb(pdb_path)
                    sequence_counter.update(sequences)  # Comptabilise les occurrences
                except Exception as e:
                    print(f"Erreur lors du traitement de {pdb_path}: {e}")

    return sequence_counter

def save_to_csv(sequence_counter, output_csv):
    """
    Sauvegarde les séquences et leur nombre d'occurrences dans un fichier CSV.
    """
    with open(output_csv, mode="w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["sequence", "nber_occurrence"])  # En-têtes

        for sequence, count in sequence_counter.items():
            writer.writerow([sequence, count])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extraction et comptage des séquences ARN depuis des fichiers PDB.")
    parser.add_argument("-i", "--input", required=True, help="Répertoire contenant les fichiers PDB.")
    parser.add_argument("-o", "--output", required=True, help="Fichier CSV de sortie.")
    args = parser.parse_args()

    sequence_counter = process_pdb_files(args.input)
    save_to_csv(sequence_counter, args.output)
