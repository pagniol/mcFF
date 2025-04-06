import os
import argparse
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
from Bio.PDB import PDBParser

def extract_sequences_from_pdb(pdb_file):
    """
    Extrait les séquences de tous les modèles dans un fichier PDB.

    :param pdb_file: chemin vers le fichier PDB
    :return: liste de séquences (une par modèle)
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    except Exception as e:
        print(f"Erreur lors de l'analyse de {pdb_file}: {e}")
        return []
    
    model_sequences = []

    for model in structure:
        chain_sequences = []
        for chain in model:
            chain_seq = ""
            for residue in chain:
                res_id = residue.get_id()[1]  # Extraction du numéro de résidu (ex: 271 dans 'A271')
                if isinstance(res_id, int):  # Vérifier que c'est bien un entier
                    chain_seq += residue.get_resname().strip()
            if chain_seq:
                chain_sequences.append(chain_seq)
        model_sequences.extend(chain_sequences)

    return model_sequences



def count_occurrences(ncm_type, ncm_sequences, query_sequences):
    """
    Compte les occurrences de chaque séquence dans un NCM donné.

    :param ncm_type: Nom du NCM ('n_m' pour double-strands ou 'p' pour single-strands).
    :param ncm_sequences: Liste des séquences contenues dans ce NCM.
    :param query_sequences: Liste des séquences recherchées.
    :return: Dictionnaire {séquence: nombre d'occurrences dans ce NCM}
    """
    # Vérifier si c'est un double-strand ou un single-strand
    if "_" in ncm_type:  # Cas des NCM de type "n_m"
        try:
            n, m = map(int, ncm_type.split("_"))
            ncm_length = n + m
        except ValueError:
            print(f"⚠️ Ignoring invalid NCM type: {ncm_type}")
            return {}
    else:  # Cas des single-strands (ex: "3", "4", etc.)
        try:
            ncm_length = int(ncm_type)
        except ValueError:
            print(f"⚠️ Ignoring invalid NCM type: {ncm_type}")
            return {}

    # Initialiser le comptage des occurrences
    sequence_counts = {seq: 0 for seq in query_sequences if len(seq) == ncm_length}

    for seq in ncm_sequences:
        for query in sequence_counts:
            sequence_counts[query] += seq.count(query)

    return sequence_counts


def main():
    parser = argparse.ArgumentParser(description="Analyse les occurrences des séquences NCM dans des fichiers PDB")
    parser.add_argument("-d", "--database", required=True, help="Répertoire contenant les NCMs")
    parser.add_argument("-s", "--sequences", required=True, help="Fichier contenant les séquences à rechercher")
    parser.add_argument("-o", "--output", required=True, help="Fichier CSV de sortie")
    parser.add_argument("-n", "--num_workers", type=int, default=4, help="Nombre de cœurs pour le multiprocessing")

    args = parser.parse_args()

    # Charger les séquences à rechercher
    with open(args.sequences, "r") as f:
        query_sequences = [line.strip() for line in f.readlines()]

    # Récupérer les types de NCM disponibles
    ncm_types = [d for d in os.listdir(args.database) if os.path.isdir(os.path.join(args.database, d))]

    # Construire les tâches pour multiprocessing
    tasks = []
    for ncm in ncm_types:
        pdb_files = [os.path.join(args.database, ncm, f) for f in os.listdir(os.path.join(args.database, ncm)) if f.endswith(".pdb")]
        ncm_sequences = []
        for pdb in pdb_files:
            ncm_sequences.extend(extract_sequences_from_pdb(pdb))
        tasks.append((ncm, ncm_sequences, query_sequences))

    # Exécution parallèle avec barre de progression
    with Pool(processes=args.num_workers) as pool:
        results = list(tqdm(pool.starmap(count_occurrences, tasks), total=len(tasks), desc="Analyse des NCMs"))

    # Fusion des résultats
    final_counts = {seq: {ncm: 0 for ncm in ncm_types} for seq in query_sequences}
    for ncm, result in zip(ncm_types, results):
        for seq, count in result.items():
            final_counts[seq][ncm] = count

    # Sauvegarde en CSV
    df = pd.DataFrame.from_dict(final_counts, orient="index")
    df.to_csv(args.output)

    print(f"✅ Analyse terminée. Résultats enregistrés dans {args.output}")


if __name__ == "__main__":
    main()
