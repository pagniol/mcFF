import os
import csv
import argparse
import multiprocessing
from tqdm import tqdm

def parse_pdb_models(pdb_file):
    """Parse un fichier PDB et retourne une liste de numéros de résidus par modèle."""
    models = []
    current_model = []
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("MODEL"):
                    current_model = []
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    res_id = line[22:26].strip()
                    if res_id.isdigit():
                        current_model.append(int(res_id))
                elif line.startswith("ENDMDL"):
                    if current_model:
                        models.append(current_model)
        return models
    except Exception as e:
        print(f"Erreur avec le fichier {pdb_file}: {e}")
        return []

def detect_junctions(pdb_path1, pdb_path2):
    """Détecte le nombre de jonctions entre deux fichiers PDB."""
    models1 = parse_pdb_models(pdb_path1)
    models2 = parse_pdb_models(pdb_path2)

    jonctions = 0
    for model1 in models1:
        last_res = model1[-1] if model1 else None
        for model2 in models2:
            first_res = model2[0] if model2 else None
            if last_res is not None and first_res is not None and last_res + 1 == first_res:
                jonctions += 1
    return jonctions

def process_ncm_pair(ncm1, ncm2, base_path, pdb_files):
    """Compte les jonctions entre deux types de NCM."""
    total_junctions = 0
    for pdb_file in pdb_files:
        pdb_path1 = os.path.join(base_path, ncm1, pdb_file)
        pdb_path2 = os.path.join(base_path, ncm2, pdb_file)

        if os.path.exists(pdb_path1) and os.path.exists(pdb_path2):
            total_junctions += detect_junctions(pdb_path1, pdb_path2)

    return f"{ncm1}-{ncm2}", total_junctions

def main(base_path, output_file=None, num_workers=4):
    """Compte les jonctions entre toutes les paires de NCMs et les enregistre."""
    ncm_types = sorted([d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))])
    pdb_files_by_ncm = {ncm: os.listdir(os.path.join(base_path, ncm)) for ncm in ncm_types}

    tasks = []
    for i in range(len(ncm_types)):
        for j in range(i + 1, len(ncm_types)):  # Évite les doublons (ex: 2_3-5_2 == 5_2-2_3)
            ncm1, ncm2 = ncm_types[i], ncm_types[j]
            common_pdbs = set(pdb_files_by_ncm[ncm1]) & set(pdb_files_by_ncm[ncm2])  # Fichiers communs
            if common_pdbs:
                tasks.append((ncm1, ncm2, base_path, common_pdbs))

    results = []
    with multiprocessing.Pool(processes=num_workers) as pool:
        for result in tqdm(pool.starmap(process_ncm_pair, tasks), total=len(tasks), desc="Analyse des jonctions"):
            results.append(result)

    # Écriture des résultats dans un CSV
    if output_file:
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Jonction", "Occurrences"])
            writer.writerows(results)
        print(f"Résultats enregistrés dans {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compteur de jonctions entre paires de NCMs")
    parser.add_argument("base_path", help="Chemin vers la base de données contenant les répertoires de NCMs")
    parser.add_argument("-o", "--output", help="Fichier de sortie CSV", default="ncm_junctions.csv")
    parser.add_argument("-n", "--num-workers", type=int, default=4, help="Nombre de processus parallèles")

    args = parser.parse_args()
    main(args.base_path, args.output, args.num_workers)
