import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import time
import multiprocessing

def process_pdb_file(pdb_file, motif_script):
    """
    Fonction pour traiter un fichier PDB en exécutant mcsearch.
    """
    command = ["/u/sagnioln/stage-E24/tools/mcsearch", motif_script, pdb_file]
    try:
        subprocess.run(command, capture_output=False, text=True, check=True)
        return f"Succès: {pdb_file}"
    except subprocess.CalledProcessError as e:
        return f"Échec: {pdb_file} - {e}"

def main(directory, motif_script, num_jobs):
    start_time = time.time()

    # Obtenir tous les fichiers PDB dans le répertoire
    pdb_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".pdb")]
    total_files = len(pdb_files)

    if total_files == 0:
        print("Aucun fichier PDB trouvé.")
        return

    print(f"Nombre total de fichiers PDB à traiter : {total_files}")
    print(f"Utilisation de {num_jobs} jobs en parallèle.")

    # Utilisation du parallélisme avec un nombre défini de processus
    with ProcessPoolExecutor(max_workers=num_jobs) as executor:
        futures = {executor.submit(process_pdb_file, pdb, motif_script): pdb for pdb in pdb_files}
        for future in tqdm(as_completed(futures), total=total_files, desc="Progression", unit="file"):
            result = future.result()
            print(result)

    end_time = time.time()
    print(f"Traitement terminé en {end_time - start_time:.2f} secondes.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recherche de motifs dans les fichiers PDB en parallèle.")
    parser.add_argument("motif_script", help="Fichier du script de motif (.mcs)")
    parser.add_argument("directory", help="Répertoire contenant les fichiers PDB")
    parser.add_argument("--num_jobs", type=int, default=multiprocessing.cpu_count(),
                        help="Nombre de processus à exécuter en parallèle (défaut: nombre de cœurs CPU)")

    args = parser.parse_args()
    main(args.directory, args.motif_script, args.num_jobs)
