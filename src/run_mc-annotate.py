import os
import sys
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_pdb_file(pdb_path, mc_annotate_executable):
    """
    Exécute la commande mc-annotate sur un fichier PDB et enregistre la sortie dans un fichier.
    
    :param pdb_path: Chemin complet vers le fichier PDB.
    :param mc_annotate_executable: Chemin vers l'exécutable mc-annotate.
    :return: Message de résultat pour le fichier traité.
    """
    command = [mc_annotate_executable, pdb_path]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        output = result.stdout
        # Générer le nom de fichier de sortie : même nom que pdb avec extension .mc-annotate
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        output_file = os.path.join(os.path.dirname(pdb_path), base_name + ".mc-annotate")
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(output)
        return f"[OK] {pdb_path} -> {output_file}"
    except subprocess.CalledProcessError as e:
        return f"[ERREUR] {pdb_path} : {e}"
    except Exception as e:
        return f"[ERREUR] {pdb_path} : {e}"

def process_directory(input_dir, mc_annotate_executable, num_workers):
    """
    Parcourt un répertoire et exécute mc-annotate sur tous les fichiers PDB en parallèle.
    
    :param input_dir: Répertoire contenant les fichiers PDB.
    :param mc_annotate_executable: Chemin vers l'exécutable mc-annotate.
    :param num_workers: Nombre de processus parallèles à utiliser.
    """
    if not os.path.isdir(input_dir):
        print(f"Erreur : {input_dir} n'est pas un répertoire valide.")
        sys.exit(1)
    
    pdb_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                 if os.path.isfile(os.path.join(input_dir, f)) and f.lower().endswith(".pdb")]
    
    total_files = len(pdb_files)
    if total_files == 0:
        print("Aucun fichier PDB trouvé dans le répertoire spécifié.")
        return

    print(f"Nombre total de fichiers PDB à traiter: {total_files}")
    print(f"Utilisation de {num_workers} processus en parallèle.")
    
    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_pdb_file, pdb, mc_annotate_executable): pdb for pdb in pdb_files}
        for future in as_completed(futures):
            results.append(future.result())
    
    # Afficher les messages de résultat
    for res in results:
        print(res)

def main():
    parser = argparse.ArgumentParser(
        description="Exécute mc-annotate sur tous les fichiers PDB d'un répertoire en parallèle."
    )
    parser.add_argument("input_dir", type=str, help="Répertoire contenant les fichiers PDB")
    parser.add_argument("mc_annotate", type=str, help="Chemin vers l'exécutable mc-annotate")
    parser.add_argument("--num_workers", type=int, default=1,
                        help="Nombre de processus parallèles à utiliser (défaut: 1)")
    args = parser.parse_args()
    
    process_directory(args.input_dir, args.mc_annotate, args.num_workers)

if __name__ == "__main__":
    main()

