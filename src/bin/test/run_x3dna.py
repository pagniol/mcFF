import subprocess
import argparse

def run_dssr(pdb_file):
    """
    Exécute la commande dssr sur le fichier PDB spécifié.

    Args:
        pdb_file (str): Chemin vers le fichier PDB.
        output_id (str): Identifiant à utiliser pour le fichier de sortie.
        output_dir (str, optional): Répertoire de sortie. Defaults to "./my-files/stage-E24/b/in/test/data/".

    Returns:
        subprocess.CompletedProcess: Objet contenant les informations sur l'exécution de la commande.
    """

    command = f"module load dssr/2.4.6 && x3dna-dssr -i={pdb_file} --idstr=short --pairs-onlys"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        print("Commande exécutée avec succès.")
        print(result.stdout)  # Affiche la sortie standard de dssr
    else:
        print("Erreur lors de l'exécution de la commande.")
        print(result.stderr)  # Affiche les messages d'erreur

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exécuter la commande dssr sur un fichier PDB.")
    parser.add_argument("pdb_file", help="Chemin vers le fichier PDB")
    args = parser.parse_args()

    run_dssr(args.pdb_file)
