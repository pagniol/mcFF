import argparse
import os
from Bio import PDB
import string

def get_existing_chain_ids(structure):
    """Récupère tous les identifiants de chaînes existants dans une structure PDB."""
    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)
    return chain_ids

def generate_unique_chain_id(existing_ids):
    """Génère un identifiant de chaîne unique qui n'est pas dans existing_ids."""
    for new_id in string.ascii_uppercase + string.ascii_lowercase + string.digits:
        if new_id not in existing_ids:
            return new_id
    raise ValueError("Pas d'identifiants de chaîne uniques disponibles.")

def convert_cif_to_pdb(cif_file, pdb_file):
    """Convertit un fichier CIF en fichier PDB et corrige les erreurs d'ID de chaîne."""
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure('structure', cif_file)

    existing_chain_ids = get_existing_chain_ids(structure)

    for model in structure:
        for chain in model:
            # Si l'ID de chaîne est trop long (> 1 char), on le remplace
            if len(chain.id) > 1:
                print(f"Identifiant de chaîne trop long détecté : {chain.id}")
                new_chain_id = generate_unique_chain_id(existing_chain_ids)
                print(f"Remplacement par un identifiant unique : {new_chain_id}")
                chain.id = new_chain_id
                existing_chain_ids.add(new_chain_id)  # Assure que cet ID ne sera plus réutilisé

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def batch_convert_cif_to_pdb(source_dir, dest_dir, issues_dir):
    """Convertit tous les fichiers CIF d'un répertoire en fichiers PDB et gère les erreurs."""
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    if not os.path.exists(issues_dir):
        os.makedirs(issues_dir)

    for filename in os.listdir(source_dir):
        if filename.endswith(".cif"):
            cif_file = os.path.join(source_dir, filename)
            pdb_file = os.path.join(dest_dir, filename.replace(".cif", ".pdb"))
            try:
                convert_cif_to_pdb(cif_file, pdb_file)
                print(f"Conversion terminée : {pdb_file}")
            except Exception as e:
                print(f"Erreur lors de la conversion de {filename}: {e}")
                issue_file = os.path.join(issues_dir, filename)
                os.rename(cif_file, issue_file)
                print(f"Fichier déplacé dans le répertoire des problèmes : {issue_file}")

def get_dotb(lst, ln):

    dotb = '.' * ln
    for (i,j) in enumerate(lst):
        dotb[i-1] = i
        dotb[j-1] = j
    return dotb

if __name__ == "__main__":
    # Définir les arguments pour l'appel via la ligne de commande
    parser = argparse.ArgumentParser(description="Convertit des fichiers CIF en fichiers PDB.")
    parser.add_argument("source_dir", help="Répertoire contenant les fichiers CIF à convertir")
    parser.add_argument("dest_dir", help="Répertoire de destination pour les fichiers PDB convertis")
    parser.add_argument("issues_dir", help="Répertoire pour stocker les fichiers qui posent problème lors de la conversion")

    # Analyser les arguments fournis
    args = parser.parse_args()

    # Appel de la fonction pour convertir les fichiers
    batch_convert_cif_to_pdb(args.source_dir, args.dest_dir, args.issues_dir)
