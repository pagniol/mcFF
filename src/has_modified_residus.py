import os
import shutil
from Bio import PDB

def contient_residus_modifies(pdb_file):
    """
    Vérifie si un fichier PDB contient des résidus non standards.

    Args:
        pdb_file (str): Chemin vers le fichier PDB.

    Returns:
        bool: True s'il contient des résidus non standards, False sinon.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Liste des résidus standards de l'ARN
    residus_standards = {'A', 'C', 'G', 'U'}

    for model in structure:
        for chain in model:
            for residue in chain:
                if not PDB.is_aa(residue, standard=True) and residue.get_resname().strip() not in residus_standards:
                    return True
    return False

def copier_fichiers_sans_residus_modifies(repertoire_entree, repertoire_sortie):
    """
    Parcourt un répertoire de fichiers PDB et copie ceux sans résidus modifiés dans un répertoire de sortie.

    Args:
        repertoire_entree (str): Chemin vers le répertoire contenant les fichiers PDB.
        repertoire_sortie (str): Chemin vers le répertoire de sortie pour les fichiers sans résidus modifiés.
    """
    if not os.path.exists(repertoire_sortie):
        os.makedirs(repertoire_sortie)

    for fichier in os.listdir(repertoire_entree):
        if fichier.endswith('.pdb'):
            chemin_complet = os.path.join(repertoire_entree, fichier)
            if not contient_residus_modifies(chemin_complet):
                shutil.copy(chemin_complet, repertoire_sortie)
                print(f"Copié: {fichier}")
            else:
                print(f"Contient des résidus modifiés: {fichier}")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python script.py <repertoire_entree> <repertoire_sortie>")
        sys.exit(1)

    repertoire_entree = sys.argv[1]
    repertoire_sortie = sys.argv[2]

    copier_fichiers_sans_residus_modifies(repertoire_entree, repertoire_sortie)
