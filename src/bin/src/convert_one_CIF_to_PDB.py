import argparse
from Bio import PDB

def convert_cif_to_pdb(input_cif, output_pdb):
    """
    Convertit un fichier CIF en fichier PDB.
    
    Args:
        input_cif (str): Chemin vers le fichier CIF d'entrée.
        output_pdb (str): Chemin pour enregistrer le fichier PDB de sortie.
    """
    try:
        # Créer un parseur pour lire le fichier CIF
        cif_parser = PDB.MMCIFParser(QUIET=True)
        structure = cif_parser.get_structure('structure', input_cif)
        
        # Utiliser un Writer pour écrire le fichier PDB
        pdb_writer = PDB.PDBIO()
        pdb_writer.set_structure(structure)
        pdb_writer.save(output_pdb)
        print(f"Conversion réussie : {input_cif} -> {output_pdb}")
    except Exception as e:
        print(f"Erreur lors de la conversion : {e}")

def main():
    parser = argparse.ArgumentParser(description="Convertir un fichier CIF en PDB.")
    parser.add_argument("input_cif", help="Chemin du fichier CIF à convertir.")
    parser.add_argument("output_pdb", help="Chemin pour enregistrer le fichier PDB converti.")
    
    args = parser.parse_args()
    convert_cif_to_pdb(args.input_cif, args.output_pdb)

if __name__ == "__main__":
    main()
