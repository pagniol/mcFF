import os
import shutil
from Bio.PDB import MMCIFParser

def extract_sequence_from_cif(file_path):
    """
    Extrait la séquence d'ARN à partir d'un fichier CIF en utilisant Biopython.
    
    :param file_path: Chemin du fichier CIF
    :return: Séquence d'ARN sous forme de chaîne
    """
    print(f"Lecture du fichier : {file_path}")
    parser = MMCIFParser(QUIET=True)
    sequence = []

    try:
        structure = parser.get_structure("RNA_structure", file_path)
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Vérifie si le résidu est un nucléotide (A, U, C, G)
                    if residue.resname in ['A', 'U', 'C', 'G']:
                        sequence.append(residue.resname)
    except Exception as e:
        print(f"Erreur lors du traitement du fichier {file_path} : {e}")
        return ""
    
    if sequence:
        sequence_str = ''.join(sequence)
        print(f"Séquence extraite : {sequence_str[:30]}")
        return sequence_str
    else:
        print("Aucune séquence trouvée.")
        return ""

def remove_redundant_sequences_keep_one(input_dir, output_dir):
    """
    Copie les fichiers CIF uniques (par séquence) dans un répertoire de sortie,
    tout en conservant une seule copie pour chaque séquence redondante.
    
    :param input_dir: Chemin du répertoire contenant les fichiers CIF
    :param output_dir: Chemin du répertoire de sortie
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    sequence_to_file = {}  # Dictionnaire pour associer une séquence à un fichier
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.cif'):
            file_path = os.path.join(input_dir, file_name)
            sequence = extract_sequence_from_cif(file_path)
            
            if sequence:  # Si une séquence a été trouvée
                if sequence not in sequence_to_file:
                    sequence_to_file[sequence] = file_name
                    shutil.copy(file_path, output_dir)  # Copie le fichier pour la première occurrence
            else:
                print(f"Aucune séquence valide trouvée pour {file_name}")
    
    # Affiche les correspondances fichier-séquence retenues
    if sequence_to_file:
        print("Fichiers retenus pour chaque séquence unique :")
        for seq, file in sequence_to_file.items():
            print(f"Sequence: {seq[:30]}... -> Fichier: {file}")
    else:
        print("Aucune séquence unique trouvée. Vérifiez les fichiers CIF dans le répertoire d'entrée.")

def main():
    """
    Point d'entrée du programme pour exécution en ligne de commande.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Éliminer les redondances de séquences dans une base de données CIF.")
    parser.add_argument("input_dir", help="Chemin du répertoire contenant les fichiers CIF en entrée")
    parser.add_argument("output_dir", help="Chemin du répertoire où sauvegarder les fichiers uniques")
    
    args = parser.parse_args()
    
    remove_redundant_sequences_keep_one(args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
