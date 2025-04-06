from Bio.PDB import PDBParser

def extract_sequence_from_pdb(pdb_file):
    """
    Extracts the amino acid sequence from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        str: The amino acid sequence.
    """

    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)

    sequence = []
    residue_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Ensure we're dealing with an amino acid residue
                if residue.get_resname() not in ['H_O', 'HOH']:
                    #residue_list.append(Residue(residue.get_full_id()).toString())
                    sequence.append(residue.get_resname()) 

    return  sequence,  residue_list

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    sequence, residue_lst = extract_sequence_from_pdb(pdb_file)
    print(sequence)
    print(f'Longueur de la séquence : {len(sequence)}')
    for elm in range(len(residue_lst)) :
        print(residue_lst[elm])