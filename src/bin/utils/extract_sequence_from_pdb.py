from Bio.PDB import *

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
    for model in structure:
        for chain in model:
            for residue in chain:
                # Ensure we're dealing with an amino acid residue
                if residue.get_resname() not in ['H_O', 'HOH']:
                    sequence.append(residue.get_resname()) 

    return ''.join(sequence) 