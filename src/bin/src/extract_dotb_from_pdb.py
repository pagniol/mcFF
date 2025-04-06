from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB import Selection

# Distance maximale pour considérer une paire
PAIRING_DISTANCE_CUTOFF = 3.0

# Appariements Watson-Crick et leurs atomes clés
PAIRING_RULES = {
    ('A', 'U'): [('N1', 'N3')],
    ('G', 'C'): [('N1', 'N3')],
    ('U', 'A'): [('N3', 'N1')],
    ('C', 'G'): [('N3', 'N1')],
}

def detect_base_pairing(residue1, residue2, neighbor_search):
    """
    Détecte si deux résidus sont appariés selon les règles définies.
    
    :param residue1: Premier résidu
    :param residue2: Deuxième résidu
    :param neighbor_search: NeighborSearch pour détecter les atomes proches
    :return: True si appariés, sinon False
    """
    base1, base2 = residue1.resname.strip(), residue2.resname.strip()
    if (base1, base2) in PAIRING_RULES:
        for atom1_name, atom2_name in PAIRING_RULES[(base1, base2)]:
            if atom1_name in residue1 and atom2_name in residue2:
                atom1 = residue1[atom1_name]
                atom2 = residue2[atom2_name]
                # Vérifie si les atomes sont suffisamment proches
                if neighbor_search.search(atom1.coord, PAIRING_DISTANCE_CUTOFF) and (
                        atom1 - atom2) <= PAIRING_DISTANCE_CUTOFF:
                    return True
    return False

def extract_rna_sequence_and_dot_bracket(pdb_file):
    """
    Extrait la séquence d'ARN et sa représentation dot-bracket d'un fichier PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA_structure", pdb_file)
    sequence, nt_indices, pairs = [], [], set()
    all_atoms = list(Selection.unfold_entities(structure, "A"))
    neighbor_search = NeighborSearch(all_atoms)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname.strip()
                if resname in ['A', 'U', 'C', 'G']:
                    sequence.append(resname)
                    nt_indices.append(residue.id[1])
                    for partner_residue in chain:
                        if partner_residue != residue:
                            if detect_base_pairing(residue, partner_residue, neighbor_search):
                                pairs.add((residue.id[1], partner_residue.id[1]))

    # Construire la dot-bracket
    sequence_length = len(sequence)
    dot_bracket = ["."] * sequence_length
    for i, nt_index in enumerate(nt_indices):
        for partner_index in [p[1] for p in pairs if p[0] == nt_index]:
            if partner_index in nt_indices:
                j = nt_indices.index(partner_index)
                dot_bracket[i] = "("
                dot_bracket[j] = ")"
    return "".join(sequence), "".join(dot_bracket)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Extrait la séquence et la structure dot-bracket d'un ARN au format PDB.")
    parser.add_argument("input_file", help="Chemin du fichier PDB en entrée")
    
    args = parser.parse_args()
    
    seq, dotb = extract_rna_sequence_and_dot_bracket(args.input_file)
    print(f"Séquence : {seq}")
    print(f"Dot-bracket : {dotb}")

if __name__ == "__main__":
    main()
