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

def extract_mcannotate_base_pairs(mc_annotate_file):
    """
    Extrait les paires de bases de la sortie MC-Annotate.
    
    :param mc_annotate_file: Chemin du fichier contenant le résultat de MC-Annotate
    :return: Liste de tuples représentant les paires de bases
    """
    base_pairs = []
    with open(mc_annotate_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # Chercher la section 'Base-pairs'
    inside_base_pairs = False
    for line in lines:
        line = line.strip()
        if line.startswith("Base-pairs"):
            inside_base_pairs = True
            
        if inside_base_pairs:
            # Une ligne vide signifie la fin de cette section
            if not line:
                continue
            # Extraire les paires de bases
            # Exemple de format de ligne: A1-A72 : G-C Ww/Ww pairing antiparallel cis XIX
            if ':' in line:
                parts = line.split(':')
                pair_info = parts[0].strip()  # A1-A72
                details = parts[1].strip()   # G-C Ww/Ww pairing antiparallel cis XIX
                base_pairs.append((pair_info, details))
    
    index_pairs = []
    for pair, details in base_pairs:
        # Séparer les deux parties de la paire, ex: 'A1-A72'
        nucleotide1, nucleotide2 = pair.split('-')
        # Extraire les indices numériques des nucléotides
        index1 = int(nucleotide1[1:])  # Supprime le 'A' pour obtenir le numéro
        index2 = int(nucleotide2[1:])  # Supprime le 'A' pour obtenir le numéro
        # Ajouter le couple (i, j) à la liste
        index_pairs.append([index1, index2, details])

    #trier les paires
    ignore_terms = [
        "O2P/Bh adjacent_5p pairing",
        "Bs/O2' adjacent_5p pairing",
        "O2'/Ww O2'/Bs inward pairing",
        "O2'/Ww pairing",
        "O2'/Hh adjacent_5p pairing",
        "O2'/Bh pairing",
        "O2P/Bh pairing",
        'Ww/O2P pairing',
        "O2'/Hh Hh/O2P Ww/O2P pairing",
        "Hh/Hh Bh/O2P pairing parallel",
        "O2'/Bs inward pairing", #[18:57]
        "Ww/Hh pairing antiparallel", #[8:14]
        "Ww/Ww pairing parallel", #[15:48] un peu contreversé
    ]

    pairs = list()
    mc_annotate_ignore_pairs = list()
    for i, j, details in index_pairs:
        # Vérifier si la description contient un terme à ignorer
        if any(term in details for term in ignore_terms):
            mc_annotate_ignore_pairs.append((i,j))
        else:
            pairs.append((i,j))

    return pairs, mc_annotate_ignore_pairs
def extract_x3dna_base_pairs(dssr_output_content):
    """
    Extrait les paires de bases à partir de la sortie x3DNA DSSR.

    :param dssr_output_content: Contenu du fichier x3DNA DSSR sous forme de chaîne
    :return: Liste des paires sous forme [(i, j, details), ...]
    """
    base_pairs = []

    # Parcourir les lignes de la sortie
    for line in dssr_output_content.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0].isdigit():
            try:
                # Extraire les indices des bases
                i = int(''.join(filter(str.isdigit, parts[1])))  # Récupérer l'indice numérique de G1, C2, etc.
                j = int(''.join(filter(str.isdigit, parts[2])))  # Même traitement pour le second résidu
                details = f"{parts[3]} {parts[4]}"  # Combiner paire et type d'appariement
                base_pairs.append((i, j))
            except ValueError:
                # Ignore les lignes mal formées
                continue

    return base_pairs

def read_file_content(file_path):
    """
    Ouvre un fichier, lit son contenu et le retourne.

    :param file_path: Chemin vers le fichier à lire
    :return: Contenu du fichier sous forme de chaîne
    :raises: FileNotFoundError si le fichier n'existe pas
    :raises: IOError en cas de problème de lecture
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Le fichier '{file_path}' est introuvable.")
    except IOError as e:
        raise IOError(f"Erreur lors de la lecture du fichier '{file_path}': {e}")

def extract_rnaview_base_pairs(rnaview_output_file):
    """
    Extrait la liste des couples (i, j) représentant les paires de bases à partir
    du fichier de sortie RNAView.
    
    :param rnaview_output_file: Chemin vers le fichier de sortie RNAView
    :return: Liste de tuples représentant les paires de bases
    """
    base_pairs = []

    # Lire le fichier RNAView
    with open(rnaview_output_file, "r") as file:
        for line in file:
            # Identifier les lignes contenant "pair-type"
            if line.strip().startswith("pair-type"):
                parts = line.split()
                # Les indices des bases sont dans les 2ème et 3ème colonnes
                i = int(parts[2])
                j = int(parts[3])
                details = str(parts[4])
                base_pairs.append([i, j, details])

    ignore_terms = [
        'HHt', # 9:23
        'SWt', #16:59
        'WSt', #[18, 55, 'WSt']
        'HWt', #[22, 46, 'HWt']
        'WWt', #15:48
    ]

    bp = list()
    for i, j, details in base_pairs:
        # Vérifier si la description contient un terme à ignorer
        if any(term in details for term in ignore_terms):
            continue  # Ignorer cette paire
        else:
            bp.append((i,j))
    return bp            
    

def generate_dot_bracket(pairs, sequence_length):
    """
    Crée une structure dot-bracket à partir des paires de bases extraites de MC-Annotate.
    
    :param pairs: Liste des paires avec leurs détails, par exemple :
                  [((1, 72), 'G-C Ww/Ww pairing antiparallel cis XIX'), ...]
    :param sequence_length: Longueur de la séquence ARN
    :return: Chaîne dot-bracket représentant la structure secondaire
    """
    # Paires à ignorer selon les critères
    ignore_terms = [
        "O2P/Bh adjacent_5p pairing",
        "Bs/O2' adjacent_5p pairing",
        "O2'/Ww O2'/Bs inward pairing",
        "O2'/Ww pairing",
        "O2'/Hh adjacent_5p pairing",
        "O2'/Bh pairing",
        "O2P/Bh pairing",
        'Ww/O2P pairing',
        "O2'/Hh Hh/O2P Ww/O2P pairing",
        "Hh/Hh Bh/O2P pairing parallel",
        "O2'/Bs inward pairing", #18:57
        "Ww/Hh pairing antiparallel", #[8, 14 ']
        'HHt', # 9:23
        'SWt', #16:59
        'WSt', #[18, 55, 'WSt']
        'HWt', #[22, 46, 'HWt']
        'WWt', #15:48
    ]

    primary_terms = [
        'Ww/Ww pairing parallel',
        "Ww/Ww pairing antiparallel"
    ]
     # Initialiser une liste pour la structure dot-bracket
    dot_bracket = ["."] * sequence_length

    for i, j in pairs:
        # Ajouter les paires à la structure dot-bracket
        # Les indices doivent être 0-based pour la liste
        if dot_bracket[i - 1] == "." and dot_bracket[j - 1] == ".":
            dot_bracket[i - 1] = "("
            dot_bracket[j - 1] = ")"
    
    return "".join(dot_bracket)

def has_small_loops(dotb):
    """
    Retourne une liste vide si la dot-bracket ne contient pas de boucle courte(<= 2 nt)
    Retourne une liste de coupe (stsrt, end) des position des courte boucles de la dot bracket

    :param dotb: Chaîne représentant une structure dot-bracket

    """
    stack = []  # Traque les positions des parenthèses ouvrantes
    pairings = {}  # Associe les indices des parenthèses ouvrantes et fermantes

    # Étape 1 : Identifier les appariements dans la dot-bracket
    for i, char in enumerate(dotb):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening = stack.pop()
                pairings[opening] = i
                pairings[i] = opening

    # Étape 2 : Identifier les petites boucles internes
    to_remove = set()  # Indices des paires à supprimer
    for i in pairings:
        if dotb[i] == '(':
            j = pairings[i]
            # Vérifier si c'est une boucle interne (contient uniquement des '.')
            if all(dotb[k] == '.' for k in range(i + 1, j)):
                loop_size = j - i - 1
                if loop_size <= 2:  # Boucle trop courte
                    to_remove.add((i,j)) #(deb, fin) de la boucle
    #On retourner la liste triée
    if to_remove:
        return sorted(to_remove)
    else:
        return list()

def main(pdb, mcannotate_f, rnaview_f, x3dna_f):

    
    sequence = extract_sequence_from_pdb(pdb)

    mcannotate_pairs, ignore_pairs = extract_mcannotate_base_pairs(mcannotate_f)

    print(f'Les paires trouvées par mc-annotate:')
    for pair in mcannotate_pairs:
        print(pair)

    rna_view_pairs = extract_rnaview_base_pairs(rnaview_f)

    print(f'Les paires trouvées par rna-view:')
    for pair in rna_view_pairs:
        print(pair)

    x3_dna_bps = extract_x3dna_base_pairs(read_file_content(x3_dna_file))
    print(f'les paires trouvées par x3DNA:')
    for pair in x3_dna_bps:
        print(pair)

    # On combine les deux liste de paires sans les redondances et en excluant les paires ignorées par mc-annotate
    pairs = mcannotate_pairs
    for elm in rna_view_pairs:
        if (elm not in ignore_pairs) and (elm not in pairs):
            pairs.append(elm)
        else:
            ignore_pairs.append(elm)

    x3_dna_only = list()
    for elt in x3_dna_bps:
        if (elt not in pairs) and (elt not in ignore_pairs):
            x3_dna_only.append(elt)
            
    dotb = generate_dot_bracket(pairs, len(sequence))

    print(f'La séquence extriate : {sequence}')
    print(f'Les données extraites sont :')
    for pair in pairs:
        print(pair)

    print('Les paires ignorées sont:')
    for ig in ignore_pairs:
        print(ig)

    print(f'la dotB est : {dotb}')

    print(f'Liste des courtes boucles: {has_small_loops(dotb)}')

    #On elève les paires impliquées dans les courtes boucles
    while has_small_loops(dotb):
        first = has_small_loops(dotb)[0]
        for idx, (i,j) in enumerate(pairs):
            if i == first[0] + 1 :
                del pairs[idx]
        dotb = generate_dot_bracket(pairs, len(sequence))        
    print(f'la Nouvelle dotb est : {dotb}')
if __name__ == '__main__':
    mcanotate_file_path = "/home/sagniol/my-files/stage-E24/bin/test/data/1ehz.annotate"
    rna_view_file_path = "/home/sagniol/my-files/stage-E24/bin/test/data/1ehz.rnaview.out"
    x3_dna_file = "/home/sagniol/my-files/stage-E24/bin/test/data/1ehz.x3dna.out"

    pdb_file = "/home/sagniol/my-files/stage-E24/bin/test/data/1evv_A.pdb"

    main(pdb_file, mcanotate_file_path, rna_view_file_path, x3_dna_file)
