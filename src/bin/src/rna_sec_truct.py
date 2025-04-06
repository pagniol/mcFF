def get_base_pairs(dot_bracket):
    stack = []
    base_pairs = {}
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                base_pairs[start] = i
                base_pairs[i] = start
    return base_pairs

def get_stems(base_pairs):
    stems = []
    current_stem = []
    
    # Trier les paires de bases par position pour l'identification séquentielle des tiges
    sorted_pairs = sorted(base_pairs.items())
    for i, j in sorted_pairs:
        if not current_stem:
            current_stem = [(i, j)]
        else:
            last_i, last_j = current_stem[-1]
            # Vérifier si cette paire est consécutive dans la tige
            if i == last_i or j == last_j -1:
                current_stem.append((i, j))
            else:
                stems.append(current_stem)
                current_stem = [(i, j)]
    
    if current_stem:
        stems.append(current_stem)
    
    return stems

def find_paired_base(dot_bracket, i):
    """
    Trouve l'indice de la base appariée à la base à l'indice i dans une structure secondaire d'ARN au format dot-bracket.

    Args:
        dot_bracket (str): La structure secondaire d'ARN au format dot-bracket.
        i (int): L'indice de la base pour laquelle on cherche la base appariée.

    Returns:
        int: L'indice de la base appariée, ou -1 si aucune base n'est appariée.
    """

    stack = []
    for j, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(j)
        elif char == ')':
            if not stack:
                raise ValueError('Unbalanced dotb [', dot_bracket, ']') # Pas '(' correspondante
            paired_index = stack.pop()
            if paired_index == i:
                return j
            elif j == i:
                return paired_index
    
    if stack:
        raise ValueError('Unbalanced dotb [', dot_bracket, ']') # pas de ')' correspondante

    return -1

def find_loops(dotb):
    """
    Retourne une liste de paires (start, end) des positions des boucles dans la structure dot-bracket.
    Chaque boucle est définie par les positions de son premier et dernier appariement de bases (bp).
    """
    loops = []
    current_loop = []
    deb = False

    i=0
    while i < len(dotb) -1 or not deb:
        if dotb[i] == '(':
            deb = True
            break
        i +=1
    
    pos = False
    while i< len(dotb)-1 and deb:
        char= dotb[i]
        if char == '(':
            pos = True
                    
        elif char == '.' and pos == True:
            if dotb[i+1] in ['.', ')']:
                current_loop.append(i)
            elif dotb[i+1] == '(':
                current_loop.clear()
                pos = False
           
        elif char == ')' and pos == True:
            loops.append((min(current_loop), max(current_loop)))
            pos = False
            current_loop.clear()
        i += 1
    return loops

def find_junctions( stems):
    """
    Retourne une liste des positions (début, fin) des jonctions dans la structure dot-bracket.
    Une jonction relie des tiges et peut être de longueur zéro.
    
    Args:
        dotb (str): La structure secondaire au format dot-bracket.
        stems (list): Liste des tiges représentées par des listes de paires de bases.
        loops (list): Liste des boucles représentées par des paires (start, end).

    Returns:
        list: Liste des positions (début, fin) des jonctions.
    """
    junctions = []
    
    # Obtenir les positions de début et de fin des tiges
    stem_positions = [(stem[0][0], stem[-1][1]) for stem in stems]

    # Trier les positions des tiges
    stem_positions.sort()

    # Itérer à travers les tiges pour trouver les jonctions entre elles
    for i in range(len(stem_positions) - 1):
        end_prev_stem = stem_positions[i][1]
        start_next_stem = stem_positions[i + 1][0]
        
        # Si la jonction est de longueur zéro
        if end_prev_stem + 1 == start_next_stem:
            junctions.append((end_prev_stem, start_next_stem))
        else:
            # Sinon, c'est une jonction de longueur non nulle
            junctions.append((end_prev_stem + 1, start_next_stem - 1))

    return junctions


# Exemple d'utilisation
sequence = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
dot_bracket = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'

base_pairs = get_base_pairs(dot_bracket)
stems = get_stems(base_pairs)
# Générer les jonctions
junctions = find_junctions(stems)

print(f'Les boucles sont :')
print(find_loops(dot_bracket))
#pb = get_base_pairs(dot_bracket)

#print(pb)
# Affichage des résultats
print('Les tiges:')
for stem in stems:
    print(stem)

print("\nJonctions trouvées :")
print(junctions)