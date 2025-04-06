class Junction:
    def __init__(self, start, end, type):
        self.start = start
        self.end = end
        self.type = type
    
    def get_tuple(self, seq):
        return (self.type ,seq[slice(self.start - 1,self.end)], list(range(self.start, self.end+1)))
    
    def __repr__(self):
        return f"Junction(start={self.start}, end={self.end} type={self.type})"

def get_base_pairs(dot_bracket):
    """
    Trouve toutes les paires de bases dans une structure d'ARN au format dot-bracket.

    Args:
        dot_bracket (str): La structure secondaire d'ARN au format dot-bracket.

    Returns:
        dict: Dictionnaire avec les indices des bases appariées.
    """
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
    """
    Identifie les tiges dans une structure d'ARN en se basant sur les paires de bases.

    Args:
        base_pairs (dict): Dictionnaire des paires de bases de la structure ARN.

    Returns:
        list: Liste des tiges, chaque tige étant une liste de paires (i, j).
    """
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
            if i == last_i + 1 or j == last_j - 1:
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
                raise ValueError('Unbalanced dotb [', dot_bracket, ']')  # Pas '(' correspondante
            paired_index = stack.pop()
            if paired_index == i:
                return j
            elif j == i:
                return paired_index
    
    if stack:
        raise ValueError('Unbalanced dotb [', dot_bracket, ']')  # Pas de ')' correspondante

    return -1

def find_all_jonctions(dot_bracket):
    """
    Fonction principale qui regroupe les jonctions trouvées par les fonctions précédentes
    et retourne un dictionnaire des différents types de jonctions.

    Args:
        dot_bracket (str): Structure secondaire d'ARN au format dot-bracket.

    Returns:
        dict: Un dictionnaire contenant les jonctions par type.
    """
    #Trouve les positions où une tige se termine sans qu'une autre ne commence immédiatement après.
    base_pairs = get_base_pairs(dot_bracket)
    stems = get_stems(base_pairs)
    
    stem_endings = []
    for stem in stems:
        # Vérifier uniquement la première paire de la tige
        first_pair = stem[0]
        paired_base = find_paired_base(dot_bracket, first_pair[0])
        
        # Vérifier si le caractère après la base complémentaire est ')' ou un '.' suivi de ')'
        if paired_base + 1 < len(dot_bracket):
            next_char = dot_bracket[paired_base + 1]
            if next_char == ')' or (next_char == '.' and paired_base + 2 < len(dot_bracket) and dot_bracket[paired_base + 2] == ')'):
                stem_endings.append((paired_base + 1, paired_base + 2, 'INC')) #

    dot_sections = []
    i = 0
    while i < len(dot_bracket):
        if dot_bracket[i] == '.':
            start = i
            # Recherche d'une séquence de '.' délimitée par le même caractère
            while i < len(dot_bracket) and dot_bracket[i] == '.':
                i += 1
            end = i
            # Vérifier si la séquence est bien délimitée par des caractères identiques
            if start > 0 and end < len(dot_bracket) and dot_bracket[start - 1] == dot_bracket[end]:
                if dot_bracket[end] == '(':
                    dot_sections.append((start + 1, end,'EXC')) # 
                else:
                    dot_sections.append((start +1, end,'EXC')) # 

        else:
            i += 1

    inversion_regions = []
    # On parcourt la structure et on cherche les inversions de type ')(' et celles avec des points.
    i = 0
    while i < len(dot_bracket) - 1:
        if dot_bracket[i] == ')' and dot_bracket[i + 1] == '(':
            inversion_regions.append((i +1, i + 2,'INC')) #
            i += 2  # Passer à l'index après la paire inversée
        elif dot_bracket[i] == ')' and dot_bracket[i + 1] == '.':
            # Vérifier si une inversion existe avec des points
            j = i + 1
            while j < len(dot_bracket) and dot_bracket[j] == '.':
                j += 1
            if j < len(dot_bracket) and dot_bracket[j] == '(':
                inversion_regions.append((i+2, j,'EXC'))  #
                i = j  # Passer à la position après la parenthèse ouvrante
            else:
                i += 1
        else:
            i += 1
            
    return inversion_regions + stem_endings + dot_sections
        

# Exemple d'utilisation
def main():
    dot_bracket = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'
    junctions = find_all_jonctions(dot_bracket)
    
    print(f'Les jonctions trouvées sont : \n {junctions}')

    print(f'Les paires de bases sont :\n {get_base_pairs(dot_bracket)}' )


if __name__ == "__main__":
    main()