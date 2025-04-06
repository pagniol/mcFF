def get_all_jonctions(dot_bracket):
    """
    Identifie toutes les jonctions d'une structure d'ARN au format dot-bracket,
    y compris les paires de bases, les tiges, les inversions, les fins de tiges
    et les sections de '.' délimitées par des caractères identiques.

    Args:
        dot_bracket (str): La structure secondaire d'ARN au format dot-bracket.

    Returns:
        list: Liste des jonctions, tiges, inversions, fins de tiges, et sections de '.'.
    """
    # Initialiser les variables
    stack = []
    base_pairs = {}
    stems = []
    inversion_regions = []
    stem_endings = []
    dot_sections = []

    # Trouver toutes les paires de bases
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                base_pairs[start] = i
                base_pairs[i] = start

    # Identifier les tiges
    current_stem = []
    sorted_pairs = sorted(base_pairs.items())
    for i, j in sorted_pairs:
        if not current_stem:
            current_stem = [(i, j)]
        else:
            last_i, last_j = current_stem[-1]
            if i == last_i + 1 or j == last_j - 1:
                current_stem.append((i, j))
            else:
                stems.append(current_stem)
                current_stem = [(i, j)]
    
    if current_stem:
        stems.append(current_stem)

    # Trouver les fins de tiges
    def find_paired_base(i):
        stack = []
        for j, char in enumerate(dot_bracket):
            if char == '(':
                stack.append(j)
            elif char == ')':
                if not stack:
                    raise ValueError('Unbalanced dotb')  
                paired_index = stack.pop()
                if paired_index == i:
                    return j
                elif j == i:
                    return paired_index
        return -1

    for stem in stems:
        first_pair = stem[0]
        paired_base = find_paired_base(first_pair[0])
        
        if paired_base + 1 < len(dot_bracket):
            next_char = dot_bracket[paired_base + 1]
            if next_char == ')' or (next_char == '.' and paired_base + 2 < len(dot_bracket) and dot_bracket[paired_base + 2] == ')'):
                stem_endings.append((paired_base, paired_base + 1))

    # Trouver les régions d'inversion (type ')....(' et autres inversions avec des points)
    i = 0
    while i < len(dot_bracket) - 1:
        if dot_bracket[i] == ')' and dot_bracket[i + 1] == '(':
            inversion_regions.append((i, i + 1))  # Région inversée simple
            i += 2
        elif dot_bracket[i] == ')' and dot_bracket[i + 1] == '.':
            # Traiter les inversions avec une séquence de '.'
            j = i + 1
            while j < len(dot_bracket) and dot_bracket[j] == '.':
                j += 1
            if j < len(dot_bracket) and dot_bracket[j] == '(':
                inversion_regions.append((i + 1, j - 1))  # Inversion avec points (.) en milieu
                i = j
            else:
                i += 1
        elif dot_bracket[i] == '.' and dot_bracket[i + 1] == ')':
            # Traiter les inversions inversées avec une séquence de '.'
            j = i + 1
            while j < len(dot_bracket) and dot_bracket[j] == '.':
                j += 1
            if j < len(dot_bracket) and dot_bracket[j] == '(':
                inversion_regions.append((i, j))  # Inversion avec points (.) en milieu
                i = j
            else:
                i += 1
        else:
            i += 1

    # Trouver les sections de '.' délimitées par des caractères identiques
    i = 0
    while i < len(dot_bracket):
        if dot_bracket[i] == '.':
            start = i
            # Recherche de la fin de la séquence de '.' délimitée par le même caractère
            while i < len(dot_bracket) and dot_bracket[i] == '.':
                i += 1
            end = i
            # Vérifier si la séquence est bien délimitée par des caractères identiques
            if start > 0 and end < len(dot_bracket) and dot_bracket[start - 1] == dot_bracket[end]:
                dot_sections.append((start + 1, end))
        else:
            i += 1

    # Retourner une liste avec toutes les jonctions, tiges, inversions, fins de tiges et sections de '.'
    all_jonctions = []
    all_jonctions.extend(inversion_regions)
    all_jonctions.extend(stem_endings)
    all_jonctions.extend(dot_sections)

    return all_jonctions

# Exemple d'utilisation
def main():
    dot_bracket = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'
    jonctions = get_all_jonctions(dot_bracket)
    print(jonctions)

if __name__ == "__main__":
    main()
    

