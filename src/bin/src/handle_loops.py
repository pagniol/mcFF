
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
                raise ValueError("Structure secondaire invalide : parenthèse fermante sans ouvrante correspondante")
            paired_index = stack.pop()
            if paired_index == i:
                return j
            elif j == i:
                return paired_index
    
    if stack:
        raise ValueError("Structure secondaire invalide : parenthèse ouvrante sans fermante correspondante")

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

# Exemple d'utilisation
dotb = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'
loops = find_loops(dotb)
print("Les boucles trouvées sont :", loops)
