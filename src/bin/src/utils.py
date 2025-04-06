def extends(lst):
    """
    Ajoute le plus petit élément de la liste moins 1 et le plus grand élément plus 1 à la liste d'origine.
    
    Args:
        liste (list): La liste d'entiers d'entrée.

    Returns:
        list: Une nouvelle liste contenant les éléments d'origine, le plus petit - 1 et le plus grand + 1.
    """
    return [min(lst) - 1, max(lst) + 1]

# Exemple d'utilisation :
#ma_liste = [3, 7, 2, 9, 5]
#resultat = ajouter_extremes(ma_liste)
#print(resultat)  # Output: [1,10]