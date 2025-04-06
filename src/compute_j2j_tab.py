import json
import argparse
import numpy as np

# Définition des NCMs d'intérêt
NCM_ORDER = [
    "3", "4", "5", "6", "2_2", "3_3", "2_3", "3_2", "2_4", "4_2", 
    "1_2", "2_1", "1_3", "3_1", "2_5", "5_2", "3_4", "4_3", "4_4", "3_5", "5_3"
]

def read_junction_counts(input_file):
    """Lit le fichier JSON et retourne un dictionnaire des comptes des jonctions."""
    junction_counts = { (ncm1, ncm2): 0 for ncm1 in NCM_ORDER for ncm2 in NCM_ORDER }
    total_pairs = 0

    with open(input_file, "r") as jsonfile:
        data = json.load(jsonfile)

        for junction, bp_dict in data.items():
            ncm1, ncm2 = junction.split("-")

            # Vérifier si les NCMs sont dans la liste d'intérêt
            if ncm1 in NCM_ORDER and ncm2 in NCM_ORDER:
                count = sum(bp_dict.values())  # Total des paires pour cette jonction
                junction_counts[(ncm1, ncm2)] += count
                total_pairs += count  # Ajouter au total général

    return junction_counts, total_pairs

def compute_probabilities(junction_counts, total_pairs):
    """Calcule la matrice des probabilités des jonctions."""
    size = len(NCM_ORDER)
    junction_probabilities = np.zeros((size, size))

    for i, ncm1 in enumerate(NCM_ORDER):
        for j, ncm2 in enumerate(NCM_ORDER):
            junction_probabilities[i, j] = junction_counts[(ncm1, ncm2)] / total_pairs if total_pairs > 0 else 0

    return junction_probabilities

def save_matrix_to_csv(output_file, matrix):
    """Sauvegarde la matrice des probabilités dans un fichier CSV."""
    with open(output_file, "w", newline="") as csvfile:
        for row in matrix:
            csvfile.write(",".join(map(str, row)) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Calcul des probabilités d'apparition des jonctions entre NCMs.")
    parser.add_argument("-i", "--input", required=True, help="Fichier JSON contenant les jonctions et leurs occurrences.")
    parser.add_argument("-o", "--output", required=True, help="Fichier de sortie pour enregistrer la matrice des probabilités.")
    
    args = parser.parse_args()

    # Lire les occurrences des jonctions
    junction_counts, total_pairs = read_junction_counts(args.input)

    # Afficher les jonctions trouvées et leur total
    print("\nJonctions trouvées et occurrences :")
    for junction, count in junction_counts.items():
        if count > 0:
            print(f"{junction}: {count}")
    print(f"\nSomme totale des paires dans les jonctions valides: {total_pairs}\n")

    # Calculer les probabilités
    junction_probabilities = compute_probabilities(junction_counts, total_pairs)

    # Sauvegarder dans un fichier CSV
    save_matrix_to_csv(args.output, junction_probabilities)

    print(f"Probabilités des jonctions enregistrées dans {args.output}")

if __name__ == "__main__":
    main()

