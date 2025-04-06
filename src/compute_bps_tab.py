import csv
import argparse
import numpy as np

def read_bp_counts(input_file):
    """Lit le fichier CSV et retourne un dictionnaire des comptes des paires de bases."""
    bases = ["A", "C", "G", "U"]
    bp_counts = { (b1, b2): 0 for b1 in bases for b2 in bases }
    total_pairs = 0

    with open(input_file, newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Ignorer l'en-tête

        for row in reader:
            pair, count, _ = row
            count = int(count)
            b1, b2 = pair.split("-")

            # Ne pas fusionner (b1, b2) et (b2, b1)
            if (b1, b2) in bp_counts:
                bp_counts[(b1, b2)] += count
            total_pairs += count

    return bp_counts, total_pairs

def compute_probabilities(bp_counts, total_pairs):
    """Calcule la matrice des probabilités des paires de bases."""
    bases = ["A", "C", "G", "U"]
    bp_probabilities = np.zeros((4, 4))

    for i, b1 in enumerate(bases):
        for j, b2 in enumerate(bases):
            bp_probabilities[i, j] = bp_counts[(b1, b2)] / total_pairs if total_pairs > 0 else 0

    return bp_probabilities

def save_matrix_to_csv(output_file, matrix):
    """Sauvegarde la matrice des probabilités dans un fichier CSV."""
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for row in matrix:
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Calcul des probabilités d'apparition des paires de bases.")
    parser.add_argument("-i", "--input", required=True, help="Fichier d'entrée contenant les paires de bases et leurs occurrences.")
    parser.add_argument("-o", "--output", required=True, help="Fichier de sortie pour enregistrer la matrice des probabilités.")
    
    args = parser.parse_args()

    # Lire les occurrences des paires
    bp_counts, total_pairs = read_bp_counts(args.input)

    # Afficher les paires trouvées et leur total
    print("\nPaires de bases trouvées et occurrences :")
    for pair, count in bp_counts.items():
        print(f"{pair}: {count}")
    print(f"\nSomme totale des paires: {total_pairs}\n")

    # Calculer les probabilités
    bp_probabilities = compute_probabilities(bp_counts, total_pairs)

    # Sauvegarder dans un fichier CSV
    save_matrix_to_csv(args.output, bp_probabilities)

    print(f"Probabilités enregistrées dans {args.output}")

if __name__ == "__main__":
    main()

