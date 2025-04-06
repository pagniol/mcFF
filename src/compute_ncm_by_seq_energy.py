import argparse
import csv
import math
from collections import defaultdict

# Ordre des colonnes souhaité
ORDERED_NCM_COLUMNS = [
    "3", "4", "5", "6", "2_2", "3_3", "2_3", "3_2", "2_4", "4_2", "1_2",
    "2_1", "1_3", "3_1", "2_5", "5_2", "3_4", "4_3"
]

def compute_probabilities(input_file, output_file, compute_energy):
    """Calcule P(NCM | seq) ou l'énergie associée et produit un CSV."""
    
    # Chargement des données depuis le fichier CSV
    with open(input_file, "r") as f:
        reader = csv.reader(f)
        header = next(reader)  # Lire la première ligne (nom des colonnes)
        
        # Filtrer les colonnes valides en conservant l'ordre demandé
        valid_ncm_columns = [col for col in ORDERED_NCM_COLUMNS if col in header]
        
        sequences = []
        data = defaultdict(lambda: defaultdict(int))  # {sequence: {NCM: count}}
        total_counts = defaultdict(int)  # {NCM: total_count}
        total_sequences = 0  # Nombre total de NCM observés
        
        for row in reader:
            seq = row[0]
            sequences.append(seq)
            for col, count in zip(header[1:], row[1:]):
                if col not in valid_ncm_columns:  
                    continue  # Ignorer les colonnes non requises
                
                count = int(count)
                data[seq][col] += count
                total_counts[col] += count
                total_sequences += count

    # Calcul des probabilités ou des énergies
    results = []

    for seq in sequences:
        seq_prob = 1.0  # Probabilité de la séquence (supposition d'une distribution uniforme des nucléotides)
        for nt in seq:
            seq_prob *= (1 / 4)  # Hypothèse de nucléotides (A, U, G, C) équiprobables
        
        row_result = [seq]  # Ligne pour le CSV
        for ncm in valid_ncm_columns:
            if total_counts[ncm] == 0:  
                P_ci_given_si = 0.0  # Éviter la division par zéro
            else:
                P_si_given_ci = data[seq][ncm] / total_counts[ncm]  # P(si | ci)
                P_ci = total_counts[ncm] / total_sequences  # P(ci)
                P_si = seq_prob  # P(si), basé sur les nucléotides
                
                P_ci_given_si = (P_si_given_ci * P_ci) / P_si if P_si > 0 else 0.0  # Bayes

            if compute_energy:
                # Calcul de l'énergie si P > 0, sinon valeur infinie
                energy = -0.616 * math.log(P_ci_given_si * 66) if P_ci_given_si > 0 else float("inf")
                row_result.append(round(energy, 6))  # Arrondi à 6 décimales
            else:
                row_result.append(round(P_ci_given_si, 6))  # Probabilité

        results.append(row_result)

    # Écriture des résultats dans le fichier de sortie CSV
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"] + valid_ncm_columns)  # Écrire l'en-tête
        writer.writerows(results)  # Écrire les données

def main():
    parser = argparse.ArgumentParser(description="Calcule P(NCM | seq) ou l'énergie associée en utilisant la formule de Bayes.")
    parser.add_argument("input_file", help="Fichier CSV contenant les séquences et occurrences des NCMs.")
    parser.add_argument("-o", "--output_file", required=True, help="Fichier de sortie CSV contenant les résultats.")
    parser.add_argument("-energy", action="store_true", help="Si activé, calcule l'énergie au lieu de la probabilité.")

    args = parser.parse_args()
    
    compute_probabilities(args.input_file, args.output_file, args.energy)
    print(f"Résultats sauvegardés dans {args.output_file}")

if __name__ == "__main__":
    main()

