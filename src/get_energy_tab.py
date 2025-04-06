import argparse
import csv
import math

def compute_energy(probability, max_value):
    """Calcule l'énergie avec la formule E = -0.616 * ln(p * 66) en gérant les erreurs."""
    try:
        probability = float(probability)  # Assure que la valeur est un float
        if probability > 0:
            energy = -0.616 * math.log(probability * 66)
            if max_value is not None:
                energy = min(energy, max_value)
            return "{:.3f}".format(energy)  # Formatage en 3 décimales
        else:
            return "inf"  # Cas où p = 0
    except (ValueError, TypeError):
        return "error"  # Pour identifier une valeur mal formatée

def process_energy_table(input_file, output_file, max_value):
    """Convertit les probabilités en énergies et sauvegarde le fichier corrigé."""
    
    with open(input_file, "r") as f:
        reader = csv.reader(f)
        data = []

        for row in reader:
            values = [float(p) if p.replace(".", "", 1).isdigit() else 0 for p in row]  # Toutes les colonnes
            energies = [compute_energy(p, max_value) for p in values]
            data.append(energies)

    # Écriture du fichier de sortie
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(data)  # Écrit toutes les lignes sans en-tête

def main():
    parser = argparse.ArgumentParser(description="Convertit un tableau de probabilités en énergies.")
    parser.add_argument("input_file", help="Fichier CSV d'entrée contenant les probabilités.")
    parser.add_argument("-o", "--output_file", required=True, help="Fichier de sortie CSV des énergies.")
    parser.add_argument("-max_value", type=float, help="Valeur maximale pour l'énergie (ex: 1, 2, 3).")

    args = parser.parse_args()
    
    process_energy_table(args.input_file, args.output_file, args.max_value)
    print(f"Résultats sauvegardés dans {args.output_file}")

if __name__ == "__main__":
    main()

