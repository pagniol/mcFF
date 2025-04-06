import json
import argparse

def compute_probabilities(input_file, output_file):
    """Calcule P(pair | hinge) à partir des données JSON en ignorant les hinges contenant 2_6 et 6_2."""
    
    # Charger les données depuis le fichier JSON
    with open(input_file, "r") as f:
        data = json.load(f)
    
    results = {}  # Stocker les probabilités calculées
    
    for hinge, pairs in data.items():
        # Ignorer les hinges contenant "2_6" ou "6_2"
        if "2_6" in hinge or "6_2" in hinge:
            continue
        
        total_pairs = sum(pairs.values())  # Nombre total de paires pour ce hinge
        
        if total_pairs == 0:
            continue  # Éviter la division par zéro
        
        # Calculer les probabilités P(pair | hinge)
        probabilities = {pair: round(count / total_pairs, 3) for pair, count in pairs.items()}
        results[hinge] = probabilities
    
    # Écrire les résultats dans le fichier de sortie
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Calcule P(pair | hinge) en ignorant les hinges contenant 2_6 et 6_2.")
    parser.add_argument("input_file", help="Fichier JSON contenant les statistiques des paires de bases par hinge.")
    parser.add_argument("-o", "--output_file", required=True, help="Fichier de sortie JSON pour les probabilités calculées.")
    
    args = parser.parse_args()
    
    compute_probabilities(args.input_file, args.output_file)
    print(f"Résultats sauvegardés dans {args.output_file}")

if __name__ == "__main__":
    main()

