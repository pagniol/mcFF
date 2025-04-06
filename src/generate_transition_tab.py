import argparse
import json

# Liste des NCMs à exclure
EXCLUDE_NCMs = {"2_6", "6_2"}

def compute_hinge_probabilities(input_file, output_file):
    # Charger le fichier JSON
    with open(input_file, "r") as f:
        data = json.load(f)

    # Filtrer les jonctions en excluant celles contenant "2_6" ou "6_2"
    filtered_data = {junc: hinges for junc, hinges in data.items() 
                     if not any(ncm in junc.split("-") for ncm in EXCLUDE_NCMs)}

    # Calculer les probabilités P(hinge | junction)
    probabilities = {}
    for junction, hinges in filtered_data.items():
        total_hinges = sum(hinges.values())  # Somme des occurrences de hinges
        
        # Normaliser chaque hinge par le total
        probabilities[junction] = {hinge: count / total_hinges for hinge, count in hinges.items()}

    # Sauvegarder les résultats dans un fichier JSON
    with open(output_file, "w") as f:
        json.dump(probabilities, f, indent=4)

    # Afficher un aperçu des résultats
    print(json.dumps(dict(list(probabilities.items())[:5]), indent=4))  # Afficher les 5 premières jonctions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculer P(hinges | junctions) à partir d'un fichier JSON.")
    parser.add_argument("input_file", help="Chemin du fichier d'entrée JSON")
    parser.add_argument("-o", "--output", required=True, help="Fichier de sortie JSON")
    args = parser.parse_args()

    compute_hinge_probabilities(args.input_file, args.output)

