import subprocess
import argparse
from Bio.PDB import *

#from rna_sec_struc_handler import rna_sec_structure

#IMPORTANT ! La structure de ce fichier est à revoir. Il y ades attributs redondants dans la classe "rna_sec_structure" !
class rna_sec_struct_analyser:
    
    def __init__(self, pdb_f):
        self.pdb_file = pdb_f
        self.seq = self.extract_sequence_from_pdb() #À enlever
        self.dotb = self.generate_dotb() # À enlever

        #self.rna = rna_sec_structure(self.seq, self.dotb)

    def extract_sequence_from_pdb(self):
        """
        Extracts the amino acid sequence from a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.

        Returns:
            str: The amino acid sequence.
        """
        parser = PDBParser()
        structure = parser.get_structure('structure', self.pdb_file)

        sequence = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Ensure we're dealing with an amino acid residue
                    if residue.get_resname() not in ['H_O', 'HOH']:
                        sequence.append(residue.get_resname()) 

        return ''.join(sequence) 
    
    def execute_command(command, capture_output=True, shell=True):
        """
        Exécute une commande système et retourne le résultat.

        :param command: Liste des arguments de la commande ou chaîne de commande complète si shell=True
        :param capture_output: Si True, capture stdout et stderr
        :param shell: Si True, exécute la commande via le shell
        :return: Tuple (stdout, stderr, returncode)
        :raises: subprocess.CalledProcessError si la commande échoue
        """
        try:
            result = subprocess.run(
                command,
                capture_output=capture_output,
                text=True,
                shell=shell,
                check=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            #
            raise(f"Erreur lors de l'éxécution de la commande: {e.stdout, e.stderr, e.returncode}")
        
    def extract_mcannotate_bps(self, mc_annotate_output):
        """
        Extrait les paires de bases de la sortie MC-Annotate.
        
        :param mc_annotate_file: Chemin du fichier contenant le résultat de MC-Annotate
        :return: Liste de tuples représentant les paires de bases
        """
        base_pairs = []

        # Chercher la section 'Base-pairs'
        inside_base_pairs = False
        for line in mc_annotate_output:
            line = line.strip()
            if line.startswith("Base-pairs"):
                inside_base_pairs = True
                
            if inside_base_pairs:
                # Une ligne vide signifie la fin de cette section
                if not line:
                    continue
                # Extraire les paires de bases
                # Exemple de format de ligne: A1-A72 : G-C Ww/Ww pairing antiparallel cis XIX
                if ':' in line:
                    parts = line.split(':')
                    pair_info = parts[0].strip()  # A1-A72
                    details = parts[1].strip()   # G-C Ww/Ww pairing antiparallel cis XIX
                    base_pairs.append((pair_info, details))
        
        index_pairs = []
        for pair, details in base_pairs:
            # Séparer les deux parties de la paire, ex: 'A1-A72'
            nucleotide1, nucleotide2 = pair.split('-')
            # Extraire les indices numériques des nucléotides
            index1 = int(nucleotide1[1:])  # Supprime le 'A' pour obtenir le numéro
            index2 = int(nucleotide2[1:])  # Supprime le 'A' pour obtenir le numéro
            # Ajouter le couple (i, j) à la liste
            index_pairs.append([index1, index2, details])

        #trier les paires
        ignore_terms = [
            "O2P/Bh adjacent_5p pairing",
            "Bs/O2' adjacent_5p pairing",
            "O2'/Ww O2'/Bs inward pairing",
            "O2'/Ww pairing",
            "O2'/Hh adjacent_5p pairing",
            "O2'/Bh pairing",
            "O2P/Bh pairing",
            'Ww/O2P pairing',
            "O2'/Hh Hh/O2P Ww/O2P pairing",
            "Hh/Hh Bh/O2P pairing parallel",
            "O2'/Bs inward pairing", 
            "Ww/Hh pairing antiparallel", 
            "Ww/Ww pairing parallel", 
        ]

        pairs = list()
        mc_annotate_ignore_pairs = list()
        for i, j, details in index_pairs:
            # Vérifier si la description contient un terme à ignorer
            if any(term in details for term in ignore_terms):
                mc_annotate_ignore_pairs.append((i,j))
            else:
                pairs.append((i,j))

        return pairs, mc_annotate_ignore_pairs
    
    def extract_x3dna_bps(self, dssr_output):
        """
        Extrait les paires de bases à partir de la sortie x3DNA DSSR.

        :param dssr_output_content: Contenu du fichier x3DNA DSSR sous forme de chaîne
        :return: Liste des paires sous forme [(i, j, details), ...]
        """
        base_pairs = []

        # Parcourir les lignes de la sortie
        for line in dssr_output.splitlines():
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                try:
                    # Extraire les indices des bases
                    i = int(''.join(filter(str.isdigit, parts[1])))  # Récupérer l'indice numérique de G1, C2, etc.
                    j = int(''.join(filter(str.isdigit, parts[2])))  # Même traitement pour le second résidu
                    #details = f"{parts[3]} {parts[4]}"  # Combiner paire et type d'appariement
                    base_pairs.append((i, j))
                except ValueError:
                    # Ignore les lignes mal formées
                    continue

        return base_pairs
    
    def extract_rnaview_bps(self,rnaview_output):
        """
        Extrait la liste des couples (i, j) représentant les paires de bases à partir
        du fichier de sortie RNAView.
        
        :param rnaview_output_file: Chemin vers le fichier de sortie RNAView
        :return: Liste de tuples représentant les paires de bases
        """

        # Exécuter RNAView
        command = f"module load rnaview/2.0.0 && rnaview --pdb {self.pdb_file}"

        base_pairs = []

        # Lire le fichier RNAView
        for line in rnaview_output:
            # Identifier les lignes contenant "pair-type"
            if line.strip().startswith("pair-type"):
                parts = line.split()
                # Les indices des bases sont dans les 2ème et 3ème colonnes
                i = int(parts[2])
                j = int(parts[3])
                details = str(parts[4])
                base_pairs.append([i, j, details])

        ignore_terms = [
            'HHt', 
            'SWt', 
            'WSt', 
            'HWt', 
            'WWt', 
        ]

        bp = list()
        for i, j, details in base_pairs:
            # Vérifier si la description contient un terme à ignorer
            if any(term in details for term in ignore_terms):
                continue  # Ignorer cette paire
            else:
                bp.append((i,j))
        return bp            
    

    def _add_pair_to_dotb(self, pair):
        (i,j) = pair
        if self.dotb[i - 1] == "." and self.dotb[j - 1] == ".":
                self.dotb[i - 1] = "("
                self.dotb[j - 1] = ")"

    def _has_small_loops(self):
        """
        Retourne une liste vide si la dot-bracket ne contient pas de boucle courte(<= 2 nt)
        Retourne une liste de couple (stsrt, end) des position des courte boucles de la dot bracket

        :param dotb: Chaîne représentant une structure dot-bracket

        """
        stack = []  # Traque les positions des parenthèses ouvrantes
        pairings = {}  # Associe les indices des parenthèses ouvrantes et fermantes

        # Étape 1 : Identifier les appariements dans la dot-bracket
        for i, char in enumerate(self.dotb):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    opening = stack.pop()
                    pairings[opening] = i
                    pairings[i] = opening

        # Étape 2 : Identifier les petites boucles internes
        to_remove = set()  # Indices des paires à supprimer
        for i in pairings:
            if self.dotb[i] == '(':
                j = pairings[i]
                # Vérifier si c'est une boucle interne (contient uniquement des '.')
                if all(self.dotb[k] == '.' for k in range(i + 1, j)):
                    loop_size = j - i - 1
                    if loop_size <= 2:  # Boucle trop courte
                        to_remove.add((i,j)) #(deb, fin) de la boucle
        #On retourner la liste triée
        if to_remove:
            return sorted(to_remove)
        else:
            return list()
    
    def generate_dotb(self):
        """
            Géenère la structure secondaire.
        """

        #1 Gets pairs from mc-annotate
        mc_command = ["/u/sagnioln/stage-E24/tools/MC-Annotate", {self.pdb_file}]
        mc_annotate_out = self.execute_command(mc_command)

        mc_pairs, mc_ignore = self.extract_mcannotate_base_pairs(mc_annotate_out)

        #2 Gets pairs from RNAVIEW
        rna_command = f"module load rnaview/2.0.0 && rnaview --pdb {self.pdb_file}"
        rna_out = self.execute_command(rna_command)
        rnaview_pairs = self.extract_rnaview_base_pairs(rna_out)


        #3 Gets pairs from x3DNA-dssr
        dssr_command = f"module load dssr/2.4.6 && x3dna-dssr -i={self.pdb_file} --idstr=short --pairs-onlys"
        dssr_out = self.execute_command(dssr_command)
        dssr_pairs = self.extract_x3dna_base_pairs(dssr_out)
        
        # On combine les deux liste de paires sans les redondances et en excluant les paires ignorées par mc-annotate
        pairs = mc_pairs
        for elm in rnaview_pairs:
            if (elm not in mc_ignore) and (elm not in pairs):
                pairs.append(elm)
            else:
                mc_ignore.append(elm)

        x3_dna_only = list()
        for (st,ed) in dssr_pairs:
            if ((st,ed) not in pairs) and ((st,ed) not in mc_ignore):
                x3_dna_only = [(idx,jdx) for idx, jdx in dssr_pairs if idx not in set(sum(mc_ignore+pairs,())) 
                            and jdx not in set(sum(mc_ignore + pairs,()))]
                
        primary_dotb = self.turn_2_dotb(pairs)

        #On elève les paires impliquées dans les courtes boucles
        while self.has_small_loops(primary_dotb):
            first = self.has_small_loops(primary_dotb)[0]
            for idx, (i,j) in enumerate(pairs):
                if i == first[0] + 1 :
                    del pairs[idx]
            primary_dotb = self._turn_2_dotb(pairs) 

            #finalement, il faut rajouter les dernieres paires eventuellemenu trouvées par dssr 

        
    def _turn_2_dotb(self, pairs_lst):
        """
        Crée une structure dot-bracket à partir d'une liste de paires [(i,j)].
        
        :param pairs: Liste des coupes des paires
        :param sequence_length: Longueur de la séquence ARN
        :return: Chaîne dot-bracket représentant la structure secondaire
        """

        # Initialiser une liste pour la structure dot-bracket
        dot_bracket = ["."] * len(self.seq)

        for i, j in pairs_lst:
            # Ajouter les paires à la structure dot-bracket
            # Les indices doivent être 0-based pour la liste
            if dot_bracket[i - 1] == "." and dot_bracket[j - 1] == ".":
                dot_bracket[i - 1] = "("
                dot_bracket[j - 1] = ")"
        
        return "".join(dot_bracket)

    def _build_and_improve_stems(self, new_pairs):
        """
        Construit les stems à partir de la structure secondaire (dot-bracket)
        et améliore les stems en ajoutant uniquement les nouvelles paires qui prolongent un stem existant.
        
        :param dotb: str, structure secondaire au format dot-bracket
        :param new_pairs: list of tuples, nouvelles paires à ajouter
        :return: list of stems, où chaque stem est une liste de paires
        """
        # Étape 1 : Construire les stems de base
        initial_stems = self.rna.stems
        
        # Étape 2 : Améliorer les stems avec les nouvelles paires
        for pair in new_pairs:
            i, j = pair
            added = False
            
            for stem in initial_stems:
                # Vérifie si la nouvelle paire prolonge le stem en début ou en fin
                if pair == (stem[-1][0] + 1, stem[-1][1] - 1):  # Prolonge en fin
                    stem.append(pair)
                    added = True
                    break
                elif pair == (stem[0][0] - 1, stem[0][1] + 1):  # Prolonge en début
                    stem.insert(0, pair)
                    added = True
                    break
            
            # Si la paire n'a pas été ajoutée, elle est ignorée
            if not added:
                #print(f"La paire {pair} n'améliore aucun stem et est ignorée.")
                pass
            
        return initial_stems


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exécuter la commande dssr sur un fichier PDB.")
    parser.add_argument("pdb_file", help="Chemin vers le fichier PDB")
    args = parser.parse_args()

    obj = rna_sec_struct_analyser(args.pdb_file)
    
    print(f'La sequence extraite est: {obj.seq}')
    print(f'La structure secondaire est : {obj.dotb}')

