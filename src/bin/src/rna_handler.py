from handle_junctions import Junction

from handle_junctions import find_all_jonctions

class RNA_structure:
    def __init__(self, seq, dotb):
        self.seq = seq
        self.dotb = dotb
        self.bps = self._get_base_pairs()
        self.junctions = self.build_junctions()
        self.stems = self.get_stems()
    
    def _get_base_pairs(self):
        fivePr = '('
        threeP = ')'
        single = '.'
        S = []      # stack of parentheses
        SS = dict() # resulting bps

        j = 0
        for c in self.dotb:
            j += 1
            if c in fivePr:
                S.append((j, c))
            elif c in threeP:
                if not S:
                    raise ValueError(f"Unbalanced dot bracket [{self.dotb}]")
                else:
                    i, _ = S.pop()
                    SS[i] = j
            elif c not in single:
                raise ValueError(f"Illegal symbol in dot bracket [{self.dotb}]")

        if S:
            raise ValueError(f"Unbalanced dotb [{self.dotb}]")
        
        return SS
    
    def _find_paired_base(self, idx):
        """
        Trouve l'indice de la base appariée à la base à l'indice i dans une structure secondaire d'ARN au format dot-bracket.

        Args:
            dot_bracket (str): La structure secondaire d'ARN au format dot-bracket.
            i (int): L'indice de la base pour laquelle on cherche la base appariée.

        Returns:
            int: L'indice de la base appariée, ou -1 si aucune base n'est appariée.
        """
        stack = []
        dot_bracket = self.dotb
        for j, char in enumerate(dot_bracket):
            if char == '(':
                stack.append(j)
            elif char == ')':
                if not stack:
                    raise ValueError('Unbalanced dotb [', dot_bracket, ']')  # Pas '(' correspondante
                paired_index = stack.pop()
                if paired_index == idx:
                    return j +1 # on ajoute un car les pbs sont comptés à partir de 1
                elif j == idx:
                    return paired_index +1 # on ajoute un car les pbs sont comptés à partir de 1
            
        if stack:
            raise ValueError('Unbalanced dotb [', dot_bracket, ']')  # Pas de ')' correspondante

        return -1
            
    def get_stems(self):
        """
        Construit les stems à partir des paires de bases appariées.
        """
        stems = []
        current_stem = []
        
        sorted_bps = sorted(self.bps.items())
        for i, j in sorted_bps:
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
        
        return stems
    
    def build_junctions(self):
        '''
        Cette fonction est à modifier le plus ôt possible.
        Il n'est question de modifier la première ligne par la fonction self.find_junctions() de cette mm class
        Actuellement, elle fonctionne correctement mais à cause du fait que les paires de base qui lui 
        sont fournies ne sont pas correct, et c'est le bordel à déboguer.
        '''
        jc = self.find_junctions()
        junc = list()
        for _, (start, end, type) in enumerate(jc):
            #cas des jonctions explicites. Exple: )...( ou )(
            if start < end and self.dotb[start] == ')' and self.dotb[end] == '(':
                if end == start +1:
                    junc.append(Junction(start,end, type))
                else:
                    junc.append(Junction(start,end, type))
            elif end == start + 1 and self.dotb[start] in '()' and self.dotb[start] == self.dotb[end] :
                junc.append(Junction(start,end, type))
            elif start < end and (start or end in '(') :
                junc.append(Junction(start,end, type))
            else:
                raise(f'Not identidied Junction {(start, end)}')

        return junc

    def find_junctions(self):
        """
        Fonction principale qui regroupe les jonctions trouvées par les fonctions précédentes
        et retourne un dictionnaire des différents types de jonctions.

        Args:
            dot_bracket (str): Structure secondaire d'ARN au format dot-bracket.

        Returns:
            dict: Un dictionnaire contenant les jonctions par type.
        """
        #Trouve les positions où une tige se termine sans qu'une autre ne commence immédiatement après.
        dot_bracket = self.dotb
        stems = self.get_stems()
        
        stem_endings = []
        for stem in stems:
            # Vérifier uniquement la première paire de la tige
            first_pair = stem[0]
            paired_base = self._find_paired_base(first_pair[0])
            
            # Vérifier si le caractère après la base complémentaire est ')' ou un '.' suivi de ')'
            if paired_base + 1 < len(dot_bracket):
                next_char = dot_bracket[paired_base + 1]
                if next_char == ')' or (next_char == '.' and paired_base + 2 < len(dot_bracket) and dot_bracket[paired_base + 2] == ')'):
                    stem_endings.append((paired_base + 1, paired_base + 2, 'INC')) #

        dot_sections = []
        i = 0
        while i < len(dot_bracket):
            if dot_bracket[i] == '.':
                start = i
                # Recherche d'une séquence de '.' délimitée par le même caractère
                while i < len(dot_bracket) and dot_bracket[i] == '.':
                    i += 1
                end = i
                # Vérifier si la séquence est bien délimitée par des caractères identiques
                if start > 0 and end < len(dot_bracket) and dot_bracket[start - 1] == dot_bracket[end]:
                    if dot_bracket[end] == '(':
                        dot_sections.append((start + 1, end,'EXC')) # 
                    else:
                        dot_sections.append((start +1, end,'EXC')) # 

            else:
                i += 1

        inversion_regions = []
        # On parcourt la structure et on cherche les inversions de type ')(' et celles avec des points.
        i = 0
        while i < len(dot_bracket) - 1:
            if dot_bracket[i] == ')' and dot_bracket[i + 1] == '(':
                inversion_regions.append((i +1, i + 2,'INC')) #
                i += 2  # Passer à l'index après la paire inversée
            elif dot_bracket[i] == ')' and dot_bracket[i + 1] == '.':
                # Vérifier si une inversion existe avec des points
                j = i + 1
                while j < len(dot_bracket) and dot_bracket[j] == '.':
                    j += 1
                if j < len(dot_bracket) and dot_bracket[j] == '(':
                    inversion_regions.append((i+2, j,'EXC'))  #
                    i = j  # Passer à la position après la parenthèse ouvrante
                else:
                    i += 1
            else:
                i += 1
                
        return inversion_regions + stem_endings + dot_sections


def main(dotb, seq):
    rna = RNA_structure(seq, dotb)
    print(f'Les paires de bases sont : \n {rna.bps}')
    print(f'Les jonctions sont :\n {rna.junctions}')

if __name__ == '__main__':
    # Test example with tRNA-PHE
    tRNA_PHE_seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
    tRNA_PHE_dotb = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'
    main(tRNA_PHE_dotb, tRNA_PHE_seq)
