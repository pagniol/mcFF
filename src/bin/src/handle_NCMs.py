import networkx as nx
import matplotlib.pyplot as plt

def getBPs(dotb):
    fivePr = '('
    threeP = ')'
    single = '.'
    S = []      # stack of parentheses
    SS = dict() # resulting bps

    j = 0
    for c in dotb:
        j += 1
        if c in fivePr:
            S.append((j, c))
        elif c in threeP:
            if S == []:
                print('Unbalanced dot bracket [', dotb, ']')
                exit(1)
            else:
                i, _ = S.pop()
                SS[i] = j
        else:
            if not c in single:
                print('Illegal symbol in dot bracket [', dotb, ']')
                exit(1)

    if S == []:
        return SS
    else:
        print('Unbalanced dotb [', dotb, ']')
        return None

def getNCMs(seq, dotb):
    M = []
    bps = getBPs(dotb)
    if isinstance(bps, str):
        return bps
    else:
        for (i, j) in bps.items():
            closed = False
            type = ""
            k = i
            ip = i
            while not closed:
                k += 1
                while dotb[k-1] == '.':
                    k += 1
                type += str(k - ip + 1)
                if bps[i] == k:
                    closed = True
                else:
                    type += "_"
                    k = bps[k]
                    ip = k
            try:
                indexes = []
                strands = type.split('_')
                s = ""
                if len(strands) == 1:
                    for offset in range(int(strands[0])):
                        indexes.append(offset + i)
                        s += seq[indexes[-1] - 1]
                elif len(strands) == 2:
                    for offset in range(int(strands[0])):
                        indexes.append(offset + i)
                        s += seq[indexes[-1] - 1]
                    for offset in range(int(strands[1]), 0, -1):
                        indexes.append(j - offset + 1)
                        s += seq[indexes[-1] - 1]
                if indexes:
                    M.append((type, s, indexes))
            except KeyError:
                M.append((type, None, []))
        M.sort(key=lambda x: x[2][0] if x[2] else float('inf'))
        return M, bps
    
def count_nt_NCM(ncms):
    return sum(len(ncm[2]) for ncm in ncms)

def get_junctions(stems):
    """
    Construit les jonctions entre les tiges.
    Chaque jonction est représentée par un tuple avec :
    - La position de début et de fin du single-strand entre les tiges
    - Les coordonnées du premier résidu de chaque tige reliée par cette jonction.
    """
    junctions = []

    if not all(isinstance(stem, list) for stem in stems):
        raise ValueError("Error: stems list contains non-list elements")

    stems = sorted(stems, key=lambda stem: stem[0][0])  # On trie par la position de début de chaque tige

    for i in range(len(stems) - 1):
        current_stem = stems[i]
        next_stem = stems[i + 1]

        current_end = current_stem[-1][0]  # Position finale de la tige courante
        next_start = next_stem[0][0]       # Position initiale de la tige suivante

        if current_end + 1 < next_start:
            junctions.append((current_end + 1, next_start - 1, 
                              (current_stem[0][0], current_stem[0][1]), 
                              (next_stem[0][0], next_stem[0][1])))
        else:
            junctions.append((current_end, current_end, 
                              (current_stem[0][0], current_stem[0][1]), 
                              (next_stem[0][0], next_stem[0][1])))
    
    return junctions

def get_stems(dotb):
    """
    Construit les stems à partir des paires de bases appariées.
    """
    bps = getBPs(dotb)
    if not bps:
        return []

    stems = []
    current_stem = []
    
    # Trier les paires de bases dans l'ordre croissant des positions
    sorted_bps = sorted(bps.items())

    for i, j in sorted_bps:
        if not current_stem:
            current_stem = [(i, j)]
        else:
            # Vérifie la continuité des paires pour regrouper dans le même stem
            last_i, last_j = current_stem[-1]
            if i == last_i + 1 and j == last_j - 1:
                current_stem.append((i, j))
            else:
                stems.append(current_stem)
                current_stem = [(i, j)]
    
    # Ajouter le dernier stem si non vide
    if current_stem:
        stems.append(current_stem)
    
    return stems

# Function to build a cycle graph based on shared nucleotides
def build_cycle_graph(seq, dotb):
    ncms, bps = getNCMs(seq, dotb)
    
    print(f'Les cycles : {len(ncms)}')
    for ncm in ncms:
        print(ncm)

    stems = get_stems(dotb)

    print(f'Nombre de nt impliqués dans les cycles : {count_nt_NCM(ncms)}')

    print(f'Les tiges sont : \n {stems}')

    junctions = get_junctions(stems)

    print(f'Les jonctions sont: \n {junctions}')
    G = nx.Graph()
    
    # Ajout des cycles comme nœuds
    for idx, (ncm_type, ncm_seq, residues) in enumerate(ncms):
        G.add_node(idx, ncm_type=ncm_type, ncm_seq=ncm_seq, residues=residues)

    # Ajout des tiges au graphe
    for i in range(len(ncms)):
        for j in range(i + 1, len(ncms)):
            if set(ncms[i][2]).intersection(ncms[j][2]):
                G.add_edge(i, j)

    # Ajout des jonctions comme nœuds avec connexions aux stems adjacents
    for junction_idx, (start, end, stem1, stem2) in enumerate(junctions, start=len(ncms)):
        G.add_node(junction_idx, ncm_type="junction", residues=(start, end))
        
        # Connexion de la jonction aux stems adjacents
        stem1_index = [idx for idx, ncm in enumerate(ncms) if ncm[2] and ncm[2][0] == stem1[0]]
        stem2_index = [idx for idx, ncm in enumerate(ncms) if ncm[2] and ncm[2][0] == stem2[0]]
        
        if stem1_index:
            G.add_edge(junction_idx, stem1_index[0])
        if stem2_index:
            G.add_edge(junction_idx, stem2_index[0])
    
    return G

def plot_cycle_graph(G):
    pos = nx.spring_layout(G)
    plt.figure(figsize=(50, 10))
    
    node_labels = {}
    for node in G.nodes():
        node_data = G.nodes[node]
        if node_data['ncm_type'] == 'junction':
            node_labels[node] = f"Junction {node_data['residues']}"
        else:
            node_labels[node] = f"{node_data['ncm_seq']}"

    nx.draw_networkx_nodes(G, pos, node_color="skyblue", node_size=500)
    nx.draw_networkx_edges(G, pos, edge_color="gray")
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8, font_weight="bold")

    plt.title("Graph of RNA Cycles with Shared Nucleotides and Junctions")
    plt.show()

# Test example with tRNA-PHE
tRNA_PHE_seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
tRNA_PHE_dotb = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'

G = build_cycle_graph(tRNA_PHE_seq, tRNA_PHE_dotb)
print(G.nodes)
print(G.edges)
plot_cycle_graph(G)
