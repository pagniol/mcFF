import networkx as nx
import matplotlib.pyplot as plt

from utils import extends

class CycleGraph:
    def __init__(self, rna_structure):
        self.rna_structure = rna_structure
        self.graph = nx.Graph()
        self.cycles, _ = self.get_ncms()

    def get_ncms(self):
        M = []
        bps = self.rna_structure.bps
        for (i, j) in bps.items():
            closed = False
            type = ""
            k = i
            ip = i
            while not closed:
                k += 1
                while self.rna_structure.dotb[k-1] == '.':
                    k += 1
                type += str(k - ip + 1)
                if bps[i] == k:
                    closed = True
                else:
                    type += "_"
                    k = bps[k]
                    ip = k
            indexes = []
            strands = type.split('_')
            s = ""
            if len(strands) == 1:
                for offset in range(int(strands[0])):
                    indexes.append(offset + i)
                    s += self.rna_structure.seq[indexes[-1] - 1]
            elif len(strands) == 2:
                for offset in range(int(strands[0])):
                    indexes.append(offset + i)
                    s += self.rna_structure.seq[indexes[-1] - 1]
                for offset in range(int(strands[1]), 0, -1):
                    indexes.append(j - offset + 1)
                    s += self.rna_structure.seq[indexes[-1] - 1]
            if indexes:
                M.append((type, s, indexes))
        M.sort(key=lambda x: x[2][0] if x[2] else float('inf'))
        return M, bps

    def build_graph(self):
        
        jnc_nodes = list()
        for elm in self.rna_structure.junctions:
            jnc_nodes.append(elm.get_tuple(self.rna_structure.seq))

        nodes = self.cycles + jnc_nodes
        for idx, (ncm_type, ncm_seq, residues) in enumerate(nodes):
            if ncm_type in ['EXC', 'INC']:
                self.graph.add_node(idx, ncm_type=ncm_type, ncm_seq=ncm_seq, residues=residues, color= '#0000FF')
            else:
                self.graph.add_node(idx, ncm_type=ncm_type, ncm_seq=ncm_seq, residues=residues, color='green')

        # Ajout des arrêtes entre les noeuds et les tiges pour compléter le graphe
        for i in range(len(nodes)):
            for j in range(i + 1, (len(nodes))):
                if set(nodes[i][2]).intersection(set(nodes[j][2])) and (nodes[i][0] and nodes[j][0] not in ['EXC', 'INC'])  :
                    self.graph.add_edge(i, j)
                if nodes[j][0] in ['EXC', 'INC'] :
                    if nodes[j][0] == 'EXC' and set(nodes[i][2]).intersection(extends(nodes[j][2])):
                        self.graph.add_edge(i,j)
                    elif set(nodes[i][2]).intersection(nodes[j][2]):
                        self.graph.add_edge(i,j)
                
                
    def plot_graph(self):
        pos = nx.spring_layout(self.graph)
        plt.figure(figsize=(20, 10))

        node_labels = {}
        for node in self.graph.nodes():
            node_data = self.graph.nodes[node]
            if node_data['ncm_type'] in ['INC', 'EXC']:
                node_labels[node] = f"JUNC:{node_data['ncm_seq']}"
            else:
                node_labels[node] = f"{node_data['ncm_seq']}"

        nx.draw_networkx_nodes(self.graph, pos, node_size=500)
        nx.draw_networkx_edges(self.graph, pos, edge_color="gray")
        nx.draw_networkx_labels(self.graph, pos, labels=node_labels, font_size=8, font_weight="bold")

        plt.title("Graph of RNA Cycles with Shared Nucleotides and Junctions")
        plt.show()

