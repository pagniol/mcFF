from rna_handler import RNA_structure
from build_graph import CycleGraph

def main():
    # Test example with tRNA-PHE
    tRNA_PHE_seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
    tRNA_PHE_dotb = '(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....'

    rna_structure = RNA_structure(tRNA_PHE_seq, tRNA_PHE_dotb)
    cycle_graph = CycleGraph(rna_structure)
    cycle_graph.build_graph()
    cycle_graph.plot_graph()

    print(f'Les paires de base sont:\n {rna_structure.bps}')
    print(f'Les tiges du graphe sont :')
    for stem in rna_structure.stems:
        print(stem)
    
    print(f'Les jonctions trouvées sont: \n {rna_structure.junctions}')

if __name__ == "__main__":
    main()
