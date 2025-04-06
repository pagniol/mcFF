def extract_mcannotate_bps(mc_annotate_output_content):
    """
    Extrait les paires de bases de la sortie MC-Annotate.

    :param mc_annotate_output_content: Contenu du fichier MC-Annotate sous forme de chaîne
    :return: Liste de tuples représentant les paires de bases
    """
    base_pairs = []

    # Parcourir les lignes du contenu
    inside_base_pairs = False
    for line in mc_annotate_output_content.splitlines():
        line = line.strip()
        if line.startswith("Base-pairs"):
            inside_base_pairs = True

        if inside_base_pairs:
            # Une ligne vide signifie la fin de cette section
            if not line:
                continue
            # Extraire les paires de bases
            if ':' in line:
                parts = line.split(':')
                pair_info = parts[0].strip()  # A1-A72
                details = parts[1].strip()   # G-C Ww/Ww pairing antiparallel cis XIX
                base_pairs.append((pair_info, details))

    # Convertir les paires en indices numériques
    index_pairs = []
    for pair, details in base_pairs:
        nucleotide1, nucleotide2 = pair.split('-')
        index1 = int(nucleotide1[1:])  # Supprime le 'A' pour obtenir le numéro
        index2 = int(nucleotide2[1:])  # Supprime le 'A' pour obtenir le numéro
        index_pairs.append([index1, index2, details])

    return index_pairs


def extract_rnaview_bps(rnaview_output_content):
    """
    Extrait les paires de bases de la sortie RNAView.

    :param rnaview_output_content: Contenu du fichier RNAView sous forme de chaîne
    :return: Liste de tuples représentant les paires de bases
    """
    base_pairs = []

    # Parcourir les lignes du contenu
    for line in rnaview_output_content.splitlines():
        # Identifier les lignes contenant "pair-type"
        if line.strip().startswith("pair-type"):
            parts = line.split()
            # Les indices des bases sont dans les 2ème et 3ème colonnes
            i = int(parts[2])
            j = int(parts[3])
            details = str(parts[4])
            base_pairs.append([i, j, details])

    return base_pairs


def extract_x3dna_bps(dssr_output_content):
    """
    Extrait les paires de bases à partir de la sortie x3DNA DSSR.

    :param dssr_output_content: Contenu du fichier x3DNA DSSR sous forme de chaîne
    :return: Liste des paires sous forme [(i, j, details), ...]
    """
    base_pairs = []

    # Parcourir les lignes de la sortie
    for line in dssr_output_content.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0].isdigit():
            try:
                # Extraire les indices des bases
                i = int(''.join(filter(str.isdigit, parts[1])))  # Récupérer l'indice numérique de G1, C2, etc.
                j = int(''.join(filter(str.isdigit, parts[2])))  # Même traitement pour le second résidu
                details = f"{parts[3]} {parts[4]}"  # Combiner paire et type d'appariement
                base_pairs.append((i, j, details))
            except ValueError:
                # Ignore les lignes mal formées
                continue

    return base_pairs
