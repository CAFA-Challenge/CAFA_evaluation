# The 'list' files use non-standard short names, while the submitted predictions use official
# taxon IDs, so we need a way to translate between the two:
taxonomy_map = {
    9606: 'HUMAN',
    3702: 'ARATH',  # Arabidopsis
    7227: 'DROME',  # Drosophila melanogaster
    10090: 'MOUSE',  # Mus musculus
    10116: 'RAT',  # Rattus norvegicus
}

ontologies = ('CCO',) #'MFO', 'BPO')
