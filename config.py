# The 'list' files use mnemonic short names, while the submitted predictions use
# taxon IDs, so we need a way to translate between the two:
taxonomy_map = {
    9606: 'HUMAN',
    3702: 'ARATH',  # Arabidopsis
    7227: 'DROME',  # Drosophila melanogaster
    10090: 'MOUSE',  # Mus musculus
    10116: 'RAT',  # Rattus norvegicus
    8355: 'XENLA',
    559292: 'YEAST',
    99287: 'SALTY',
    160488: 'PSEPK',
    243273: 'MYCGE',
    243232: 'METJA',
    85962: 'HELPY',
    284812: 'SCHPO',
    44689: 'DICDI',
}

ontologies = ('CCO',) #'MFO', 'BPO')
