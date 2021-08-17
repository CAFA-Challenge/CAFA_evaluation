import pandas as pd


if __name__ == "__main__":

    source_file = "./data/speclist.txt"

    kingdoms = {
        'A': 'archaebacteria',
        'B' : 'bacteria',
        'E': 'eukarya',
        'V': 'viruses and phages',
        'O': 'other',
    }

    parsed_data = {}

    with open(source_file, "r") as source_handle:
        '''
        =======================
        (1) Real organism codes
        =======================

        Code    Taxon    N=Official (scientific) name
                Node     C=Common name
                         S=Synonym
        _____ _ _______  ____________________________________________________________
        AADNV V  648330: N=Aedes albopictus densovirus (isolate Boublik/1994)
                         C=AalDNV
        AAV2  V   10804: N=Adeno-associated virus 2
        
        '''
        lines = "".join(source_handle.readlines())

        start_match_str = '''=======================
(1) Real organism codes
======================='''

        end_match_str = '''-----------------------------------------------------------------------
Copyrighted by the UniProt Consortium, see '''

        end_match_str = '''=======================================================================
(2) "Virtual" codes that regroup organisms at a certain taxonomic level
'''



        start_index = lines.find(start_match_str) + len(start_match_str) + 1
        end_index = lines.find(end_match_str)
        lines = [l.strip() for l in lines[start_index: end_index].strip().split("\n")]
        lines = lines[4:]

        current_organism = None
        last_organism = None
        name = None
        common_name = None
        synonoym = None

        for line in lines[:20]:

            if line[1] != '=':
                mnemonic, kingdom, ncbi_id = line.split(":")[0].split()
                current_organism = mnemonic
                name = line.split(":")[-1].split("=")[-1].strip()
                common_name = None
                synonoym = None

                if last_organism is not None:
                    last_organism = mnemonic

                print(mnemonic, kingdoms.get(kingdom, kingdom), ncbi_id, common_name, synonoym)

            elif current_organism is not None:
                _split = line.split()

                if _split[0] == 'C':
                    common_name = _split[-1]
                elif _split[0] == 'S':
                    synonoym = _split[-1]




        #for line in lines:

