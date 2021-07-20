import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


if __name__ == "__main__":
    NAMESPACE = "cellular_component"
    NAMESPACE_SHORT = "CCO"

    NAMESPACE = "biological_process"
    NAMESPACE_SHORT = "BPO"

    NAMESPACE = "molecular_function"
    NAMESPACE_SHORT = "MFO"

    obo_filepath = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    output_filepath = f"./data/propagation/propagation_map_df_{NAMESPACE_SHORT}.pkl"

    optional_attrs = ["relationship", "replaced_by", "consider"]
    optional_relationships = {'part_of', }

    dag = GODag(obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None)

    # by using the item_id field and not 'term', we translate term alt-ids to canonical ids:
    valid_namespace_terms = {dag[term].item_id for term in dag if dag[term].namespace == NAMESPACE}
    valid_namespace_terms = sorted(list(valid_namespace_terms))

    prop_map = pd.DataFrame(data=0, index=valid_namespace_terms, columns=valid_namespace_terms)
    term_subdag = GoSubDag(valid_namespace_terms, dag, relationships=optional_relationships, prt=None)

    for term in valid_namespace_terms:
        ancestors = term_subdag.rcntobj.go2parents[term]
        ancestors.add(term)
        prop_map.loc[term, ancestors] = 1

    prop_map.to_pickle(output_filepath)