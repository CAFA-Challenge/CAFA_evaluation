from pathlib import Path
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


if __name__ == "__main__":
    import sys
    import yaml

    config_filepath = sys.argv[1]

    with open(config_filepath, "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        output_directory_str = config.get("propagation_map_directory")
        output_directory = Path(output_directory_str)
        output_directory.mkdir(exist_ok=True, parents=True)

        obo_filepath = config.get('obo_filepath')
        ontologies = config.get('ontologies')

        optional_attrs = ["relationship", "replaced_by", "consider"]
        optional_relationships = {'part_of', }
        dag = GODag(obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None)

        for namespace in ontologies:
            namespace_long = namespace.get("long_name")
            namespace_short = namespace.get("short_name")
            print(f"PROCESSING {namespace_long}/{namespace_short}")
            output_filepath = output_directory / f"propagation_map_df_{namespace_short}.pkl"

            # by using the item_id field and not 'term', we translate term alt-ids to canonical ids:
            valid_namespace_terms = {dag[term].item_id for term in dag if dag[term].namespace == namespace_long}
            valid_namespace_terms = sorted(list(valid_namespace_terms))

            prop_map = pd.DataFrame(data=0, index=valid_namespace_terms, columns=valid_namespace_terms)
            term_subdag = GoSubDag(valid_namespace_terms, dag, relationships=optional_relationships, prt=None)

            for term in valid_namespace_terms:
                ancestors = term_subdag.rcntobj.go2parents[term]
                ancestors.add(term)
                prop_map.loc[term, ancestors] = 1

            print(f"\tWRITING {output_filepath}\n")
            prop_map.to_pickle(output_filepath)