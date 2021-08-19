
# CAFA Evaluation Example Walkthrough

This document contains instructions for running the evaluation in the `example` directory. 

The example is based on CAFA3 benchmarks. The raw benchmark data has been pre-parsed into a json format and can
be found in the `example/benchmark` directory.

There is a separate file for each ontology/species/evidence-type combination and the file contents have this general form:
```
{
    "benchmark_taxon": "ARATH", 
    "benchmark_taxon_id": 3702, 
    "benchmark_ontology": "cco", 
    "benchmark_ontology_term_count": 3905, 
    "protein_annotations": {
        "T37020000125": ["GO:0005575", "GO:0005622", "GO:0005623", "GO:0005737", "GO:0005739",...]
        "T37020000185": ["GO:0005575", "GO:0005622", "GO:0005623", "GO:0005737", "GO:0005773",...],
        ...
    } 
}
```

The annotation data here has been propagated.

