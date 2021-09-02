
# CAFA Evaluation Example Walkthrough

This document contains instructions for running the evaluation in the `example` directory. 

The example is based on CAFA3 benchmarks. The raw benchmark data has been pre-parsed into a json format and can
be found in the `example/benchmark` directory.

Alternatively, these benchmark files can be generated with `raw_benchmark_parser.py`

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


There are also pandas DataFrames (one per ontology namespace) that are used for propagating ontology terms saved as pkl files.
Unfortunately, these files are too large for github, so those files are stored elsewhere: `https://iastate.box.com/v/cafa-evaluation-example-data/propagatation`

Next, there are pandas DataFrames (one per ontology namespace) representing the information content of the benchmark
annotations. That data is in the `example/dag_ia` directory of this git repository.

##Steps

1. Parse "raw" CAFA prediction files into intermediate json-formatted files. 

    This step produces one file per species/ontology pair. Additionally, the predicted "leaf" annotations are propagated.

    The output json files contain nested dictionaries mapping proteins to terms and thresholds with this form:
    ```json
    {
        "T37020000125": {
            "GO:0005575": 0.04, 
            "GO:0005576": 0.01, 
            "GO:0005618": 0.01, 
            "GO:0005622": 0.02, 
            "GO:0005623": 0.02, 
            "GO:0005634": 0.01, 
            "GO:0005654": 0.01,
        }, 
        "T37020000185": {
            "GO:0005575": 0.97,
            "GO:0005576": 0.01, 
            "GO:0005618": 0.01, 
            "GO:0005622": 0.97, 
            "GO:0005623": 0.97,
        },
        "T37020000188": {
            "GO:0000325": 0.02, 
            "GO:0005575": 0.22, 
            "GO:0005576": 0.01,
        }
    }
    ```

2. Compute evaluation metrics from the json data generated in step #1 by running `evaluate_species_prediction.py`

    This step creates pandas DataFrames containing per-protein evaluation metrics.
    
    There is one DataFrame per taxon/ontology pair.
    
    The DataFrames have this form:
 
    |                        | ontology   |   taxon_id | taxon   |   tp_ia |   fp_ia |   fn_ia |   benchmark_ia |   weighted_precision |   weighted_recall |   tp |   fp |   fn |      tn |   precision |   recall |
    |:-----------------------|:-----------|-----------:|:--------|--------:|--------:|--------:|---------------:|---------------------:|------------------:|-----:|-----:|-----:|--------:|------------:|---------:|
    | ('T96060000375', 0.01) | cco        |       9606 | HUMAN   | 1.95999 | 54.1254 | 0       |        1.95999 |            0.0349466 |                 1 |   10 |   62 |    0 | 3831.86 |    0.138889 |        1 |
    | ('T96060000375', 0.02) | cco        |       9606 | HUMAN   | 1.95999 | 32.3119 | 0       |        1.95999 |            0.0571894 |                 1 |   10 |   35 |    0 | 3858.78 |    0.222222 |        1 |
    | ('T96060000375', 0.03) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
    | ('T96060000375', 0.04) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
    | ('T96060000375', 0.05) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
    | ('T96060000375', 0.06) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
    | ('T96060000375', 0.07) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
    | ('T96060000375', 0.08) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
    | ('T96060000375', 0.09) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
    | ('T96060000375', 0.1)  | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
    
    
3. Finally, the data from the previous step can be combined to provide general metrics on a per-ontology basis using `evaluate_cross_species.py`

    |   threshold | ontology   | average_precision |   average_recall |   average_weighted_precision |   average_weighted_recall |
    |------------:|:-----------|--------------------:|-----------------:|---------------------------:|--------------------------:|
    |        0.01 | cco        |            0.273477  |        0.209426  |                 0.168031  |                0.118232   |
    |        0.02 | cco        |            0.27435   |        0.199523  |                 0.168737  |                0.107476   |
    |        0.03 | cco        |            0.273529  |        0.191298  |                 0.168631  |                0.0996621  |
    |        0.04 | cco        |            0.270859  |        0.185562  |                 0.167624  |                0.0951871  |
    |        0.05 | cco        |            0.267388  |        0.18051   |                 0.16605   |                0.0916171  |
    |        0.06 | cco        |            0.263533  |        0.175753  |                 0.164141  |                0.0883728  |
    |        0.07 | cco        |            0.259436  |        0.171403  |                 0.161765  |                0.0853578  |
    |        0.08 | cco        |            0.255555  |        0.167408  |                 0.159595  |                0.082586   |
    |        0.09 | cco        |            0.251594  |        0.163661  |                 0.157325  |                0.080092   |
    |        0.1  | cco        |            0.247779  |        0.160089  |                 0.155134  |                0.077756   |
    |        0.11 | cco        |            0.243962  |        0.156662  |                 0.152863  |                0.0755321  |
    |        0.12 | cco        |            0.240168  |        0.153611  |                 0.150614  |                0.0736488  |
    |        0.13 | cco        |            0.23655   |        0.150693  |                 0.148466  |                0.0718578  |
    |        0.14 | cco        |            0.233022  |        0.147986  |                 0.146355  |                0.0701898  |
    |        0.15 | cco        |            0.229759  |        0.145576  |                 0.14437   |                0.0687319  |
    |        0.16 | cco        |            0.226621  |        0.143347  |                 0.142366  |                0.0674593  |
    |        0.17 | cco        |            0.223707  |        0.141273  |                 0.140468  |                0.0662037  |
