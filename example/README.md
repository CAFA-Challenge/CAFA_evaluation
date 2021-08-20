
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

    This step produces one file per species and ontology pair. Additionally, the predicted annotations are propagated.

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
            ...
        }, 
        "T37020000185": {
            "GO:0005575": 0.97,
            "GO:0005576": 0.01, 
            "GO:0005618": 0.01, 
            "GO:0005622": 0.97, 
            "GO:0005623": 0.97,
            ...
        },
        "T37020000188": {
            "GO:0000325": 0.02, 
            "GO:0005575": 0.22, 
            "GO:0005576": 0.01,
            ...
        }
    }
    ```

2. 
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
|                        | ontology   |   taxon_id | taxon   |   tp_ia |   fp_ia |   fn_ia |   benchmark_ia |   weighted_precision |   weighted_recall |   tp |   fp |   fn |      tn |   precision |   recall |
|========================+============+============+=========+=========+=========+=========+================+======================+===================+======+======+======+=========+=============+==========+
| ('T96060000375', 0.01) | cco        |       9606 | HUMAN   | 1.95999 | 54.1254 | 0       |        1.95999 |            0.0349466 |                 1 |   10 |   62 |    0 | 3831.86 |    0.138889 |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.02) | cco        |       9606 | HUMAN   | 1.95999 | 32.3119 | 0       |        1.95999 |            0.0571894 |                 1 |   10 |   35 |    0 | 3858.78 |    0.222222 |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.03) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.04) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.05) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.06) | cco        |       9606 | HUMAN   | 1.95999 |  0      | 0       |        1.95999 |            1         |                 1 |   10 |    0 |    0 | 3893    |    1        |        1 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.07) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.08) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.09) | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+
| ('T96060000375', 0.1)  | cco        |       9606 | HUMAN   | 0       |  0      | 1.95999 |        1.95999 |            0         |                 0 |    0 |    0 |   10 | 3895    |    0        |        0 |
+------------------------+------------+------------+---------+---------+---------+---------+----------------+----------------------+-------------------+------+------+------+---------+-------------+----------+




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