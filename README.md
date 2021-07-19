
# CAFA Predictions Evaluation 

This repository is a collection of code for evaluating CAFA predictions

To use the code, the steps are:
1. Parse raw benchmark data into json-formatted per-species, per-ontology (CCO, MFO, BPO, etc) files.
This step uses `raw_benchmark_parser.py`
   
2. Parse raw prediction data into json-fomatted files on a per-species and per-ontology basis.
This step uses `raw_prediction_parser.py`
   
3. Evaluate the prediction data on a per-species, per-ontology basis against the benchmark data using `evaluate_species_prediction.py`

4. Evaluate the prediction data per-ontology across all species using `evaluate_all_species_prediction.py`