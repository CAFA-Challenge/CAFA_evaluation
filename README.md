
# CAFA Predictions Evaluation 

This repository is a collection of code for evaluating CAFA predictions

To use the code, the steps are:
1. Obtain or generate propagation Pandas DataFrame files based on the relevant CAFA challenge obo file. These files are pickled
DataFrames, one per relevant GO ontology namespace that are used for propagating both benchmark and prediction annotation data.
To generate these files yourself, see `generate_propagation_map_dataframe.py`.

2. Obtain or generate weighted DAG files (one per GO ontology namespace). These are pickled Pandas DataFrame files containing
the Information Content data for the relevant DAG nodes. To generate your own weighted DAG files, see `generate_information_content_matrices.py`

3. Parse raw benchmark data into json-formatted per-species, per-ontology (CCO, MFO, BPO, etc) files.
This step uses `raw_benchmark_parser.py`
   
4. Parse raw prediction data into json-fomatted files on a per-species and per-ontology basis.
This step uses `raw_prediction_parser.py`
   
5. Evaluate the prediction data on a per-species, per-ontology basis against the benchmark data using `evaluate_species_prediction.py`

6. Evaluate the prediction data per-ontology across all species using `evaluate_all_species_prediction.py`

7. Generate cross-species evaluation metrics using `evaluate_cross_species.py`