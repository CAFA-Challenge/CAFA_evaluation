
# CAFA Predictions Evaluation 

Test data are available in FigShare: https://figshare.com/account/projects/123691/articles/16695214

This repository is a collection of code for evaluating CAFA predictions

See also the [walkthrough example](https://github.com/CAFA-Challenge/CAFA_evaluation/blob/master/walkthrough.md)


To use the code, the steps are:
1. Obtain or generate propagation Pandas DataFrame files based on the relevant CAFA challenge obo file. These files are pickled
DataFrames, one per relevant GO ontology namespace that are used for propagating both benchmark and prediction annotation data.
To generate these files yourself, see `generate_propagation_map_dataframe.py`. This file can be called directly from the
   shell with a yaml configuration file. See `obo_parser_config.yml` as an example of the expected configuration.

2. Obtain or generate weighted DAG files (one per GO ontology namespace). These are pickled Pandas DataFrame files containing
the Information Content data for the relevant DAG nodes. To generate your own weighted DAG files, see `generate_information_content_matrices.py`. 
   This file uses the same configuration yaml as step 1.

3. Parse raw benchmark data into json-formatted per-species, per-ontology (CCO, MFO, BPO, etc) files.
This step uses `raw_benchmark_parser.py`. This step requires a yaml configuration file as well. See `parser_config.yml` for 
   an example of the necessary configuration keys and values
   
4. Parse raw prediction data into json-fomatted files on a per-species and per-ontology basis.
This step uses `raw_prediction_parser.py` in conjunction with the same configuration file from step #3. 
   
5. Evaluate the prediction data on a per-species, per-ontology basis against the benchmark data using `evaluate_species_prediction.py`

6. Generate cross-species evaluation metrics using `evaluate_cross_species.py`



