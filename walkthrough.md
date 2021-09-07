# Example Usage

## Create foundation data
1. Obtain an obo file specific to the CAFA challenge of interest. An obo file
for CAFA 3 is available at https://iastate.box.com/v/cafa-evaluation-example-data

2. Obtain "raw" benchmark data for the relevant CAFA challenge. For CAFA3, this data is available in the `raw_benchmark` 
folder at https://iastate.box.com/v/cafa-evaluation-example-data

3. Create a yaml configuration file used for parsing the obo file as well as the benchmark files. The configuration file should
have this general form:
```yaml
obo_filepath: ./data/go_cafa3.obo
dag_directory: ./data_v3/dag_ia/
propagation_map_directory: ./data_v3/propagation/

ontologies:
  - short_name: CCO
    long_name: cellular_component
    raw_benchmark_filepath: ./data/benchmark/raw/leafonly_CCO.txt
    propagation_df_filepath: ./data_v3/propagation/propagation_map_df_CCO.pkl

  - short_name: MFO
    long_name: molecular_function
    raw_benchmark_filepath: ./data/benchmark/raw/leafonly_MFO.txt
    propagation_df_filepath: ./data_v3/propagation/propagation_map_df_MFO.pkl

  - short_name: BPO
    long_name: biological_process
    raw_benchmark_filepath: ./data/benchmark/raw/leafonly_BPO.txt
    propagation_df_filepath: ./data_v3/propagation/propagation_map_df_BPO.pkl
```

4. Use `generate_propagation_map_dataframe.py` to generate propagation data. Given a yaml file named `my_obo_parser.yml`,
call the python file like so: `python generate_propagation_map_dataframe.py my_obo_parser.yml`. This will create pickled Pandas
DataFrame files at the locations specified by `propagation_df_filepath:` keys in the yaml file.

5. Generate Information Content (IC) data by running `python generate_information_content_matrics.py ./my_obo_parser.yml`. This will 
generate one file per ontology in the location specified by the `dag_directory` key in the yaml file. The generated files
are, again, pickled Pandas DataFrame objects.


## Create intermediate benchmark and prediction data structures
1. Create an second yaml configuration file with this general form:
```yaml
prediction_filepath: data/predictions
benchmark_filepath: data/parsed_benchmark
predictor_group_name: ExamplePredictionsLab1
model_id: 1
ontologies: [CCO, BPO, MFO]

benchmark_json_directory: ./data/parsed_benchmark
dag_ic_directory: ./data_v3/dag_ia
raw_predictions_directory: /media/airport/cafa3_submissions/ExamplePredictionsLab1
predictions_json_directory: ./data_v3/parsed_predictions/ExamplePredictionsLab1/json
predictions_dataframes_directory: ./data_v3/parsed_predictions/ExamplePredictionsLab1/dataframes

# delimiter needs to be quoted:
prediction_file_delimiter: "\t"
prediction_file_header_line_count: 4
knowledge_type: 1
propagation_df_directory: ./data/propagation

raw_benchmark_path:  ./data/benchmark/raw
propagation_map_directory: ./data_v3/propagation/
```

1. Parse raw benchmark data into json-formatted per-species, per-ontology (CCO, MFO, BPO, etc) files.
This step uses `raw_benchmark_parser.py` and a yaml configuration file as well. See the sample yaml immediately above
for an example of the necessary configuration keys and values. The syntax for running this step is simply `python raw_benchmark_parser.py <my_config.yml>`.

    Once this has run, the directory that corresponds with the `benchmark_json_directory` key in the yaml configuration should be populated
  with several json files.

2. Now we need to parse the "raw" prediction data into intermediate json files, similar to what we did with the benchmark data in the previous step. 

    Run `python raw_prediction_parser.py <my_config.yml>`. 
    
## Compute Evaluation Metrics

1. Run `evaluate_species_prediction.py`. This creates a pickled Pandas DataFrame for every prediction json file that has a corresponding (matching species and ontology) benchmark json file. These files will be created in a directory corresponding with the `predictions_dataframes_directory` key in the configuration file. The syntax for calling this file is `python evaluate_species_prediction <my_config.yml>`.

    The created DataFrames have the form where the keys are a protein ID and threshold pairing:
    ```
    |                         | ontology   |   taxon_id | taxon   |   tp_ia |   fp_ia |   fn_ia |   benchmark_ia |   ru |      mi |   weighted_precision |   weighted_recall |   tp |   fp |   fn |      tn |   precision |   recall |
    |:------------------------|:-----------|-----------:|:--------|--------:|--------:|--------:|---------------:|-----:|--------:|---------------------:|------------------:|-----:|-----:|-----:|--------:|------------:|---------:|
    | ('T100900000782', 0.01) | cco        |      10090 | MOUSE   |  23.523 | 23.4959 |       0 |         23.523 |    0 | 23.4959 |             0.500288 |                 1 |   54 |   25 |    0 | 3824.32 |    0.683544 |        1 |
    | ('T100900000782', 0.02) | cco        |      10090 | MOUSE   |  23.523 | 12.9745 |       0 |         23.523 |    0 | 12.9745 |             0.644511 |                 1 |   54 |   14 |    0 | 3835.21 |    0.794118 |        1 |
    | ('T100900000782', 0.03) | cco        |      10090 | MOUSE   |  23.523 | 10.8047 |       0 |         23.523 |    0 | 10.8047 |             0.685249 |                 1 |   54 |    9 |    0 | 3840.14 |    0.857143 |        1 |
    | ('T100900000782', 0.04) | cco        |      10090 | MOUSE   |  23.523 | 10.8047 |       0 |         23.523 |    0 | 10.8047 |             0.685249 |                 1 |   54 |    9 |    0 | 3840.14 |    0.857143 |        1 |
    | ('T100900000782', 0.05) | cco        |      10090 | MOUSE   |  23.523 | 10.8047 |       0 |         23.523 |    0 | 10.8047 |             0.685249 |                 1 |   54 |    9 |    0 | 3840.14 |    0.857143 |        1 |
    ```

2. Once a DataFrame for each species/ontology pair has been generated, `evaluate_cross_species.py` can be used to compute per-ontology, cross-species metrics. This file uses the same yaml configuration file as `evaluate_species_prediction.py` in the previous step. 