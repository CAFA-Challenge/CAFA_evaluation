from operator import itemgetter
from typing import Iterable
from goatools.obo_parser import GODag
import numpy as np
import pandas as pd
from collections import namedtuple
from typing import Optional, NamedTuple


class ConfusionMatrix(NamedTuple):
    """ A namedtuple for storing evaluation metrics. """

    # Wrapper for get_confusion_matrix return data:
    TP: int
    FP: int
    TN: int
    FN: int
    precision: float
    recall: float
    fmeasure: float


cm_metrics = namedtuple(
    "cm_metrics", ["TP", "FP", "TN", "FN", "precision", "recall", "fmeasure"]
)


def threshold_range():
    """ generator for float values from .00 -> 1 inclusive """
    for i in range(0, 101):
        yield i / 100.00


# TODO: this needs a rename to be specific that it operates on a 1-dim array, we
# need a second function that operates on an entire matrix as well
def get_confusion_matrix(
    benchmark: np.array,
    predicted: np.array,
    confidence_threshold: Optional[float] = None,
) -> namedtuple:
    """ Calculates true-positive, false-positive, true-negative, false-negative,
    precision, recall and fmeasure values from two binary np.arrays
    """
    if confidence_threshold is not None:
        predicted[predicted < confidence_threshold] = 0
        # .astype(np.int8)

    # assert len(predicted) == 3905
    # assert len(benchmark) == 3905

    # Trying to figure out the differences between my code and the Matlab code,
    # what would npos and nneg be (see pfp_confmat.m):
    npos = np.sum(predicted != 0)
    nneg = np.sum(predicted == 0)
    # print(f"\nTHRESHOLD: {probability_threshold}, NPOS: {npos} AND NNEG: {nneg}")

    tp = np.sum(np.logical_and(predicted != 0, benchmark == 1))
    fp = np.sum(np.logical_and(predicted != 0, benchmark == 0))
    tn = np.sum(np.logical_and(predicted == 0, benchmark == 0))
    fn = np.sum(np.logical_and(predicted == 0, benchmark == 1))

    # assert (tp + fp + tn + fn) == 3904

    precision = np.nan_to_num(tp / (tp + fp))
    recall = tp / (tp + fn)
    fmeasure = 0  # np.nan_to_num((2 * precision * recall) / (precision + recall))

    # print(f"TRUE POSITIVE IS {int(tp)}")

    ret = ConfusionMatrix(
        int(tp),
        int(fp),
        int(tn),
        int(fn),
        round(precision, 3),
        round(recall, 3),
        round(fmeasure, 3),
    )
    # print(ret)
    # print("=========================")

    return ret


def remove_zero_annotation_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    _dataframe = dataframe.copy()
    nonzero_mask = _dataframe.sum(axis="index") > 0
    cols_of_interest = nonzero_mask.index[nonzero_mask].tolist()
    _dataframe = _dataframe.loc[:, cols_of_interest]
    return _dataframe


def parse_raw_annotation_line(line: str) -> tuple:
    """ Parses a line of raw annotation text.

    Notes
    -------------
    The input string is expected to be one of two formats:
    1. protein_id   go_term_id
    2. protein_id   go_term_id  confidence_level
    """
    line_split = line.rstrip().split()

    if len(line_split) == 3:
        # line is of the second form:
        line_split[2] = float(line_split[2])

    return tuple(line_split)


def parse_annotation_file(annotation_filepath: str, sort: bool = False) -> list:
    """ Parse a txt file where each line is a protein ID - GO ID pair """
    annotations = []

    with open(annotation_filepath, "r") as annotation_handle:
        annotations = set(
            [parse_raw_annotation_line(line) for line in annotation_handle]
        )

    # list is needed for sorting. For a consistent API, always return a list
    annotations = list(annotations)

    if sort:
        annotations.sort(key=itemgetter(0))

    return annotations


def get_annotation_dataframe(
    propagation_map_filepath: str,
    dag: GODag,
    annotation_pairs: Iterable,
    pare_df: bool = False,
) -> pd.DataFrame:
    """Creates and returns a pandas DataFrame representing propagated protein
    annotations.

    the resulting DataFrame has this form:
                       GO:0000123  GO:0000124  GO:0000126  GO:0000137  GO:0000139
    T100900000026           0           0           0           0           0
    T100900000115           0           0           0           0           0
    T100900000116           0           0           0           0           0
    T100900000167           0           0           0           0           0
    T100900000453           0           0           0           0           0
    """
    proteins = sorted(list({protein for protein, _ in annotation_pairs}))

    annotation_df = pd.DataFrame(
        data=0, index=proteins, columns=sorted(list(dag.keys()))
    )
    propagation_map = pd.read_pickle(propagation_map_filepath)
    # propagation_map = pd.read_pickle("./cco_propagation_map.pkl")

    # loop the leaf-only groundtruth, updating the dataframe with propagated
    # annotations (toggle zero vals to 1 where appropriate):
    for protein, node in annotation_pairs:
        # TODO: Consider alt-IDs
        try:
            propagation_df_row = propagation_map.loc[node, :]
            propagation_mask = propagation_df_row == 1
            # column names (GO IDs) that are ancestors of our leaf GO ID:
            propagation_mask_cols = propagation_mask.index[propagation_mask].tolist()
            annotation_df.loc[protein, propagation_mask_cols] = 1
        except KeyError:
            pass

    if pare_df:
        # TODO: For the sake of expediency, temporarily pare the terms down to those
        # TODO: that are actually used by our limited protein dataset:
        term_mask = annotation_df.sum(axis="index") > 0
        cols_of_interest = term_mask.index[term_mask].tolist()
        annotation_df = annotation_df.loc[:, cols_of_interest]

    return annotation_df


def propagate_annotation_dataframe(
    annotation_dataframe: pd.DataFrame,
    propagation_map_filepath: str,
    dag: GODag,
    annotation_pairs: Iterable,
    pare_df: bool = False,
) -> pd.DataFrame:
    """Creates and returns a pandas DataFrame representing propagated protein
    annotations.

    the resulting DataFrame has this form:
                       GO:0000123  GO:0000124  GO:0000126  GO:0000137  GO:0000139
    T100900000026           0           0           0           0           0
    T100900000115           0           0           0           0           0
    T100900000116           0           0           0           0           0
    T100900000167           0           0           0           0           0
    T100900000453           0           0           0           0           0
    """
    proteins = sorted(list({protein for protein, _ in annotation_pairs}))

    annotation_df = pd.DataFrame(
        data=0, index=proteins, columns=sorted(list(dag.keys()))
    )
    propagation_map = pd.read_pickle(propagation_map_filepath)
    # propagation_map = pd.read_pickle("./cco_propagation_map.pkl")

    # loop the leaf-only groundtruth, updating the dataframe with propagated
    # annotations (toggle zero vals to 1 where appropriate):
    for protein, node in annotation_pairs:
        # TODO: Consider alt-IDs
        try:
            propagation_df_row = propagation_map.loc[node, :]
            propagation_mask = propagation_df_row == 1
            # column names (GO IDs) that are ancestors of our leaf GO ID:
            propagation_mask_cols = propagation_mask.index[propagation_mask].tolist()
            annotation_df.loc[protein, propagation_mask_cols] = 1
        except KeyError:
            pass

    if pare_df:
        # TODO: For the sake of expediency, temporarily pare the terms down to those
        # TODO: that are actually used by our limited protein dataset:
        term_mask = annotation_df.sum(axis="index") > 0
        cols_of_interest = term_mask.index[term_mask].tolist()
        annotation_df = annotation_df.loc[:, cols_of_interest]

    return annotation_df


def dataframe_to_binary(dataframe: pd.DataFrame, threshold) -> pd.DataFrame:
    """ dataframe values equal to or greater than threshold are converted
    to 1, all others are converted to zero.

    parameters
    -------------
    dataframe
        a dataframe of protein->GO annotations.
        Dataframe indices should be proteins and columns should be GO terms
        with floating point numbers as cell values
    theshold
        floating point cutoff for converting all dataframe values to 0 or 1

    returns
    ----------
    dataframe
        pandas DataFrame of binary (zero or 1) data
    """
    df = pd.DataFrame(
        data=0,
        index=dataframe.index.tolist(),
        columns=dataframe.columns.tolist(),
        dtype=int,
    )
    df.mask(dataframe >= threshold, other=1, inplace=True)
    return df


def get_prediction_coverage(prediction: pd.DataFrame, benchmark: pd.DataFrame) -> float:
    """ Computes the prediction coverage as percentage of benchmark proteins present
    in the prediction. This is described in figure 2 of
    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1835-8

    This function expects two pandas DataFrames with proteins as indices and
    GO terms as columns.
    """
    benchmark_protein_count = benchmark.shape[0]
    # count proteins from the benchmark that also exist in the prediction and have
    # at least one non-zero prediction value:
    coverage_proteins = [
        protein
        for protein in benchmark.index
        if protein in prediction.index and any(prediction.loc[protein, :] > 0)
    ]
    return len(coverage_proteins) / benchmark_protein_count
