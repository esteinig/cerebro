import pandas

from pathlib import Path


def read_qc_table(file: Path, remove_ntc: bool = False, min_ercc_constructs: int = None, sep: str = "\t") -> pandas.DataFrame:

    df = pandas.read_csv(file, sep=sep, header=0)

    if min_ercc_constructs:
        df = df.loc[df["ercc_constructs"] >= min_ercc_constructs, :]

    if remove_ntc: 
        df = df.loc[~df["id"].str.contains("NTC"), :]

    return df

def get_vircov_ercc_from_results(result_dir: Path, umi: bool = False, after: bool = False) -> pandas.DataFrame:

    """
    Obtain the Vircov alignment summary for ERCC/EDCC results from
    the standard quality control the UMI quality control subworkflows
    """


LAPUTA_MEDIUM = [
    '#14191F',
    '#1D2645',
    '#403369',
    '#5C5992',
    '#AE93BE',
    '#B4DAE5',
    '#F0D77B',
]

YESTERDAY_MEDIUM = [
    '#061A21',
    '#132E41',
    '#26432F',
    '#4D6D93',
    '#6FB382',
    '#DCCA2C',
    '#92BBD9',
]

LAPUTA_MEDIUM.reverse()
YESTERDAY_MEDIUM.reverse()