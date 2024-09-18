import typer
import warnings
import pandas

from pathlib import Path

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

CONTROLS = {
    "Saccharomyces": ["Saccharomyces"],
    "Pneumocystis": ["Pneumocystis", "NW"],
    "MHV": ["MHV"],
    "Vaccinia": ["Vaccinia"],
    "Truepera": ["Truepera"],
    "Allobacillus": ["Allobacillus"],
    "Imtechella": ["Imtechella"],
    "T4": ["T4-Monash"],
    "PhiX": ["PhiX"]
}

def map_taxon(taxon):
    for key, values in CONTROLS.items():
        if any(taxon.startswith(v) for v in values):
            return key
    return None  # If no match is found

@app.command()
def summarize(
    controls: Path = typer.Option(
        ..., help="Controls record table"
    ),
    quality: Path = typer.Option(
        ..., help="Quality control table"
    ),
    output: Path = typer.Option(
        ..., help="Output summary table"
    ),
):
    """ Summarize controls per sample"""

    df = pandas.read_csv(controls, sep="\t", header=0)
    qc = pandas.read_csv(quality, sep="\t", header=0)

    df['control_group'] = df['reference'].apply(map_taxon)
    df = df.groupby(['id', 'control_group'])['alignments'].sum().reset_index()

    df = df.merge(qc[['id', 'input_reads']], on='id', how='left')
    df['rpm'] = (df['alignments'] / df['input_reads']) * 1e6

    df.to_csv(output, sep="\t", index=False)
