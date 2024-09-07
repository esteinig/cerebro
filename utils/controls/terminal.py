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
    table: Path = typer.Option(
        ..., help="Controls record table"
    ),
    output: Path = typer.Option(
        ..., help="Output summary table"
    ),
):
    """ Summarize controls per sample"""

    df = pandas.read_csv(table, sep="\t", header=0)

    df['control_group'] = df['reference'].apply(map_taxon)
    df = df.groupby(['id', 'control_group'])['alignments'].sum().reset_index()

    df.to_csv(output, sep="\t", index=False)
