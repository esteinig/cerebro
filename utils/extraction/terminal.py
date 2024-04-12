import typer
import warnings

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM

from matplotlib import pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

@app.command()
def plot_extraction_panel(
    taxa: Path = typer.Option(
        ..., help="Taxon table filtered"
    ),
):
    
    """
    Plot the extraction validation experiment
    """

    


