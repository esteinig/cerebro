import typer
import pandas
import warnings

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM


from matplotlib import pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

@app.command()
def plot_group_metric(
    taxa: Path = typer.Option(
        ..., help="Taxa overview summary table"
    ),
    group: str = typer.Option(
        "sample_id", help="Grouping variable to plot by"
    ),
    metric: str = typer.Option(
        "rpm", help="Metric to plot"
    ),
    top: int = typer.Option(
        None, help="Extract the `top` taxa by `top-metric`"
    ),
    top_metric: str = typer.Option(
        None, help="Extract the top taxa by `top-metric` (total_rpm, kmer-rpm, alignment_rpm, contigs)"
    ),
):
    
    """
    Plot taxa metrics across a grouping variable
    """

    df = pandas.read_csv(taxa, sep=",", header=0)

    fig1, axes1 = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
            
    ax1 = axes1
    
    df_taxa = df[df["name"].isin(["Entamoeba dispar", "Escherichia coli"])]
    sns.stripplot(x=group, y=metric, hue="name", data=df_taxa, palette=YESTERDAY_MEDIUM, ax=ax1, size=8, alpha=0.7)
    ax1.set_xlabel(f"\n{group}")
    ax1.set_ylabel(f"{metric}\n")
    ax1.tick_params(axis='x', labelsize=8, rotation=45)
    legend = ax1.get_legend()
    legend.set_title("Taxa") 
    sns.despine()
    ax1.grid(False)




    # fig1_suptitle = "Quality control"
    # fig1.suptitle(fig1_suptitle, fontsize=16, fontweight="bold")
    fig1.savefig("taxa_summary.png", dpi=300,  transparent=False)
    # ax2 = axes1[1]
    # sns.barplot(x="sample_group", y="phage_coverage_percent", hue="phage_spike", data=df, ax=ax2, palette=YESTERDAY_MEDIUM)
    # ax2.set_title(f"\nT4 Coverage")
    # ax2.set_xlabel("\n")
    # ax2.set_ylabel("Coverage (%)")
    # legend = ax2.get_legend()
    # legend.set_title(None) 
    # sns.despine()
    # ax2.grid(False)


# def get_top_taxa() -> List[pandas.DataFrame]:

#     """
#     Helper function to obtain highest ranked taxa for each sample
#     """