import typer
import pandas
import warnings

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM


from matplotlib import pyplot as plt

import seaborn as sns
import sys

from typing import List, Optional

warnings.filterwarnings("ignore")

def stdin_callback(value: Optional[Path]) -> Path:
    return value if value else Path('/dev/stdin')

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
    names: List[str] = typer.Option(
        None, help="Taxon names to plot"
    )
):
    
    """
    Plot taxa metrics across a grouping variable
    """

    df = pandas.read_csv(taxa, sep=",", header=0)

    fig1, axes1 = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
            
    ax1 = axes1
    
    df_taxa = df[df["name"].isin(names)] if names else df

    sns.stripplot(x=group, y=metric, hue="name", data=df_taxa, palette=YESTERDAY_MEDIUM, ax=ax1, size=8, alpha=0.7)
    ax1.set_xlabel(f"\n{group}")
    ax1.set_ylabel(f"{metric}\n")

    ax1.tick_params(axis='x', labelsize=8, rotation=45)
    legend = ax1.get_legend()
    legend.set_title("Taxa") 
    sns.despine()
    ax1.grid(False)

    fig1.savefig("taxa_species_summary.png", dpi=300,  transparent=False)

@app.command()
def plot_taxon(
    taxa: Path = typer.Argument(
        Path("/dev/stdin"), help="Taxa overview summary"
    ),
    tax_name: str = typer.Option(
        None, help="Taxon name to plot"
    ),
    taxid: str = typer.Option(
        None, help="Taxon identifier to plot"
    ),
    group: str = typer.Option(
        "run_date", help="Grouping variable (x-axis)"
    ),
    metric: str = typer.Option(
        "rpm", help="Metric to plot"
    ),
    hue: str =typer.Option(
        "sample_tag", help="Categories to plot"
    ),
    panel: str = typer.Option(
        None, help="Panel variable"
    ),
    highlight: List[str] = typer.Option(
        None, help="Highlight these samples/libraries (sample_id)"
    ),
    highlight_color: str = typer.Option(
        "black", help="Color for highlighted data point edges"
    ),
):
    
    """
    Plot taxa metrics across a grouping variable
    """


    df = pandas.read_csv(taxa, sep=",", header=0)

    panels = len(df[panel].unique()) if panel else 1
    fig1, axes = plt.subplots(nrows=1, ncols=panels, figsize=(6*panels,6))
            
    if taxid:
        df_taxa = df[df["taxid"] == taxid]
    elif tax_name:
        df_taxa = df[df["name"] == tax_name]
    else:
        tops = []
        for g, group_df in df.groupby("cerebro_id"):
            sort_df = group_df.sort_values(metric, ascending=False)
            top_df = sort_df[0:20]
            tops.append(top_df)
            df_taxa = pandas.concat(tops).reset_index(drop=True)
    
    if panel:
        for (i, (p, panel_data)) in enumerate(df_taxa.groupby(panel)):
            ax = axes[i]

            if highlight:
                mask = panel_data['sample_id'].isin(highlight)
                a = sns.stripplot(x=group, y=metric, hue=hue, data=panel_data[~mask], palette=YESTERDAY_MEDIUM, ax=ax, size=8 if taxid or tax_name else 2, alpha=0.7, jitter=True, legend=True if i == panels-1 else None)
                sns.stripplot(x=group, y=metric, hue=hue, data=panel_data[mask], color=highlight_color, edgecolor=highlight_color, linewidth=1, ax=a, size=8 if taxid or tax_name else 2, alpha=0.7, jitter=True, legend=None)
            else:
                sns.stripplot(x=group, y=metric, hue=hue, data=panel_data, palette=YESTERDAY_MEDIUM, ax=ax, size=8, alpha=0.7, jitter=True, legend=True if i == panels-1 else None)

            ax.set_title(f"{p}")
            ax.set_xlabel(f"\n{group}")
            ax.set_ylabel(f"{metric}\n")
            ax.tick_params(axis='x', labelsize=8, rotation=45)
            if i == panels-1:
                legend = ax.get_legend()
                legend.set_title("") 
                sns.move_legend(ax, bbox_to_anchor=(1, 1.02), loc='upper left')
            sns.despine()
            ax.grid(False)
    else:
        ax = axes
        sns.stripplot(x=group, y=metric, hue=hue, data=df_taxa, palette=YESTERDAY_MEDIUM, ax=axes, size=8, alpha=0.7)
        ax.set_title(f"{tax_name if tax_name else taxid} ({metric})")
        ax.set_xlabel(f"\n{group}")
        ax.set_ylabel(f"{metric}\n")
        
        ax.tick_params(axis='x', labelsize=8, rotation=45)
        legend = ax.get_legend()
        legend.set_title("") 
        sns.despine()
        ax.grid(False)

    plt.suptitle(tax_name if tax_name else taxid if taxid else "Top 10")
    fig1.savefig("taxon_summary.png", dpi=300,  transparent=False)    

@app.command()
def plot_heatmap(
    taxa: Path = typer.Option(
        ..., help="Taxa overview summary table"
    ),
    group: str = typer.Option(
        "group", help="Grouping variable to plot by"
    ),
    metric: str = typer.Option(
        "rpm", help="Metric to plot as color scale (total_rpm, kmer-rpm, alignment_rpm, contigs)"
    ),
    top: int = typer.Option(
        20, help="Extract the `top` ranked taxa by `top-metric` (total_rpm, kmer-rpm, alignment_rpm, contigs)"
    ),
    top_metric: str = typer.Option(
        "rpm", help="Extract the top taxa by `top-metric` (total_rpm, kmer-rpm, alignment_rpm, contigs)"
    ),
    subset: str = typer.Option(
        None, help="Sample name substring to subset"
    ),
    controls: str = typer.Option(
        "NTC", help="Substring in group column to differentiate negative template controls"
    ),
):
    """
    Plot taxa heatmap for a set of sample taxonomic profiles returned from Cerebro API
    """


    df = pandas.read_csv(taxa, sep=",", header=0)

    tops = []
    for group, group_df in df.groupby("cerebro_id"):
        sort_df = group_df.sort_values(top_metric, ascending=False)
        top_df = sort_df[0:top]
        tops.append(top_df)
    
    top_df = pandas.concat(tops).reset_index(drop=True)

    
    taxon_presence = {}
    for name_taxid in top_df["name_taxid"].unique():
        for sample_name in top_df["sample_name"].unique():
            search = df[(df["name_taxid"] == name_taxid) & (df["sample_name"] == sample_name)]
            if search.empty:
                metric_value = 0
            else:
                metric_value = search[metric].values[0]

            if sample_name not in taxon_presence.keys():
                taxon_presence[sample_name] = {name_taxid: metric_value}
            else:
                taxon_presence[sample_name][name_taxid] = metric_value

    heatmap_data = []
    for col_name, rows in taxon_presence.items():
        heatmap_data.append(
            pandas.DataFrame([[k, v] for k,v in rows.items()], columns=["taxon", col_name]).set_index("taxon", drop=True)
        )
    
    df = pandas.concat(heatmap_data, axis=1)
    df = df.reindex(sorted(df.columns), axis=1)

    if subset:
        df = df.drop(columns=[c for c in df.columns if subset not in c])
    
    print(df)
    
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(24,18))
    p = sns.heatmap(df, ax=ax1, square=False, cmap="crest", robust=True, xticklabels=1, yticklabels=1)

    fig1.savefig("taxa_heatmap.png", dpi=300,  transparent=False)





