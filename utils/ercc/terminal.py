import typer
import warnings
import pandas as pd
import importlib.resources

from pathlib import Path
from ..utils import YESTERDAY_MEDIUM

import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from typing import List

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

def load_reference(ercc_reference: Path) -> pd.DataFrame:
    if ercc_reference is None:
        with importlib.resources.path('utils.assets', 'ercc.tsv') as ercc_path:
            ercc_reference = ercc_path
            if not ercc_reference.exists():
                raise ValueError(f"ERCC reference path does not exist: {ercc_reference}")
    return pd.read_csv(ercc_reference, sep="\t", header=0)

def clean_id(identifier: str) -> str:
    """Removes trailing string after the last underscore in control id."""
    return identifier.rsplit("_", 1)[0]

def subset_controls_by_metadata(controls: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Subset controls data to only include IDs present in metadata, after cleaning IDs."""
    controls["cleaned_id"] = controls["id"].apply(clean_id)
    return controls[controls["cleaned_id"].isin(metadata["id"])]

def merge_controls_with_metadata(controls: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Merge controls with metadata on cleaned id."""
    controls = controls.drop(columns="id")  # Drop original id to avoid conflicts in merge
    return pd.merge(controls, metadata, left_on="cleaned_id", right_on="id")

def compute_group_r2(ercc_data: pd.DataFrame, group_column: str) -> (dict, list):
    r2_values = {}
    hue_order = sorted(ercc_data[group_column].unique())
    for group in hue_order:
        group_data = ercc_data[ercc_data[group_column] == group]
        if len(group_data) > 1:
            slope, intercept, r_value, p_value, std_err = linregress(
                group_data["concentration"], group_data["alignments"]
            )
            r2_values[group] = r_value ** 2
        else:
            r2_values[group] = None
    hue_order_with_r2 = [
        f"{group} (R² = {r2_values[group]:.2f})" if r2_values[group] is not None else f"{group} (R² = N/A)"
        for group in hue_order
    ]
    ercc_data[f"{group_column}_with_r2"] = ercc_data[group_column].apply(
        lambda g: f"{g} (R² = {r2_values[g]:.2f})" if r2_values[g] is not None else f"{g} (R² = N/A)"
    )
    return r2_values, hue_order_with_r2, ercc_data

def plot_ercc_scatter(ercc_data: pd.DataFrame, hue_order_with_r2: list, constructs_detected: int, total_count: int, group_col: str, ax):
    # Plot scatterplot
    sns.scatterplot(
        data=ercc_data, x="concentration", y="alignments", hue=f"{group_col}_with_r2",
        hue_order=hue_order_with_r2, palette=YESTERDAY_MEDIUM, ax=ax, s=30
    )
    
    # Compute and display total R²
    if len(ercc_data) > 1:
        overall_slope, overall_intercept, overall_r_value, _, _ = linregress(
            ercc_data["concentration"], ercc_data["alignments"]
        )
        total_r2 = overall_r_value ** 2
    else:
        total_r2 = None
    
    summary_text = f"{constructs_detected}/{total_count}, R² = {total_r2:.2f}" if total_r2 is not None else f"{constructs_detected}/{total_count}, R² = N/A"
    ax.text(
        0.95, 0.05, summary_text, transform=ax.transAxes,
        fontsize=12, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray')
    )
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("Alignments (log10)")
    ax.set_xlabel("Concentration (log10)")
    legend = ax.get_legend()
    if legend:
        legend.set_title(None)
        sns.move_legend(ax, "upper left")

def plot_coverage_bar(controls: pd.DataFrame, ax):

    sns.barplot(
        data=controls, x="coverage", y="reference", orient="h",
        palette=YESTERDAY_MEDIUM, ax=ax, legend=False
    )
    ax.set_ylabel(None)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Coverage (%)")
    ax.axvline(80, color="gray", linestyle="--", linewidth=1)  # Vertical line at 80% coverage

@app.command()
def plot_ercc_controls(
    qc_controls: Path = typer.Option(..., help="Quality control table for internal controls"),
    ercc_reference: Path = typer.Option(None, help="ERCC reference data table"),
):
    
    """ERCC quality control plots """

    reference = load_reference(ercc_reference)
    controls = pd.read_csv(qc_controls, sep="\t", header=0)
    
    for library, data in controls.groupby("id"):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

        constructs = data[data["class"] == "ercc"]
        constructs_detected = constructs[constructs["alignments"] > 0]

        print(f"Detected {len(constructs_detected)}/{len(reference)} synthetic constructs (ERCC) for library {library}")

        ercc_data = constructs.merge(reference, on="reference", how="left")
        _, hue_order_with_r2, ercc_data = compute_group_r2(ercc_data, "group")

        plot_ercc_scatter(ercc_data, hue_order_with_r2, len(constructs_detected), len(reference), "group", axes[0])

        controls_data = data[data["class"] == "control"].sort_values(by="reference")
        plot_coverage_bar(controls_data, axes[1])

        fig.tight_layout()
        fig.savefig(f"{library}.png", dpi=300, transparent=False)


@app.command()
def plot_qc_overview(
    qc_reads: Path = typer.Option(..., help="Quality control summary table"),
    metadata: Path  = typer.Option(..., help="Experiment metadata table"),
    column: str = typer.Option("ercc_constructs", help="Data to show for comparison"),
    experiment: str = typer.Option("pool", help="Subset by metadata column for experiment"),
    log_scale: bool = typer.Option(False, help="Log scale for plot"),
    title: str = typer.Option("ERCC Constructs", help="Plot titles"),
    output: Path  = typer.Option("ercc_constructs.png", help="POutput plot file"),
):
    
    """Quality control summary overview from experiment """

    qc = pd.read_csv(qc_reads, sep="\t", header=0)
    metadata = pd.read_csv(metadata, sep="\t", header=0)
    
    # Remove the sample identifier from the sequencing library
    qc["id"] = qc["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    # Merge metadata with QC table
    merged_data = qc.merge(metadata, on="id", how="left")

    # Subset to experiment
    merged_data = merged_data[merged_data["experiment"] == experiment]

    # Remove the repeat identifier from the sequencing library
    merged_data["id"] = merged_data["id"].str.replace(r"__P[0-9]+", "", regex=True)
    merged_data["id"] = merged_data["id"].str.replace(r"__RPT[0-9]+", "", regex=True)

    plot_grouped_qc(merged_data=merged_data, column=column, log_scale=log_scale, title=title, output=output, hue="label")


@app.command()
def plot_pools_qubit_reads(
    qc_reads: Path = typer.Option(..., help="Quality control summary table"),
    metadata: Path  = typer.Option(..., help="Experiment metadata table"),
    experiment: str = typer.Option("pool", help="Subset by metadata column for experiment"),
):
    """Pooling quality control with attention to Qubit values """


    qc = pd.read_csv(qc_reads, sep="\t", header=0)
    metadata = pd.read_csv(metadata, sep="\t", header=0)
    
    # Remove the sample identifier from the sequencing library
    qc["id"] = qc["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    # Merge metadata with QC table
    merged_data = qc.merge(metadata, on="id", how="left")
    
    # Subset to experiment
    merged_data = merged_data[merged_data["experiment"] == experiment]

    # Remove the repeat identifier from the sequencing library
    merged_data["id"] = merged_data["id"].str.replace(r"__P[0-9]+", "", regex=True)
    merged_data["id"] = merged_data["id"].str.replace(r"__RPT[0-9]+", "", regex=True)


    # Convert "host" to a numeric column if it is continuous
    merged_data["host_qubit"] = pd.to_numeric(merged_data["host_qubit"], errors="coerce")
    
    print(merged_data)

    fig, axes = plt.subplots(
        nrows=1, 
        ncols=2, 
        figsize=(20, 12)
    )

    for i, nucleic_acid in enumerate(("dna", "rna")):

        # Ensure a reasonable number of unique markers
        if merged_data["host_spike"].nunique() > 10:
            raise ValueError(f"Too many unique marker values")

        data = merged_data[merged_data["nucleic_acid"] == nucleic_acid]

        ax = axes[i]

        sns.scatterplot(
            x="host_qubit", y=f"input_reads", hue="label", style="host_spike",
            data=data, hue_order=["P1", "P2"],
            ax=ax, palette="deep", edgecolor="black"
        )

        sns.regplot(
            x="host_qubit", y=f"input_reads",
            data=data, scatter=False, ax=ax,
            color="black", line_kws={"linestyle": "dashed"}
        )

        ax.set_xlabel("Library Qubit (ng/ul)")
        ax.set_ylabel(f"Input reads (n)\n")

        ax.set_title(f"{nucleic_acid.upper()}")
        ax.set_ylim(0)

        legend = ax.get_legend()
        if legend:
            legend.set_title(None)

    fig.suptitle(f"Input reads vs. library concentration\n", fontsize=18)
    fig.tight_layout()
    fig.savefig(f"input_reads_qubit_correlation.png", dpi=300, transparent=False)


@app.command()
def plot_internal_controls_overview(
    qc_controls: Path = typer.Option(..., help="Quality control module table controls"),
    metadata: Path  = typer.Option(..., help="Experiment metadata table"),
    column: str = typer.Option("coverage", help="Data to show for comparison"),
    experiment: str = typer.Option("pool", help="Subset by metadata column for experiment"),
    log_scale: bool = typer.Option(False, help="Log scale for plot"),
    title: str = typer.Option("Genome coverage (%)", help="Plot titles"),
    output: Path  = typer.Option("internal_controls.png", help="Output plot file"),
    dna_control: str = typer.Option("T4-DNA", help="Internal DNA control"),
    rna_control: str = typer.Option("MS2-RNA", help="Internal RNA control"),
):
    
    """Controls overview from experiment """

    controls = pd.read_csv(qc_controls, sep="\t", header=0)
    metadata = pd.read_csv(metadata, sep="\t", header=0)
    
    # Remove ERCC from controls table
    controls = controls[~controls["reference"].str.startswith("ERCC")]

    # Remove the sample identifier from the sequencing library
    controls["id"] = controls["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    # Merge metadata with controls table
    merged_data = controls.merge(metadata, on="id", how="left")

    merged_data = merged_data[merged_data["experiment"] == experiment]

    # Remove the repeat identifier from the sequencing library
    merged_data["id"] = merged_data["id"].str.replace(r"__P[0-9]+", "", regex=True)
    merged_data["id"] = merged_data["id"].str.replace(r"__RPT[0-9]+", "", regex=True)
    
    # Output data table for reference checks

    merged_data.to_csv("internal_controls.csv", index=False, header=True)

    plot_grouped_qc(merged_data=merged_data, column=column, log_scale=log_scale, title=title, output=output, hue="label", dna_phage_control=dna_control, rna_phage_control=rna_control)


@app.command()
def plot_positive_controls_overview(
    species: Path = typer.Option(..., help="Species output table"),
    metadata: Path  = typer.Option(..., help="Experiment metadata table"),
    experiment: str = typer.Option("pool", help="Subset by metadata column for experiment"),
    log_scale: bool = typer.Option(False, help="Log scale for plot"),
    title: str = typer.Option("Genome coverage (%)", help="Plot titles"),
    output: Path  = typer.Option("internal_controls.png", help="Output plot file"),
):
    
    """Controls overview from experiment """

    species = pd.read_csv(species, sep="\t", header=0)
    metadata = pd.read_csv(metadata, sep="\t", header=0)

    positive_controls = [
        "Imtechella halotolerans", "Truepera radiovictrix", "Allobacillus halotolerans",
        "Orthopoxvirus vaccinia", "Betacoronavirus muris",
        "Saccharomyces cerevisiae", "Pneumocystis jirovecii"
    ]
    
    # Remove the sample identifier from the sequencing library
    species["id"] = species["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    # Merge metadata with controls table
    merged_data = species.merge(metadata, on="id", how="left")
    merged_data = merged_data[merged_data["experiment"] == experiment]

    # Remove the repeat identifier from the sequencing library
    merged_data["id"] = merged_data["id"].str.replace(r"__P[0-9]+", "", regex=True)
    merged_data["id"] = merged_data["id"].str.replace(r"__RPT[0-9]+", "", regex=True)

    plot_grouped_positive_controls(merged_data=merged_data, log_scale=log_scale, title=title, output=output, positive_controls=positive_controls)

def plot_grouped_positive_controls(merged_data: pd.DataFrame, title: str, output: Path, log_scale: bool = False, positive_controls: List[str] = None):

    custom_colors = [
        '#A0C4FF',  # Soft light blue
        '#6495ED',  # Cornflower blue
        '#4169E1',  # Royal blue (slightly more saturated)
        '#8068dd',  # Muted dark mauve/purple
        '#c1a6eb',  # Soft muted purple
        '#E377C2',  # Soft pinkish-red
        '#FF99CC'   # Soft pastel pink
    ]

    classifiers = ["kraken", "bracken", "metabuli", "ganon"]
    merged_data = merged_data[merged_data["name"].isin(positive_controls)]
    
    # Separate DNA and RNA data
    dna_data = merged_data[merged_data["nucleic_acid"] == "dna"]
    rna_data = merged_data[merged_data["nucleic_acid"] == "rna"]
    
    # Create figure with extra space at bottom for legend
    fig, axes = plt.subplots(nrows=len(classifiers), ncols=2, figsize=(12, 20), gridspec_kw={'bottom': 0.15})
    
    # Store handles and labels for legend
    handles_list = []
    labels_list = []
    
    for i, classifier in enumerate(classifiers):
        if log_scale:
            dna_data[f"{classifier}_rpm"] = np.log10(dna_data[f"{classifier}_rpm"])
            rna_data[f"{classifier}_rpm"] = np.log10(rna_data[f"{classifier}_rpm"])
        
        # DNA
        ax = axes[i][0]
        dna_bp = sns.barplot(
            x="label", y=f"{classifier}_rpm", hue="name",
            data=dna_data, hue_order=positive_controls,
            ax=ax, palette=custom_colors
        )
        sns.stripplot(
            x="label", y=f"{classifier}_rpm", hue="name",
            data=dna_data, hue_order=positive_controls,
            ax=ax, palette=custom_colors, dodge=True,
            edgecolor="black", linewidth=2, legend=False
        )
        ax.set_title(f"\nDNA ({classifier.capitalize()})")
        ax.set_xlabel("\n")
        ax.set_ylabel(f"{classifier.capitalize()} RPM {'(log10)' if log_scale else ''}\n")
        ax.set_ylim(0)
        
        # RNA
        ax = axes[i][1]
        rna_bp = sns.barplot(
            x="label", y=f"{classifier}_rpm", hue="name",
            data=rna_data, hue_order=positive_controls,
            ax=ax, palette=custom_colors
        )
        sns.stripplot(
            x="label", y=f"{classifier}_rpm", hue="name",
            data=rna_data, hue_order=positive_controls,
            ax=ax, palette=custom_colors, dodge=True,
            edgecolor="black", linewidth=2, legend=False
        )
        ax.set_title(f"\nRNA ({classifier.capitalize()})")
        ax.set_xlabel("\n")
        ax.set_ylabel(f"{classifier.capitalize()} RPM {'(log10)' if log_scale else ''}\n")
        ax.set_ylim(0)
        
        # Collect handles and labels from the first subplot (DNA)
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles_list.extend(handles)
            labels_list.extend(labels)
    
    # Remove individual legends
    for ax in axes.flat:
        if ax.get_legend():
            ax.get_legend().remove()
    
    # Create a single legend at the bottom of the figure
    fig.legend(
        handles_list, 
        labels_list, 
        loc='lower center', 
        bbox_to_anchor=(0.5, 0.08),  # Adjust this to position the legend
        ncol=3,  # Arrange legend in one row
        title=None
    )
    
    # Adjust layout and save the figure
    plt.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches='tight')
    fig.suptitle("Positive Controls\n", fontsize=24)
    plt.close(fig)


def plot_grouped_qc(merged_data: pd.DataFrame, title: str, output: Path, column: str = "ercc_constructs", hue: str = "label", log_scale: bool = False, dna_phage_control: str = None, rna_phage_control: str = None, positive_controls: List[str] = None):
    
    if log_scale:
        merged_data[column] = np.log10(merged_data[column])

    # Separate DNA and RNA data
    dna_data = merged_data[merged_data["nucleic_acid"] == "dna"]
    rna_data = merged_data[merged_data["nucleic_acid"] == "rna"]

    # If phage is requested use only the provided phage references

    if dna_phage_control:
        dna_data = dna_data[dna_data["reference"] == dna_phage_control]


    if rna_phage_control:
        rna_data = rna_data[rna_data["reference"] == rna_phage_control]

    hue_order = sorted(merged_data[hue].unique())

    # Create grouped bar plot for DNA
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

    sns.barplot(
        x="host",
        y=column,
        hue=hue,
        hue_order=hue_order,
        data=dna_data,
        ax=axes[0],
        palette="Blues",
        dodge=True,
    )
    sns.stripplot(
        x="host",
        y=column,
        hue=hue,
        hue_order=hue_order,
        data=dna_data,
        ax=axes[0],
        color="black",
        dodge=True,
        edgecolor="gray",
        linewidth=0.5,
        alpha=0.7,
        legend=None,
    )
    
    if column == "ercc_constructs":
        axes[0].axhline(60, color="gray", linestyle="--", linewidth=1)  # Horizontal line
        axes[0].set_ylim(0, 92)

    axes[0].set_title("DNA")
    axes[0].set_ylabel(title)
    axes[0].set_xlabel(None)

    legend = axes[0].get_legend()
    if legend:
        legend.set_title(None)

    # Create grouped bar plot for RNA
    sns.barplot(
        x="host",
        y=column,
        hue=hue,
        hue_order=hue_order,
        data=rna_data,
        ax=axes[1],
        palette="Greens",
        dodge=True,
    )
    sns.stripplot(
        x="host",
        y=column,
        hue=hue,
        hue_order=hue_order,
        data=rna_data,
        ax=axes[1],
        color="black",
        dodge=True,
        edgecolor="gray",
        linewidth=0.5,
        alpha=0.7,
        legend=None,
    )

    if column == "ercc_constructs":
        axes[1].axhline(60, color="gray", linestyle="--", linewidth=1)  # Horizontal line
        axes[1].set_ylim(0, 92)


    axes[1].set_title("RNA")
    axes[1].set_ylabel(title)
    axes[1].set_xlabel(None)

    legend = axes[1].get_legend()
    if legend:
        legend.set_title(None)

    # Adjust layout and save the figure
    fig.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close(fig)


