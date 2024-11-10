import typer
import warnings
import pandas as pd
import importlib.resources

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM

from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import linregress

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

@app.command()
def plot_detroit_spikes(
    qc_table: Path = typer.Option(
        ..., help="Quality control table including ERCC/EDCC data"
    ),
    min_ercc: int = typer.Option(
        None, help="Minimum constructs to include sample from table"
    ),
    biomass_column: str = typer.Option(
        "total_biomass", help="Biomass column to plot on y-axis"
    ),
    detroit_column: str = typer.Option(
        "detroit_input", help="Biomass of Detroit cells input to sample"
    ),
    biomass_label: str = typer.Option(
        "Host Biomass", help="Biomass super title"
    ),
    pipeline_variant: str = typer.Option(
        None, help="Pipeline variant description"
    ),
    constructs: int = typer.Option(
        92, help="Total distinct constructs in the ERCC/EDCC input"
    ),
    prefix: str = typer.Option(
        "detroit_", help="Plot prefix"
    )
):
    
    """
    Plot the Detroit cell spike-in experiment quality control results
    """

    df = read_qc_table(qc_table, remove_ntc=True, min_ercc_constructs=min_ercc)

    df[biomass_column] = df[biomass_column].astype(float)
    df["biomass_error"] = df[biomass_column] - df[detroit_column]
    df["nucleic_acid"] = ["DNA" if "DNA" in identifier else "RNA" for identifier in df.id]

    fig1, axes1 = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
    fig2, axes2 = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
    fig3, axes3 = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
    
    panel_indices = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for (i, (ercc_input_mass, group_df)) in enumerate(df.groupby("ercc_input_mass")):
        
        ax1 = axes1[panel_indices[i][0], panel_indices[i][1]]
        sns.barplot(x=detroit_column, y="ercc_constructs", hue="nucleic_acid", data=group_df, ax=ax1, palette=YESTERDAY_MEDIUM)
        ax1.set_title(f"\nInput: {ercc_input_mass} pg")
        ax1.set_xlabel("Detroit cell input biomass (pg)")
        ax1.set_ylabel("Constructs detected\n")
        ax1.set_ylim(0, constructs)
        legend = ax1.get_legend()
        legend.set_title(None) 
        sns.despine()
        ax1.grid(False)

        ax2 = axes2[panel_indices[i][0], panel_indices[i][1]]
        sns.barplot(x=detroit_column, y=biomass_column, hue="nucleic_acid", data=group_df, ax=ax2, palette=YESTERDAY_MEDIUM)
        sns.stripplot(x=detroit_column, y=biomass_column, hue="nucleic_acid", data=group_df, ax=ax2, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)
        ax2.set_title(f"\nInput: {ercc_input_mass} pg")
        ax2.set_xlabel("Detroit cell input biomass (pg)")
        ax2.set_ylabel("Estimated biomass (pg)\n")
        legend = ax2.get_legend()
        legend.set_title(None) 
        sns.despine()
        ax2.grid(False)

        ax3 = axes3[panel_indices[i][0], panel_indices[i][1]]
        sns.barplot(x=detroit_column, y="biomass_error", hue="nucleic_acid", data=group_df, ax=ax3, palette=YESTERDAY_MEDIUM)
        sns.stripplot(x=detroit_column, y="biomass_error", hue="nucleic_acid", data=group_df, ax=ax3, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)
        ax3.set_title(f"\nInput: {ercc_input_mass} pg")
        ax3.set_xlabel("Detroit cell input biomass (pg)")
        ax3.set_ylabel("Delta estimated - cell input biomass (pg)\n")
        legend = ax3.get_legend()
        legend.set_title(None) 
        sns.despine()
        ax3.grid(False)

    fig4, axes4 = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

    sns.stripplot(x=detroit_column, y="ercc_percent", hue="ercc_input_mass", data=df, palette=YESTERDAY_MEDIUM, ax=axes4, size=16)
    axes4.set_xlabel("\nDetroit cell input biomass (pg)")
    axes4.set_ylabel("Proportion of total reads (%)\n")
    legend = axes4.get_legend()
    legend.set_title("Input (pg)") 
    sns.despine()
    axes4.grid(False)


    fig1_suptitle = "Total constructs" + f"\n{pipeline_variant}" if pipeline_variant else ""
    fig1.suptitle(fig1_suptitle, fontsize=16, fontweight="bold")
    fig2_suptitle = biomass_label + f"\n{pipeline_variant}" if pipeline_variant else ""
    fig2.suptitle(fig2_suptitle, fontsize=16, fontweight="bold")
    fig3_suptitle = biomass_label+" ERROR" + f"\n{pipeline_variant}" if pipeline_variant else ""
    fig3.suptitle(fig3_suptitle, fontsize=16, fontweight="bold")
    fig4_suptitle = f"ERCC/EDCC Proportions" + f"\n{pipeline_variant}" if pipeline_variant else ""
    fig4.suptitle(fig4_suptitle, fontsize=16, fontweight="bold")


    fig1.savefig(f"{prefix}_constructs_detected.png", dpi=300,  transparent=False)
    fig2.savefig(f"{prefix}_biomass_estimated.png", dpi=300,  transparent=False)
    fig3.savefig(f"{prefix}_biomass_error.png", dpi=300,  transparent=False)
    fig4.savefig(f"{prefix}_constructs_proportions.png", dpi=300,  transparent=False)

    plt.close()


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