import typer
import warnings

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM

from matplotlib import pyplot as plt
import seaborn as sns

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


