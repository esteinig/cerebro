import typer
import pandas
import warnings
import numpy as np
import seaborn as sns
from pathlib import Path
from dataclasses import dataclass
from matplotlib import pyplot as plt
from typing import List, Optional, Union, Dict
from ..utils import read_qc_table, YESTERDAY_MEDIUM

from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

pandas.set_option('display.max_columns', None)
warnings.filterwarnings("ignore")

def stdin_callback(value: Optional[Path]) -> Path:
    return value if value else Path('/dev/stdin')

app = typer.Typer(add_completion=False)


@app.command()
def plot_species_reference(
    species: Path = typer.Argument(
        ..., help="Pathogen detection table (species)"
    ),
    reference: Path = typer.Option(
        ..., help="Reference table of species and taxid"
    ),
):

    """
    Heatmap of reference classifications for a single sample
    """

    species = pandas.read_csv(species, sep="\t", header=0)
    reference = pandas.read_csv(reference, sep="\t", header=0)

    ref_calls = reference.merge(species[["taxid", "kraken_reads", "bracken_reads", "metabuli_reads", "ganon_reads", "kmcp_reads", "sylph_reads"]], on="taxid", how="left")

    # Drop 'taxid' and '_reads' suffix from columns
    ref_calls = ref_calls.drop(columns=["taxid", "reads"])
    ref_calls.columns = [col.replace("_reads", "") for col in ref_calls.columns]

    # Set 'name' as the index (assuming 'name' is the first column)
    ref_calls = ref_calls.set_index("name")

    # Scale values to be percentages of 668 reads
    ref_calls = (ref_calls / 688) * 100

    plt.figure(figsize=(12, 10))

    sns.heatmap(
        ref_calls, 
        cmap="BuGn", 
        vmin=0, 
        vmax=100, 
        annot=True, 
        fmt=".1f",
    )

    # Customize plot labels and display
    plt.xlabel("")
    plt.ylabel("")
    plt.title("CipherDB v0.1.0 - CNS Syndromic (Rank: Species)")

    plt.tight_layout()
    plt.savefig("cns_true.png", dpi=300, transparent=False)


    # Find distinct species in 'species' that are not in 'reference' by 'taxid'
    distinct_species = species[~species["taxid"].isin(reference["taxid"])]

    # Count the number of distinct species
    num_distinct_species = distinct_species["taxid"].nunique()
    print(f"Number of distinct species in 'species' not in 'reference': {num_distinct_species}")

    # Select only the relevant columns for the heatmap and drop 'taxid' and '_reads' suffix from columns
    distinct_species = distinct_species[["name", "kraken_reads", "bracken_reads", "metabuli_reads", "ganon_reads", "kmcp_reads", "sylph_reads"]]
    distinct_species.columns = [col.replace("_reads", "") for col in distinct_species.columns]


    # Set 'name' as the index for the heatmap
    distinct_species = distinct_species.set_index("name")
    
    # Remove species with fewer than 5 reads across all classifiers
    distinct_species = distinct_species[(distinct_species >= 5).any(axis=1)]

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        distinct_species, 
        cmap="OrRd", 
        annot=True, 
        fmt=".0f",
    )

    plt.xlabel("")
    plt.ylabel("")
    plt.title("CipherDB v0.1.0 - CNS Syndromic (Rank: Species)")

    plt.tight_layout()
    plt.savefig("cns_false.png", dpi=300, transparent=False)



@app.command()
def plot_qc_pools(
    qc_reads: Path = typer.Argument(
        ..., help="Quality control data table"
    ),
    metadata: Path = typer.Option(
        ..., help="Reference metadata table"
    ),
    y_label: str = typer.Option(
        "Synthetic constructs detected (n)", help="Y-axis label"
    ),
    column: str = typer.Option(
        "ercc_constructs", help="Y-axis value column"
    ),
    title: str = typer.Option(
        "ERCC Constructs", help="Subplot title"
    ),
    output: Path = typer.Option(
        "ercc_pool.png", help="Output plot"
    ),
    experiment: str = typer.Option(
        "pool", help="Experiment column subset"
    ),
    label_title: str = typer.Option(
        "Pooling strategy", help="Legent title for label column"
    ),
):
    """
    Comparison of quality control data pooling ERCC/EDCC
    """

    
    controls = pandas.read_csv(qc_reads, sep="\t", header=0)
    metadata = pandas.read_csv(metadata, sep="\t", header=0)

    controls["id"] = controls["id"].str.replace(r'(_[^_]*)$', '', regex=True)

    controls_metadata = controls.merge(metadata, on="id", how="left")
    controls_metadata = controls_metadata[controls_metadata["experiment"] == experiment]

    controls_metadata_dna = controls_metadata[controls_metadata["nucleic_acid"] == "dna"]
    controls_metadata_rna = controls_metadata[controls_metadata["nucleic_acid"] == "rna"]

    fig1, axes = plt.subplots(nrows=1, ncols=2, figsize=(24,12))

    ax1 = axes[0]
    ax2 = axes[1]

    hue_order_dna = sorted(controls_metadata_dna["label"].unique())
    hue_order_rna = sorted(controls_metadata_rna["label"].unique())

    sns.barplot(
        x="host", y=column, hue="label",
        data=controls_metadata_dna, hue_order=hue_order_dna, 
        ax=ax1, palette=YESTERDAY_MEDIUM
    )

    sns.stripplot(
        x="host", y=column, 
        hue="label", data=controls_metadata_dna, 
        hue_order=hue_order_dna, ax=ax1, 
        palette=YESTERDAY_MEDIUM, dodge=True, 
        edgecolor="black", linewidth=2, 
        legend=None
    )

    sns.barplot(
        x="host", y=column, hue="label",
        data=controls_metadata_rna, hue_order=hue_order_rna, 
        ax=ax2, palette=YESTERDAY_MEDIUM
    )

    sns.stripplot(
        x="host", y=column, 
        hue="label", data=controls_metadata_rna, 
        hue_order=hue_order_rna, ax=ax2, 
        palette=YESTERDAY_MEDIUM, dodge=True, 
        edgecolor="black", linewidth=2, 
        legend=None
    )


    ax1.set_title(f"\n{title} (DNA)")
    ax1.set_xlabel("\n")
    ax1.set_ylabel(f"{y_label}\n")
    ax1.set_ylim(0)


    legend = ax1.get_legend()
    if legend:
        legend.set_title(label_title)
        sns.move_legend(ax1, "upper right")

    ax2.set_title(f"\n{title} (RNA)")
    ax2.set_xlabel("\n")
    ax2.set_ylabel(f"{y_label}\n")
    ax2.set_ylim(0)


    legend = ax2.get_legend()
    if legend:
        legend.set_title(label_title)
        sns.move_legend(ax2, "upper right")


    fig1.savefig(output, dpi=300, transparent=False)


@app.command()
def plot_twist_comparison(
    viruses: Path = typer.Argument(
        ..., help="Pathogen detection table for Zeptometrix Viruses (species)"
    ),
    metadata: Path = typer.Option(
        ..., help="Reference metadata table"
    ),
):

    """
    Comparison of virus detection TWIST vs CNS
    """

    names = {
        "Simplexvirus humanalpha1": "HSV-1",
        "Simplexvirus humanalpha2": "HSV-2",
        "Roseolovirus humanbeta6a": "HHV-6a",
        "Roseolovirus humanbeta6b": "HHV-6b",
        "Cytomegalovirus humanbeta5": "CMV",
        "Parechovirus ahumpari": "PEV",
        "Enterovirus betacoxsackie": "EV",
        "Varicellovirus humanalpha3": "VZV"
    }

    viruses = pandas.read_csv(viruses, sep="\t", header=0)
    metadata = pandas.read_csv(metadata, sep="\t", header=0)

    viruses_metadata = viruses.merge(metadata, on="id", how="left")

    viruses_metadata["name"] = viruses_metadata["name"].replace(names)

    print(viruses_metadata["name"])
    
    for panel, data in viruses_metadata.groupby("panel"):

        fig1, axes = plt.subplots(nrows=5, ncols=2, figsize=(20,40))

        viruses_dna = data[data["nucleic_acid"] == "dna"]
        viruses_rna = data[data["nucleic_acid"] == "rna"]

        for i, classifier in enumerate(["kraken", "bracken", "metabuli", "ganon", "kmcp"]):

            ax1 = axes[i][0]

            viruses_dna[f"{classifier}_rpm"] = np.log10(viruses_dna[f"{classifier}_rpm"])
            viruses_rna[f"{classifier}_rpm"] = np.log10(viruses_rna[f"{classifier}_rpm"])

            sns.barplot(
                x="name", y=f"{classifier}_rpm", hue="protocol",
                data=viruses_dna, hue_order=["cns", "twist"], 
                ax=ax1, palette=YESTERDAY_MEDIUM
            )

            sns.stripplot(
                x="name", y=f"{classifier}_rpm", 
                hue="protocol", data=viruses_dna, 
                hue_order=["cns", "twist"], ax=ax1, 
                palette=YESTERDAY_MEDIUM, dodge=True, 
                edgecolor="black", linewidth=2, 
                legend=None
            )

            ax2 = axes[i][1]

            sns.barplot(
                x="name", y=f"{classifier}_rpm", hue="protocol",
                data=viruses_rna, hue_order=["cns", "twist"], 
                ax=ax2, palette=YESTERDAY_MEDIUM
            )

            sns.stripplot(
                x="name", y=f"{classifier}_rpm", 
                hue="protocol", data=viruses_rna, 
                hue_order=["cns", "twist"], ax=ax2, 
                palette=YESTERDAY_MEDIUM, dodge=True, 
                edgecolor="black", linewidth=2, 
                legend=None
            )


            ax1.set_title(f"\nDNA {panel.capitalize()} ({classifier.capitalize()})")
            ax1.set_xlabel("\n")
            ax1.set_ylabel(f"{classifier.capitalize()} RPM (log10)\n")
            ax1.set_ylim(0)

            # ax1.tick_params(axis='x', rotation=45)

            ax2.set_title(f"\nRNA {panel.capitalize()} ({classifier.capitalize()})")
            ax2.set_xlabel("\n")
            ax2.set_ylabel(f"{classifier.capitalize()} RPM(log10)\n")
            ax2.set_ylim(0)
            

            # ax2.tick_params(axis='x', rotation=45)
            
            legend = ax1.get_legend()
            if legend:
                legend.set_title(None)  


            legend = ax2.get_legend()
            if legend:
                legend.set_title(None)  

        fig1.savefig(f"twist_cns_zepto_{panel}.png", dpi=300, transparent=False)
        

@app.command()
def plot_pools(
    viruses: Path = typer.Argument(
        ..., help="Pathogen detection table for Viruses (species)"
    ),
    metadata: Path = typer.Option(
        ..., help="Reference metadata table"
    ),
    experiment: str = typer.Option(
        "pool", help="Experiment column subset"
    ),
    output: str = typer.Option(
        "pools.png", help="Plot output"
    ),
    hsv_species: str = typer.Option(
        "Simplexvirus humanalpha1", help="Plot output"
    ),
):

    """
    Comparison of pooling strategies (pathogen detection)
    """

    viruses = pandas.read_csv(viruses, sep="\t", header=0)
    metadata = pandas.read_csv(metadata, sep="\t", header=0)

    # Remove the sample identifier from the sequencing library
    viruses["id"] = viruses["id"].str.replace(r'(_[^_]*)$', '', regex=True)

    viruses_dna = viruses[viruses["id"].str.contains("__DNA__")]
    viruses_rna = viruses[viruses["id"].str.contains("__RNA__")]

    mve = viruses_rna[viruses_rna["name"] == "Orthoflavivirus murrayense"]
    hsv1 = viruses_dna[viruses_dna["name"] == hsv_species]
    
    mve_metadata = mve.merge(metadata, on="id", how="left")
    hsv1_metadata = hsv1.merge(metadata, on="id", how="left")

    mve_metadata = mve_metadata[mve_metadata["experiment"] == experiment]
    hsv1_metadata = hsv1_metadata[hsv1_metadata["experiment"] == experiment]

    fig1, axes = plt.subplots(nrows=6, ncols=2, figsize=(12,24))
    
    for i, classifier in enumerate(["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]):

        ax1 = axes[i][0]

        sns.barplot(
            x="host", y=f"{classifier}_rpm", hue="label",
            data=mve_metadata, hue_order=["P1", "P2"], 
            ax=ax1, palette=YESTERDAY_MEDIUM
        )

        sns.stripplot(
            x="host", y=f"{classifier}_rpm", 
            hue="label", data=mve_metadata, 
            hue_order=["P1", "P2"], ax=ax1, 
            palette=YESTERDAY_MEDIUM, dodge=True, 
            edgecolor="black", linewidth=2, 
            legend=None
        )

        ax2 = axes[i][1]

        sns.barplot(
            x="host", y=f"{classifier}_rpm", hue="label",
            data=hsv1_metadata, hue_order=["P1", "P2"], 
            ax=ax2, palette=YESTERDAY_MEDIUM
        )

        sns.stripplot(
            x="host", y=f"{classifier}_rpm", 
            hue="label", data=hsv1_metadata, 
            hue_order=["P1", "P2"], ax=ax2, 
            palette=YESTERDAY_MEDIUM, dodge=True, 
            edgecolor="black", linewidth=2, 
            legend=None
        )
        
        ax1.set_title(f"\nOrthoflavivirus murrayense ({classifier.capitalize()})")
        ax1.set_xlabel("\n")
        ax1.set_ylabel(f"{classifier.capitalize()} RPM\n")
        ax1.set_ylim(0)


        ax2.set_title(f"\n{hsv_species} ({classifier.capitalize()})")
        ax2.set_xlabel("\n")
        ax2.set_ylabel(f"{classifier.capitalize()} RPM\n")
        ax2.set_ylim(0)
        legend = ax1.get_legend()
        if legend:
            legend.set_title(None)  


        legend = ax2.get_legend()
        if legend:
            legend.set_title(None)  

    fig1.savefig(f"{output}", dpi=300, transparent=False)


@app.command()
def plot_simulation_evaluation(
    summaries: List[Path] = typer.Argument(
        ..., help="Reference metadata table"
    ),
    pathogens: Path = typer.Option(
        ..., help="Pathogen detection table"
    ),
    match_on: str = typer.Option(
        "name", help="Match expected taxa ('taxid' from Cipher simulation) on this column in the pathogen table"
    ),
    output: Path = typer.Option(
        "simulation_evaluation.png", help="Plot output"
    ),
):

    """
    Evaluate simulations based on 
    """

    # Prepare data

    pathogens = pandas.read_csv(pathogens, sep="\t", header=0)
    
    summary_data = []
    for summary in summaries:
        sim_id = summary.name.split(".")[0]
        sim_name, replicate, coverage = sim_id.split("_")
        
        df = pandas.read_csv(summary, sep="\t", header=0)

        df["sim_id"] = [sim_id for _ in df.iterrows()]
        df["sim_name"] = [sim_name for _ in df.iterrows()]
        df["replicate"] = [replicate for _ in df.iterrows()]
        df["coverage"] = [coverage for _ in df.iterrows()]
        
        summary_data.append(df)

    summary = pandas.concat(summary_data)
    summary['coverage'] = summary['coverage'].apply(convert_coverage)  

    df = summary.merge(
        pathogens, 
        left_on=["sim_id", "taxid"], 
        right_on=["id", match_on], 
        how="left", 
        suffixes=("_sim", "_run")
    )

    df_pathogen = df[df["role"] != "host"]


    # Create a limit of detection by (x-axis = coverage, y-axis = classifier RPM)

    fig1, axes = plt.subplots(nrows=6, ncols=1, figsize=(16,24))
    
    for i, classifier in enumerate(["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]):

        ax1 = axes[i]

        df_pathogen[f"{classifier}_rpm"] = np.log10(df_pathogen[f"{classifier}_rpm"])

        sns.barplot(
            x="coverage", y=f"{classifier}_rpm", hue="id_sim",
            data=df_pathogen, hue_order=sorted(df_pathogen["id_sim"].unique()), 
            ax=ax1, palette=YESTERDAY_MEDIUM
        )

        sns.stripplot(
            x="coverage", y=f"{classifier}_rpm", hue="id_sim", 
            data=df_pathogen, hue_order=sorted(df_pathogen["id_sim"].unique()),
            ax=ax1, palette=YESTERDAY_MEDIUM, 
            dodge=True, edgecolor="black", linewidth=2, legend=None
        )

        ax1.invert_xaxis()
        
        ax1.set_title(f"\n{classifier.capitalize()}")
        ax1.set_xlabel("\n")
        ax1.set_ylabel(f"{classifier.capitalize()}  RPM (log10)\n")
        ax1.set_ylim(0)


        legend = ax1.get_legend()
        if legend:
            legend.set_title(None)  

    fig1.suptitle("Limit of Detection (simulation with Cipher)", fontsize=24)
    fig1.tight_layout(rect=[0, 0, 1, 0.99])
    fig1.savefig(output, dpi=300, transparent=False)

    # Extract unique replicates and coverage values
    unique_replicates = sorted(df_pathogen["replicate"].unique())
    unique_coverages = sorted(df_pathogen["coverage"].unique(), reverse=True)

    # Preprocess Data for Heatmap Plot
    heatmap_data = []

    for classifier in ["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]:
        classifier_data = []
        for sim_id in df_pathogen["id_sim"].unique():
            sim_data = df_pathogen[df_pathogen["id_sim"] == sim_id]
            sim_data = sim_data.copy()
            
            # Compute the proportion of classifier reads over total reads
            sim_data[f"{classifier}_proportion"] = (
                sim_data[f'{classifier}_reads'] / sim_data["reads"]
            )
            
            # Prepare data for heatmap
            pivot_data = sim_data.pivot_table(
                values=f"{classifier}_proportion",
                index="replicate",
                columns="coverage",
                aggfunc="mean"
            )
            
            # Reindex to ensure all replicates and coverage values are present
            pivot_data = pivot_data.reindex(
                index=unique_replicates, 
                columns=unique_coverages, 
                fill_value=0
            )
            
            # Replace zero values with NaN
            pivot_data = pivot_data.replace(np.nan, 0)
            
            classifier_data.append(pivot_data)
        
        heatmap_data.append(classifier_data)
        
    # Create the Heatmap Panel
    n_classifiers = len(["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"])
    n_simulations = len(df_pathogen["id_sim"].unique())

    fig2, axes = plt.subplots(nrows=n_classifiers, ncols=n_simulations, figsize=(16, 24), sharex=True, sharey=True)

    for i, classifier in enumerate(["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]):
        for j, sim_id in enumerate(df_pathogen["id_sim"].unique()):
            ax = axes[i, j]
            
            # Select heatmap data for this classifier and simulation
            heatmap_df = heatmap_data[i][j]
            
            # Create a custom annotation array
            annot = np.where(heatmap_df.values == 0, None, heatmap_df.values)  # Set None for 0, otherwise use the value
            annot = [[f"{v:.2f}" if v else "" for v in row] for row in annot]
            
            sns.heatmap(
                heatmap_df,
                cmap=sns.color_palette("Greens", as_cmap=True),
                cbar=None, # (j == n_simulations - 1),  # Add colorbar only to the last column
                ax=ax,
                linewidths=2,  
                linecolor="white",  
                annot=annot,
                fmt="",
                vmin=0,  # Set minimum value of color scale
                vmax=1   # Set maximum value of color scale
            )
            
            
            if i == 0:
                ax.set_title(f"{sim_id}\n", fontsize=14, fontweight="bold")  # Title for each simulation
            if j == 0:
                ax.set_ylabel(f"{classifier.capitalize()}\n", fontsize=14, fontweight="bold")  # Label for each classifier
                ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
            else:
                ax.set_ylabel(None)

            if i == 5:
                ax.set_xlabel("\nDepth of coverage (x)", fontsize=12, fontweight="bold")  # Add X-axis label for coverage values
                ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
            else:
                ax.set_xlabel(None)

    # Set the overall figure title
    fig2.suptitle("Expected reads recovered (fraction)", fontsize=18)

    # Adjust layout
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(f"{output.stem}_heatmap{output.suffix}", dpi=300, transparent=False)


    sim_taxid_names = summary["taxid"].unique().tolist()

    df = pathogens.copy()

    df_data = []
    for sim_id, data in df.groupby("id"):
        sim_name, replicate, coverage = sim_id.split("_")
        data["sim_name"] = [sim_name for _ in data.iterrows()]
        data["replicate"] = [replicate for _ in data.iterrows()]
        data["coverage"] = [coverage for _ in data.iterrows()]
        df_data.append(data)
    
    df = pandas.concat(df_data)
    df['coverage'] = df['coverage'].apply(convert_coverage)  
    df_pathogen = df.copy() # df[df["name"] != "Homo sapiens"]

    # Initialize DataFrame to store results
    heatmap_data = []

    # Extract unique replicates and coverage values
    unique_replicates = sorted(df_pathogen["replicate"].unique())
    unique_coverages = sorted(df_pathogen["coverage"].unique(), reverse=True)

    # Loop through each classifier
    incorrect_taxa_data = []
    for classifier in ["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]:
        classifier_results = []
        classifier_incorrect_taxa = []
        
        for replicate in unique_replicates:
            replicate_results = []
            replicate_incorrect_taxa = []
            for coverage in unique_coverages:

                # Filter data for this replicate and coverage
                subset = df_pathogen[
                    (df_pathogen["replicate"] == replicate) & (df_pathogen["coverage"] == coverage) & (df_pathogen[f"{classifier}_reads"] > 0)
                ]
                
                # Get incorrect taxa:
                incorrect_taxa = subset[~subset["name"].isin(sim_taxid_names)]

                # Count the number of incorrect taxa
                incorrect_taxa_count = incorrect_taxa.shape[0]

                # Placeholder if no incorrect taxa called
                if incorrect_taxa.empty:
                    incorrect_taxa = pandas.DataFrame({
                        "coverage": [coverage],
                        "replicate": [replicate],
                        "name": [np.nan],
                        f"{classifier}_reads": [0],
                        f"{classifier}_rpm": [0]
                    })

                replicate_incorrect_taxa.append(incorrect_taxa)
                replicate_results.append(incorrect_taxa_count)
            
            classifier_incorrect_taxa.append(pandas.concat(replicate_incorrect_taxa))
            classifier_results.append(replicate_results)
        
        incorrect_taxa_df = pandas.concat(classifier_incorrect_taxa)
        classifier_df = pandas.DataFrame(classifier_results, columns=unique_coverages, index=unique_replicates)

        incorrect_taxa_data.append(incorrect_taxa_df)
        heatmap_data.append(classifier_df)

    heatmap_agg = pandas.concat(heatmap_data)
    highlight_names = ["Homo sapiens", "Simplexvirus"] 

    fig3, axes = plt.subplots(nrows=n_classifiers, ncols=2, figsize=(24, 24))

    for i, classifier in enumerate(["kraken", "bracken", "metabuli", "ganon", "kmcp", "sylph"]):
        ax = axes[i][0]
            
        # Select heatmap data for this classifier and simulation
        heatmap_df = heatmap_data[i]


        # Create a custom annotation array
        annot = np.where(heatmap_df.values == 0, None, heatmap_df.values)  # Set None for 0, otherwise use the value
        annot = [[f"{v:.0f}" if v else ""  for v in row] for row in annot]
        
        sns.heatmap(
            heatmap_df,
            cmap=sns.color_palette("Reds", as_cmap=True),
            cbar=True,  # Add colorbar only to the last column
            ax=ax,
            linewidths=2,  
            linecolor="white",  
            annot=annot,
            fmt="",
            vmin=0,  # Set minimum value of color scale
            vmax=heatmap_agg.max(None)   # Set maximum value of color scale
        )
        
        
        ax.set_title(f"{classifier.capitalize()}\n", fontsize=14, fontweight="bold")
        ax.set_ylabel(f"Unexpected taxa at rank (n)\n", fontsize=12, fontweight="bold")  # Label for each classifier
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")

        if i == 5:
            ax.set_xlabel("\nDepth of coverage (x)", fontsize=12, fontweight="bold")  # Add X-axis label for coverage values
            ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        else:
            ax.set_xlabel(None)

        # Incorrect taxa plot
        incorrect_taxa_df = incorrect_taxa_data[i]
        ax1 = axes[i][1]
        
        incorrect_taxa_df[f"{classifier}_rpm"] = np.log10(incorrect_taxa_df[f"{classifier}_rpm"])

        sns.barplot(
            x="coverage", y=f"{classifier}_rpm", hue="replicate",
            data=incorrect_taxa_df, hue_order=sorted(incorrect_taxa_df["replicate"].unique(), reverse=True), 
            ax=ax1, palette=YESTERDAY_MEDIUM
        )

        stripplot = sns.stripplot(
            x="coverage", y=f"{classifier}_rpm", hue="replicate", 
            data=incorrect_taxa_df, hue_order=sorted(incorrect_taxa_df["replicate"].unique(), reverse=True),
            ax=ax1, palette=YESTERDAY_MEDIUM, 
            dodge=True, edgecolor="black", linewidth=2, legend=None
        )

        ax1.invert_xaxis()
                
        ax1.set_title(f"{classifier.capitalize()}\n", fontsize=14, fontweight="bold")
        ax1.set_ylabel(f"{classifier.capitalize()}  RPM (log10)\n")
        ax1.set_ylim(0)

        if i == 5:
            ax1.set_xlabel("\nDepth of coverage (x)", fontsize=12, fontweight="bold")  # Add X-axis label for coverage values
            ax1.set_xticklabels(ax.get_xticklabels(), rotation=0)
        else:
            ax1.set_xlabel(None)

        legend = ax1.get_legend()
        if legend:
            legend.set_title(None)  

    fig3.tight_layout()
    fig3.savefig(f"{output.stem}_heatmap_error_{output.suffix}", dpi=300, transparent=False)


# Convert the coverage column to floats
def convert_coverage(value: str):
    leading_zeros = len(value) - len(value.lstrip('0'))
    if leading_zeros < 1:  # Handle normal integers
        return float(value)
    else:  # Handle the '01', '001', etc. format
        converted = float(value) / (10 ** leading_zeros)
        return converted

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
def plot_heatmap_toptax(
    taxa: Path = typer.Option(
        ..., help="Taxa overview summary table"
    ),
    group: str = typer.Option(
        "cerebro_id", help="Grouping variable to plot by"
    ),
    metric: str = typer.Option(
        "rpm", help="Metric to plot as color scale (total_rpm, kmer_rpm, alignment_rpm, contigs)"
    ),
    top: int = typer.Option(
        20, help="Extract the `top` ranked taxa by `top-metric` (total_rpm, kmer_rpm, alignment_rpm, contigs)"
    ),
    top_metric: str = typer.Option(
        "rpm", help="Extract the top taxa by `top-metric` (total_rpm, kmer_rpm, alignment_rpm, contigs)"
    ),
    subset: str = typer.Option(
        None, help="Sample name substring to subset"
    ),
    controls: str = typer.Option(
        "NTC", help="Substring in group column to differentiate negative template controls"
    ),
):
    """
    Plot taxa heatmap for a set of sample taxonomic profiles and the top abundant taxa returned from Cerebro API
    """


    df = pandas.read_csv(taxa, sep=",", header=0)

    # Top taxa across the group variable
    tops = []
    for group, group_df in df.groupby(group):
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


@app.command()
def plot_heatmap_targets(
    taxa: Path = typer.Option(
        ..., help="Taxa overview summary table"
    ),
    targets: str = typer.Option(
         "me_control_1.txt,me_control_2.txt", help="Target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) "
    ),
    target_labels: Optional[str] = typer.Option(
        None, help="Off-target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) - checked for presence in samples"
    ),
    vmax_rpm: Optional[int] = typer.Option(
        None, help="RPM for color scale limit - any value above will be at the highest end of the color scale"
    ),
    vmax_bp: Optional[int] = typer.Option(
        None, help="BP for color scale limit - any value above will be at the highest end of the color scale"
    ),
    prefix: str = typer.Option(
        "", help="Prefix for plots and data table output files"
    ),
    transpose: bool = typer.Option(
        False, help="Transpose the heatmap to sample columns"
    ),
    substrings: bool = typer.Option(
        False, help="Target file specifies substrings in organism name to extract instead of taxonomic identifiers"
    ),
    plot_size: str = typer.Option(
        "36,30", help="Plot size as comma-delimited string e.g. 36,30"
    ),
    title: str = typer.Option(
        "Target Detection", help="Plot title"
    ),
    exclude_tag_substring: str = typer.Option(
        "", help="Exclude these samples where substring is in the (joined) sample tag (comma-delimited list as string)"
    ),
    group_by: str = typer.Option(
        "cerebro_id", help="Group taxonomic identification by values in this column"
    ),

):
    """
    Plot target and offtarget taxa heatmap for a set of sample taxonomic profiles returned from Cerebro API
    """

    # Read target labels
    if target_labels:
        label_map = dict()
        for tl in target_labels.split(","):
            sample_id, label = tl.split("::")
            label_map[sample_id] = label
    else:
        label_map = None

    # Read the taxon profiles
    df = pandas.read_csv(taxa, sep=",", header=0)

    print(df)

    # Filter if requested
    if exclude_tag_substring:
        for substring in exclude_tag_substring.split(','):
            df = df[~df['sample_tag'].str.contains(substring.strip())]

    print(df)

    # Check and get target/offtarget list file paths
    target_paths = get_file_paths(targets=targets)
    
    # Read the files into dictionaries
    label_targets = get_taxa_to_include(files=target_paths, substrings=substrings, df=df)

    # Create the data for targets
    target_classifications = extract_classifications(
        df=df, 
        group=group_by, 
        label_targets=label_targets 
    )

    # Shape the target data into the final 
    # heatmap format and plot (quantitative)
    create_heatmap(
        library_results=target_classifications,
        qualitative=False, 
        label_map=label_map, 
        vmax_rpm=vmax_rpm, 
        vmax_bp=vmax_bp, 
        prefix=prefix,
        transpose=transpose,
        plot_size=plot_size
    )

    # Shape the target data into the final 
    # heatmap format and plot (qualitative)
    create_heatmap(
        library_results=target_classifications,
        qualitative=True, 
        label_map=label_map, 
        vmax_rpm=vmax_rpm, 
        vmax_bp=vmax_bp, 
        prefix=prefix,
        transpose=transpose,
        plot_size=plot_size,
        title=title
    )

@dataclass
class QualitativeResult:
    label: str
    present: bool = False
    kmer: bool = False
    alignment: bool = False
    assembly: bool = False

@dataclass
class QuantitativeResult:
    label: str
    present: bool = False
    total_rpm: Optional[float] = None
    kmer_rpm: Optional[float] = None
    alignment_rpm: Optional[float] = None
    remap_rpm: Optional[float] = None
    assembly_contigs: Optional[int] = None
    assembly_contigs_bases: Optional[int] = None


@dataclass
class LibraryData:
    organism:  Optional[str]
    organism_category:  Optional[str]
    organism_domain: Optional[str]
    copies_ul:  Optional[str]
    extraction: Optional[str]
    machine: Optional[str]
    extraction_machine: Optional[str]
    nucleic_acid: Optional[str]
    panel: Optional[str]

@dataclass
class LibraryResult:
    group: str
    sample_id: str
    sample_tag: str
    library_data: Optional[LibraryData]
    qualitative: QualitativeResult
    quantitative: QuantitativeResult

    def get_nucleic_acid(self) -> Optional[str]:

        if "DNA" in self.sample_tag:
            nucleic_acid = "DNA"
        elif "RNA" in self.sample_tag:
            nucleic_acid = "DNA"
        else:
            nucleic_acid = None

        return nucleic_acid
    
    def get_control_tag(self) -> Optional[str]:

        if "NTC" in self.sample_tag:
            control_tag = "NTC"
        elif "ENV" in self.sample_tag:
            control_tag = "ENV"
        else:
            control_tag = None

        return control_tag
    

    def get_primer_tag(self) -> Optional[str]:

        if "BIO" in self.sample_tag:
            primer_tag = "BIO"
        elif "NEB" in self.sample_tag:
            primer_tag = "NEB"
        else:
            primer_tag = None

        return primer_tag


@dataclass
class LibraryDataFrame:
    results: List[LibraryResult]
    dataframe: Optional[pandas.DataFrame]

    def get_dataframe(self, qualitative: bool = False):

        data = []
        for result in self.results:
            row = {
                'sample_id': result.sample_id,
                'sample_tag': result.sample_tag,
                'group': result.group,
            }
            if qualitative:
                row.update(result.qualitative.__dict__)
            else:
                row.update(result.quantitative.__dict__)

            if result.library_data:
                row.update(result.library_data.__dict__)

            data.append(row)
        
        return pandas.DataFrame(data)
    
def extract_classifications(
    df: pandas.DataFrame, 
    group: str,
    label_targets: Dict[str, List[int]], 
    panel: str = None
) -> List[LibraryResult]:
    
    classifications = []
    for grp, library in df.groupby(group):
        for label, taxa in label_targets.items():
            print(grp, label)

            # Is this target in the library profile?
            
            target_profile = library.loc[library['taxid'].isin(taxa), :]
            

            # Libraries only ever have a unique sample identifier
            # but for now we want to group by model and extract it
            # may change on testing
            sample_id = library.sample_id.unique()[0]
            sample_tag = library.sample_tag.unique()[0]

            if library.empty:
                raise ValueError("Empty library!")
            
            try :
                organism = library["organism"].unique()[0]
            except KeyError:
                organism = None
            try :
                organism_category = library["organism_category"].unique()[0]
            except KeyError:
                organism_category = None
            try :
                organism_domain = library["organism_domain"].unique()[0]
            except KeyError:
                organism_domain = None
            try :
                copies_ul = library["copies_ul"].unique()[0]
            except KeyError:
                copies_ul = None
            try :
                extraction = library["extraction"].unique()[0]
            except KeyError:
                if "ST" in sample_id:
                    extraction = "sonication"
                elif "BL" in sample_id:
                    extraction = "baseline"
                elif "BT" in sample_id:
                    extraction = "beads"
                else:
                    extraction = None

            try :
                machine = library["machine"].unique()[0]
            except KeyError:
                if "ST" in sample_id:
                    machine = "tanbead"
                elif "BL" in sample_id:
                    machine = "tanbead"
                elif "BT" in sample_id:
                    machine = "tanbead"
                else:
                    machine = None

            try :
                extraction_machine = library["extraction_machine"].unique()[0]
            except KeyError:
                if "ST" in sample_id:
                    extraction_machine = "sonication_tanbead"
                elif "BL" in sample_id:
                    extraction_machine = "baseline_tanbead"
                elif "BT" in sample_id:
                    extraction_machine = "beads_tanbead"
                else:
                    extraction_machine = None

            try :
                nucleic_acid = library["nucleic_acid"].unique()[0]
            except KeyError:
                if "DNA" in sample_tag:
                    nucleic_acid = "DNA"
                elif "RNA" in sample_tag:
                    nucleic_acid = "RNA"
                else:
                    nucleic_acid = None

            try:
                panel = library["panel"].unique()[0]
                if isinstance(panel, float) and np.isnan(panel):
                    panel = None 
                    for s in ("-S500-", "-S501-", "-S504-", "-S505-", "-S508-", "-S509-"):
                        if s in sample_id:
                            panel = "ZM1"
                    for s in ("-S502-", "-S503-", "-S506-", "-S507-", "-S510-", "-S511-"):
                        if s in sample_id:
                            panel = "ZM2"
            except KeyError:
                panel = None 
                for s in ("-S500-", "-S501-", "-S504-", "-S505-", "-S508-", "-S509-"):
                    if s in sample_id:
                        panel = "ZM1"
                for s in ("-S502-", "-S503-", "-S506-", "-S507-", "-S510-", "-S511-"):
                    if s in sample_id:
                        panel = "ZM2"

            print(panel, sample_id, extraction, extraction_machine)

            if (
                organism, 
                organism_category, 
                organism_domain, 
                copies_ul, 
                extraction, 
                machine, 
                extraction_machine, 
                nucleic_acid
            ) == (None, None, None, None, None, None, None, None):
                library_data = None
            else:
                library_data = LibraryData(
                    organism=organism,
                    organism_category=organism_category,
                    organism_domain=organism_domain,
                    copies_ul=copies_ul,
                    extraction=extraction,
                    machine=machine,
                    extraction_machine=extraction_machine,
                    nucleic_acid=nucleic_acid,
                    panel=panel
                )

            # Quantitative and qualitiative results
            if not target_profile.empty:
                
                qualititative = QualitativeResult(
                    label=label,
                    present=not target_profile.empty,
                    kmer=any(target_profile.kmer),
                    alignment=any(target_profile.alignment),
                    assembly=any(target_profile.assembly)
                )
                quantitative = QuantitativeResult(
                    label=label,
                    present=not target_profile.empty,
                    total_rpm=sum(target_profile.rpm),
                    kmer_rpm=sum(target_profile.rpm_kmer),
                    alignment_rpm=sum(target_profile.rpm_alignment),
                    remap_rpm=sum(target_profile.rpm_alignment),  # change here plz
                    assembly_contigs=sum(target_profile.contigs),
                    assembly_contigs_bases=sum(target_profile.contigs_bases)
                )
            else:
                quantitative, qualititative = \
                    QuantitativeResult(label=label), \
                          QualitativeResult(label=label)

            classifications.append(LibraryResult(
                group=grp,
                sample_id=sample_id,
                sample_tag=sample_tag,
                library_data=library_data,
                quantitative=quantitative,
                qualitative=qualititative
            ))

    return classifications

def get_file_paths(targets: str) -> List[Path]:

    target_paths = []
    for target_file in targets.strip().split(","):
        if Path(target_file).exists():
            target_paths.append(Path(target_file))
        else:
            print(f"Target file {target_file} could not be found")

    return target_paths


def get_taxa_to_include(files: List[Path], substrings: bool = False, df: Optional[pandas.DataFrame] = None) -> Dict[str, List[int]]:

    taxa = dict()
    for file in files:
        with file.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue

                content = line.strip().split("\t")

                for c in content:
                    if not c:
                        continue
                
                identifiers = content[0].strip()
                if identifiers.replace(" ", "") == "":
                    print("Empty identifier found in file, skipping...")
                    continue

                if substrings:
                    if df is None:
                        raise ValueError("You need to supply the taxonomic profile dataframe for substring extractions")
                    
                    # We get the taxids of all possible substring organism matches in the database
                    # and assign the name as label
                    organism_substrings = [name.strip() for name in identifiers.split(",")] # substring organism names
                    regex_pattern = '|'.join(organism_substrings)
                    matches = df.loc[df['name'].str.contains(regex_pattern, na=False), :]
                    print(f"Identifiers: {identifiers}")
                    for i, row in matches.iterrows():
                        label = row['name']
                        taxid = row['taxid']

                        # Matches may have multiple taxonomic identifiers since organism names are not ensured to be unique
                        if label in taxa.keys():
                            # Only unique taxonomic identifiers since there are almost certainly duplicates across many samples
                            if taxid not in taxa[label]:
                                taxa[label] += [taxid]
                        else:
                            taxa[label] = [taxid]
                else:
                    tax = [int(tid.strip()) for tid in identifiers.split(",")]  # taxids
                    label = content[1].strip()
                    taxa[label] = tax

    return taxa
                

def create_heatmap(
    library_results: List[LibraryResult], 
    qualitative: bool = False, 
    label_map: Optional[Dict[str, str]] = None,
    vmax_rpm: Optional[int] = None,
    vmax_bp: Optional[int] = None,
    prefix: str = "",
    transpose: bool = False,
    title: str = "Target Detection",
    plot_size: str = "36,30"
):

    data: Dict[str, List[QuantitativeResult]] = dict()
    for result in library_results:
        
        clean_tag = result.sample_tag
        for substr in ('-S', '-POS', '-PS'):
            clean_tag = clean_tag.replace(substr, "")

        if label_map:
            sample_id = label_map.get(result.sample_id)
            if not sample_id:
                sample_id = result.sample_id
        else:
            sample_id = result.sample_id

        sample_name = f"{sample_id}-{clean_tag}"
        if sample_name in data:
            data[sample_name] += [result.qualitative if qualitative else result.quantitative]
        else:
            data[sample_name] = [result.qualitative if qualitative else result.quantitative]
    
    if qualitative:
        plot_qualitative_heatmap(data=data, prefix=prefix, transpose=transpose, plot_size=plot_size, title=title)
    else:
        plot_quantitative_heatmap(data=data, vmax_rpm=vmax_rpm, vmax_bp=vmax_bp, prefix=prefix, transpose=transpose, plot_size=plot_size, title=title)

def create_qualitative_dataframe(
    data: Dict[str, List[QualitativeResult]],
    field: str = "present"
):
    column_names = ["Library"]
    rows = []
    for sample_name, quantitative_results in data.items():
        row = [sample_name]

        labels_seen = []
        for result in quantitative_results:
            # Column headers if not already present
            if result.label not in column_names:
                column_names.append(result.label)

            if result.label in labels_seen:
                continue
            
            row.append(getattr(result, field))
            labels_seen.append(result.label)

        rows.append(row)

    print(rows)
    print(column_names)
    df = pandas.DataFrame(rows, columns=column_names).set_index("Library").sort_index()
    
    return df


def create_quantitative_dataframe(
    data: Dict[str, List[QuantitativeResult]],
    field: str = "present"
):  
    
    column_names = ["Library"]
    rows = []
    for sample_name, quantitative_results in data.items():
        row = [sample_name]
            
        labels_seen = []
        for result in quantitative_results:
            # Column headers if not already present
            if result.label not in column_names:
                column_names.append(result.label)
            
            if result.label in labels_seen:
                continue
            
            row.append(getattr(result, field))
            labels_seen.append(result.label)
        rows.append(row)

    print(rows)
    print(column_names)
    df = pandas.DataFrame(rows, columns=column_names).set_index("Library").sort_index()
    
    return df


def plot_qualitative_heatmap(
    data: Dict[str, List[QualitativeResult]], 
    prefix: str = "", 
    outdir: Path = Path.cwd(), 
    transpose: bool = False,
    title: str = "Target Detection (Qualitative)",
    plot_size: str = "36,30"
):


    dfs = [
        create_qualitative_dataframe(data=data, field="present"),
        create_qualitative_dataframe(data=data, field="kmer"),
        create_qualitative_dataframe(data=data, field="alignment"),
        create_qualitative_dataframe(data=data, field="assembly")
    ]

    fig1, axes = plt.subplots(nrows=2, ncols=2, figsize=[int(s.strip()) for s in plot_size.split(",")])
    
    axes = axes.flat

    colors = [
        ["#f9f9f9", "#58a787"], ["#f9f9f9", "#4c849a"], ["#f9f9f9", "#c1bd38"], ["#f9f9f9", "#5d5686"]
    ]

    classifications = [
        "Combined", "K-mer", "Alignment", "Assembly"
    ]

    for i, df in enumerate(dfs):
        if transpose:
            dataframe = df.transpose()
        else:
            dataframe = df.copy()

        x_labels = dataframe.columns.tolist()
        ax = axes[i]
        
        p = sns.heatmap(dataframe, ax=ax, square=False, cmap=colors[i], linewidths=2, robust=True, xticklabels=1, yticklabels=1, cbar=False)
        
        ax.xaxis.tick_top()
        ax.set_xticklabels(x_labels, rotation=90, ha="center")
        ax.set_xlabel(f"{classifications[i]}", fontsize=24)
        ax.tick_params(axis='x', labelsize="large")  
        ax.tick_params(axis='y', labelsize="large") 

    fig1.suptitle(title, fontsize=36)

    fig1.savefig(outdir / f"{prefix}qualitative.pdf", dpi=300,  transparent=False)

    for i, df in enumerate(dfs):
        df["Metric"] = [classifications[i] for _ in df.iterrows()]
    
    pandas.concat(dfs).to_csv(outdir / f"{prefix}qualitative.tsv", sep="\t")

def plot_quantitative_heatmap(
    data: Dict[str, List[QuantitativeResult]],  
    vmax_rpm: Optional[int] = None,
    vmax_bp: Optional[int] = None,
    outdir: Path = Path.cwd(),
    prefix: str = "",
    transpose: bool = False,
    title: str = "Target Detection (Quantitative)",
    plot_size: str = "36,30"
):

    dfs = [
        create_quantitative_dataframe(data=data, field="total_rpm"),
        create_quantitative_dataframe(data=data, field="kmer_rpm"),
        create_quantitative_dataframe(data=data, field="alignment_rpm"),
        create_quantitative_dataframe(data=data, field="assembly_contigs_bases")
    ]

    fig1, axes = plt.subplots(nrows=2, ncols=2, figsize=[int(s.strip()) for s in plot_size.split(",")])
    
    axes = axes.flat

    colors = [
        ["#f9f9f9", "#58a787"], ["#f9f9f9", "#4c849a"], ["#f9f9f9", "#c1bd38"], ["#f9f9f9", "#5d5686"]
    ]

    classifications = [
        "Combined RPM", "K-mer RPM", "Alignment RPM", "Assembly BP"
    ]

    for i, df in enumerate(dfs):
        
        if transpose:
            dataframe = df.transpose()
        else:
            dataframe = df.copy()

        dataframe.replace(to_replace=0, value=np.nan, inplace=True)

        x_labels = dataframe.columns.tolist()
        ax = axes[i]

        rgb_colors = [to_rgb(color) for color in colors[i]]

        # Create a custom colormap
        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", rgb_colors)

        dataframe = dataframe.sort_index()
        
        dataframe = dataframe[dataframe.columns].astype(float)  # important

        p = sns.heatmap(
            dataframe, 
            ax=ax, 
            square=False, 
            cmap=custom_cmap, 
            linewidths=2, 
            robust=True, 
            xticklabels=1, 
            yticklabels=1, 
            cbar=True, 
            vmax=vmax_rpm if i < 3 else vmax_bp
        )
        
        ax.xaxis.tick_top()
        ax.set_xticklabels(x_labels, rotation=90, ha="center")
        ax.set_xlabel(f"{classifications[i]}", fontsize=24)
        ax.tick_params(axis='x', labelsize="large")  
        ax.tick_params(axis='y', labelsize="large") 

        cbar = ax.collections[0].colorbar

        # Get the current tick labels
        labels = [item.get_text() for item in cbar.ax.get_yticklabels()]

        # Modify the color scale labels
        if labels:
            labels[-1] = f'{vmax_rpm}+ rpm' if i < 3 else f'{vmax_bp}+ bp'
            labels[0] = f'0 rpm' if i < 3 else f'0 bp'

            cbar.set_ticks(cbar.get_ticks())
            cbar.set_ticklabels(labels)

    plt.subplots_adjust(hspace=0.4)

    fig1.suptitle(title, fontsize=36)

    fig1.savefig(outdir / f"{prefix}quantitative.pdf", dpi=300,  transparent=False)

    for i, df in enumerate(dfs):
        df.replace(to_replace=0, value=np.nan, inplace=True)
        df["Metric"] = [classifications[i] for _ in df.iterrows()]
    
    pandas.concat(dfs).to_csv(outdir / f"{prefix}quantitative.tsv", sep="\t")