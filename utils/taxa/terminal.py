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
    species: Path = typer.Argument(
        ..., help="Pathogen detection table for species"
    ),
    metadata: Path = typer.Option(
        ..., help="Reference metadata table"
    ),
    plot_species: str = typer.Option(
        ..., help="List of species names to plot, comma-separated"
    ),
    experiment: str = typer.Option(
        "pool", help="Experiment column subset"
    ),
    output: str = typer.Option(
        "pools.png", help="Plot output"
    ),
    qubit: bool = typer.Option(
        False, help="Host values are not categorical but continuous measurements from Qubit"
    ),
):
    """
    Comparison of pooling strategies (pathogen detection)
    """

    classifiers = ["kraken", "bracken", "metabuli", "ganon"]

    plot_species = [s.strip() for s in plot_species.split(",")]

    species = pandas.read_csv(species, sep="\t", header=0)
    metadata = pandas.read_csv(metadata, sep="\t", header=0)

    # Remove the sample identifier from the sequencing library
    species["id"] = species["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    species = species.merge(metadata, on="id", how="left")

    for nucleic_acid in ("dna", "rna"):

        fig, axes = plt.subplots(
            nrows=len(classifiers), 
            ncols=len(plot_species), 
            figsize=(6 * len(plot_species), 20)
        )

        species_nucleic_acid = species[species["nucleic_acid"] == nucleic_acid]
        species_nucleic_acid = species_nucleic_acid[species_nucleic_acid["experiment"] == experiment]

        for col_index, species_name in enumerate(plot_species):
            species_data = species_nucleic_acid[species_nucleic_acid["name"] == species_name]

            for i, classifier in enumerate(classifiers):
                ax = axes[i][col_index]

                if qubit:

                    # Convert "host" to a numeric column if it is continuous
                    species_data["host_qubit"] = pandas.to_numeric(species_data["host_qubit"], errors="coerce")
                    
                    print(f"Species: {species_name} Classifier: {classifier}")
                    print(species_data)
                    print(f"=============================================================")

                    # Scatterplot for correlation
                    sns.scatterplot(
                        x="host_qubit", y=f"{classifier}_rpm", hue="label",
                        data=species_data, hue_order=["P1", "P2"],
                        ax=ax, palette="deep", edgecolor="black"
                    )

                    # Get colors from the scatterplot palette
                    palette = dict(zip(["P1", "P2"], sns.color_palette("deep", 2)))

                    # Separate regression lines for P1 and P2
                    for label in ["P1", "P2"]:
                        subset = species_data[species_data["label"] == label]
                        sns.regplot(
                            x="host_qubit", y=f"{classifier}_rpm",
                            data=subset, scatter=False, ax=ax,
                            color=palette[label], line_kws={"linewidth": 2, "alpha": 0.8}
                        )


                    ax.set_xlabel("Library Qubit (ng/ul)")
                else:

                    all_hosts = species_nucleic_acid["host_spike"].unique()

                    sns.barplot(
                        x="host_spike", y=f"{classifier}_rpm", hue="label",
                        data=species_data, hue_order=["P1", "P2"],
                        ax=ax, palette=YESTERDAY_MEDIUM,
                        order=all_hosts  # Add this to ensure all hosts are shown
                    )

                    sns.stripplot(
                        x="host_spike", y=f"{classifier}_rpm", hue="label",
                        data=species_data, hue_order=["P1", "P2"],
                        ax=ax, palette=YESTERDAY_MEDIUM, dodge=True,
                        edgecolor="black", linewidth=2, legend=None,
                        order=all_hosts  # Add this to ensure all hosts are shown
                    )

                    ax.set_xlabel("\n")
                
                ax.set_title(f"\n{species_name} ({classifier.capitalize()})")
                ax.set_ylabel(f"{classifier.capitalize()} RPM\n")
                ax.set_ylim(0)

                legend = ax.get_legend()
                if legend:
                    legend.set_title(None)

        fig.suptitle(f"Targets ({nucleic_acid.upper()} libraries)\n", fontsize=18)
        fig.tight_layout()
        fig.savefig(f"{nucleic_acid}_{output}", dpi=300, transparent=False)


@app.command()
def plot_targets(
    species: Path = typer.Argument(
        ..., help="Pathogen detection table for species"
    ),
    output: str = typer.Option(
        "taxa_detection.png", help="Plot taxon detection output"
    ),
    log_scale: bool = typer.Option(False, help="Log scale for plot"),
    ids: str = typer.Option(None, help="Sample identifier start strings to subset dataset"),
    plot_species: str = typer.Option(
        None, help="List of species names to plot, comma-separated"
    ),
    plot_labels: str = typer.Option(
        None, help="List of species labels to plot, comma-separated"
    ),
    exclude_species: str = typer.Option(
        None, help="List of species to exclude, comma-separated"
    )
):
    """
    Simple taxa detection plot across classifiers
    """

    classifiers = ["kraken", "bracken", "metabuli", "ganon", "vircov"]

    species = pandas.read_csv(species, sep="\t", header=0)

    # Remove the sample identifier from the sequencing library
    species["id"] = species["id"].str.replace(r"(_[^_]*)$", "", regex=True)
    
    if ids:
        ids = tuple([s.strip() for s in ids.split(",")])
        species = species[species["id"].str.startswith(ids)]

    if plot_species:
        plot_species = [s.strip() for s in plot_species.split(",")]
    else:
        plot_species = [
            'Aspergillus niger', 
            'Cryptococcus neoformans', 
            "Haemophilus influenzae", 
            "Mycobacterium tuberculosis", 
            "Streptococcus pneumoniae",
            "Toxoplasma gondii",
            'Simplexvirus humanalpha1', 
            'Orthoflavivirus murrayense',
        ]

    if plot_labels:
        plot_labels = [s.strip() for s in plot_labels.split(",")]
    else:
        plot_labels = [
            'ANIG', 
            'CNEO', 
            "HINF", 
            "MTB", 
            "SPNEUMO",
            "TOXO",
            'HSV-1', 
            'MVEV',
        ]

    if len(plot_labels) != len(plot_species):
        raise ValueError("Label and species designations are not of equal length")

    if exclude_species:
        exclude_species = [s.strip() for s in exclude_species.split(",")]

        new_plot_species = []
        exclude_indices = []
        for i, sp in enumerate(plot_species):
            if sp not in exclude_species:
                new_plot_species.append(sp)
            else:
                exclude_indices.append(i)
        
        plot_species = new_plot_species.copy()
        plot_labels = [l for i, l in enumerate(plot_labels) if i not in exclude_indices]

    print(plot_species, plot_labels)

    fig, axes = plt.subplots(
        nrows=len(classifiers), 
        ncols=2, 
        figsize=(6 * 2, 20)
    )

    print(species)

    for ni, nucleic_acid in enumerate(("DNA", "RNA")):
        
        species_nucleic_acid = species[species["id"].str.contains(nucleic_acid)]
        species_data = species_nucleic_acid[species_nucleic_acid["name"].isin(plot_species)]

        for i, classifier in enumerate(classifiers):
            ax = axes[i][ni]

            if log_scale:
                species_data[f"{classifier}_rpm"] = np.log10(species_data[f"{classifier}_rpm"])

            sns.barplot(
                x="name", y=f"{classifier}_rpm", hue=None,
                data=species_data, hue_order=None,
                ax=ax, palette=YESTERDAY_MEDIUM,
                order=plot_species  # Add this to ensure all species are shown
            )

            sns.stripplot(
                x="name", y=f"{classifier}_rpm", hue=None,
                data=species_data, hue_order=None,
                ax=ax, palette=YESTERDAY_MEDIUM, dodge=True,
                edgecolor="black", linewidth=2, legend=None,
                order=plot_species  # Add this to ensure all species are shown
            )

            ax.set_title(f"\n{classifier.capitalize()}")
            ax.set_xlabel("\n")
            ax.set_ylabel(f"{classifier.capitalize()} RPM\n")
            ax.set_ylim(0)

            ax.set_xticklabels(plot_labels, rotation=45, ha="right")

            legend = ax.get_legend()
            if legend:
                legend.set_title(None)

    fig.suptitle(f"Target species (LOD)\n", fontsize=18)
    fig.tight_layout()
    fig.savefig(output, dpi=300, transparent=False)


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


import re
from pathlib import Path
from typing import Optional, List

import typer
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

NAMES_OF_INTEREST = [
    "Truepera radiovictrix",
    "Imtechella halotolerans",
    "Allobacillus halotolerans",
    "Betacoronavirus muris",
    "Orthopoxvirus vaccinia",
    "Saccharomyces cerevisiae",
    "Pneumocystis jirovecii",
]

BACTERIA = {
    "Truepera radiovictrix",
    "Imtechella halotolerans",
    "Allobacillus halotolerans",
}

VIRUSES = {
    "Betacoronavirus muris",
    "Orthopoxvirus vaccinia",
}

EUKARYOTES = {
    "Saccharomyces cerevisiae",
    "Pneumocystis jirovecii",
}

DOMAIN_COLORS = {
    "bacteria": "tab:blue",
    "virus": "tab:orange",
    "eukaryote": "tab:green",
}


def detect_domain(name: str) -> str:
    if name in BACTERIA:
        return "bacteria"
    if name in VIRUSES:
        return "virus"
    if name in EUKARYOTES:
        return "eukaryote"
    return "bacteria"


def parse_substrings(arg: Optional[str]) -> List[str]:
    if not arg:
        return []
    return [s.strip() for s in arg.split(",") if s.strip()]


@app.command()
def plot_grouped(
    csv_path: Path = typer.Argument(..., exists=True, readable=True),
    output: Path = typer.Option("tool_panels.png", "--output", "-o"),
    group1_substrings: Optional[str] = typer.Option(
        None, "--group1", help="Comma-separated list defining group1 membership by id substring match."
    ),
    group2_substrings: Optional[str] = typer.Option(
        None, "--group2", help="Comma-separated list defining group2 membership by id substring match."
    ),
    group1_name: Optional[str] = typer.Option(
        "group1", "--name1", help="Group1 name"
    ),
    group2_name: Optional[str] = typer.Option(
        "group2", "--name2", help="Group2 name"
    ),
):  
    LEGEND_NAMES = {
        "group1": group1_name,
        "group2": group2_name,
    }
    
    df = pd.read_csv(csv_path)

    df = df[df["name"].isin(NAMES_OF_INTEREST)].copy()
    if df.empty:
        raise typer.Exit("No rows match the selected names.")

    TOOLS = [
        "kraken_reads",
        "bracken_reads",
        "metabuli_reads",
        "ganon_reads",
    ]

    # Filter available columns to the allowed tool list
    read_cols = [c for c in TOOLS if c in df.columns]

    if not read_cols:
        raise typer.Exit("None of the specified tool read columns were found in the CSV.")

    g1 = parse_substrings(group1_substrings)
    g2 = parse_substrings(group2_substrings)

    def match_all(substrings, value):
        """Return True if every substring appears in the value."""
        return all(sub in value for sub in substrings)

    df["id_str"] = df["id"].astype(str)
    df["group"] = None

    # Assign groups using AND logic
    if g1:
        df.loc[df["id_str"].apply(lambda x: match_all(g1, x)), "group"] = "group1"

    if g2:
        df.loc[df["id_str"].apply(lambda x: match_all(g2, x)), "group"] = "group2"

    # Keep only rows that matched at least one group
    df = df[df["group"].notna()].copy()


    if df.empty:
        raise typer.Exit("No rows matched either group's substring rules.")

    grouped = df.groupby(["name", "group"], as_index=False)[read_cols].sum()

    tools = read_cols
    n_tools = len(tools)

    fig, axes = plt.subplots(
        n_tools,
        1,
        figsize=(12, 3 * n_tools),
        sharex=True,
        constrained_layout=True,
    )

    if n_tools == 1:
        axes = [axes]

    x_names = [n for n in NAMES_OF_INTEREST if n in grouped["name"].unique()]
    x_positions = np.arange(len(x_names))

    groups = ["group1", "group2"]
    present_groups = [g for g in groups if g in grouped["group"].unique()]
    bar_width = 0.8 / len(present_groups)
    hatches = ["", "//"]

    for ax, tool_col in zip(axes, tools):
        pivot = grouped.pivot(index="name", columns="group", values=tool_col)
        pivot = pivot.reindex(x_names)

        for i, grp in enumerate(present_groups):
            offset = (i - (len(present_groups) - 1) / 2) * bar_width
            heights = pivot[grp].fillna(0).values
            colors = [DOMAIN_COLORS[detect_domain(n)] for n in x_names]

            ax.bar(
                x_positions + offset,
                heights,
                width=bar_width,
                label=LEGEND_NAMES[grp],
                color=colors,
                edgecolor="black",
                hatch=hatches[i],
            )

        ax.set_ylabel(tool_col)
        ax.set_xlabel(None)
        ax.set_title(tool_col)
        ax.set_xticks(x_positions)
        ax.set_yscale("log")
        ax.set_xticklabels(x_names, rotation=45, ha="right")
        ax.grid(axis="y", linestyle="--", alpha=0.4)

    axes[0].legend(title="Group")
    axes[-1].set_xlabel("Taxon")

    fig.savefig(output, dpi=300)
    plt.close(fig)

    grouped.to_csv("grouped.csv", index=False)
