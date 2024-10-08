import typer
import warnings
import json 

import typer
import pandas
import warnings
import numpy as np
import seaborn as sns

from typing import Optional, List
from pathlib import Path
from dataclasses import dataclass
from matplotlib import pyplot as plt
from matplotlib.axes import Axes


from ..utils import read_qc_table, YESTERDAY_MEDIUM
from ..taxa.terminal import get_file_paths, get_taxa_to_include, extract_classifications, LibraryResult, LibraryDataFrame

# pandas.set_option('display.max_columns', None)
warnings.filterwarnings("ignore")

EXTRACTION_MACHINE_PALETTE = [
    "#D8511D",
    "#FEB424",
    "#FCC893",
    "#444E7E",
    "#8087AA",
    # "#B7ABBC",
    # "#F9ECE8",
    # "#FD8700",
]

HUE_ORDER_EXPERIMENT = [
    "Baseline - EZ1",
    "Sonication - EZ1",
    "Sonication - TanBead",
    "Bead Beating - EZ1",
    "Bead Beating - TanBead"
]


HUE_ORDER_SIMULATION = [
    "SIMULATION",
]

RNA_ORGANISMS = [
    "Murray Valley encephalitis virus",
    "Enterovirus B"
]

def stdin_callback(value: Optional[Path]) -> Path:
    return value if value else Path('/dev/stdin')

app = typer.Typer(add_completion=False)

@app.command()
def plot_experiment(
    taxa_data: Path = typer.Option(
        ..., help="Taxon table filtered"
    ),
    meta_data: Path = typer.Option(
        ..., help="Taxon table filtered"
    ),
    targets: str = typer.Option(
        ..., help="Target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) "
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
    exclude_organism_substring: str = typer.Option(
        "", help="Exclude these samples where substring is in the organism column values (comma-delimited list as string)"
    ),
    exclude_organism_category_substring: str = typer.Option(
        "", help="Exclude these samples where substring is in the organism_category column values (comma-delimited list as string)"
    ),
    exclude_organism_domain_substring: str = typer.Option(
        "", help="Exclude these samples where substring is in the organism_domain column values (comma-delimited list as string)"
    ),
    sep: str = typer.Option(
        "\t", help="QC and metadata table delimiter"
    ),
):
    
    """
    Plot the extraction validation experiment
    """

    group = "sample_id" # target classifications by sample for testing
    

    # Read the taxon profiles
    taxa = pandas.read_csv(taxa_data, sep=sep, header=0)
    meta = pandas.read_csv(meta_data, sep=sep, header=0)

    df = merge_meta_data(meta, other=taxa, other_column="id")

    # Filter if requested
    if exclude_organism_substring:
        for substring in exclude_organism_substring.split(','):
            df = df[~df['organism'].str.contains(substring.strip())]
    
    # Filter if requested
    if exclude_organism_category_substring:
        for substring in exclude_organism_category_substring.split(','):
            df = df[~df['organism_category'].str.contains(substring.strip())]

    # Filter if requested
    if exclude_organism_domain_substring:
        for substring in exclude_organism_domain_substring.split(','):
            df = df[~df['organism_domain'].str.contains(substring.strip())]

    print(df)


    # Check and get target/offtarget list file paths
    target_paths = get_file_paths(targets=targets)
    
    for i, target_path in enumerate(target_paths):
        
        print(target_path)

        # Read the files into dictionaries
        label_targets = get_taxa_to_include(files=[target_path], substrings=substrings, df=df)

        # Extract the data for targets
        target_classifications = extract_classifications(
            df=df, 
            group=group, 
            label_targets=label_targets 
        )

        print(target_classifications)

        plot_comparison(target_classifications=target_classifications, organism="Mycobacterium tuberculosis", dna_only=True, output="test_mycobacterium_tuberculosis.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Enterovirus B", rna_only=True, output="test_enterovirus_b.pdf", remap_rpm=True)        
        plot_comparison(target_classifications=target_classifications, organism="Murray Valley encephalitis virus", rna_only=True, output="test_murray_valley.pdf", remap_rpm=True)
        plot_comparison(target_classifications=target_classifications, organism="Aspergillus niger", dna_only=True, output="test_aspergillus_niger.pdf", zoom_range=["0.1", "0.01"])
        plot_comparison(target_classifications=target_classifications, organism="Candida albicans", dna_only=True, output="test_candida_albicans.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Cryptococcus neoformans", dna_only=True, output="test_cryptococcus_neoformans.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Listeria monocytogenes", dna_only=True, output="test_listeria_monocytogenes.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Haemophilus influenzae", dna_only=True, output="test_haemophilus_influenzae.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Toxoplasma gondii", dna_only=True, output="test_toxoplasma_gondii.pdf")
        plot_comparison(target_classifications=target_classifications, organism="Human mastadenovirus C", dna_only=True, output="test_mastadenovirus_c.pdf", remap_rpm=True)
        plot_comparison(target_classifications=target_classifications, organism="Human alphaherpesvirus 1", dna_only=True, output="test_hsv1.pdf", remap_rpm=True)

        # plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="pellet")
        # plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction")
        # plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="machine")
        # plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction_machine_primer")


@app.command()
def plot_experiment_v1(
    taxa: Path = typer.Option(
        ..., help="Taxon table filtered"
    ),
    targets: str = typer.Option(
         "me_control_1.txt,me_control_2.txt", help="Target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) "
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
    exclude_id_substring: str = typer.Option(
        "", help="Exclude these samples where substring is in the sample identifier (comma-delimited list as string)"
    ),
    sep: str = typer.Option(
        "\t", help="Taxa file delimiter"
    ),
):
    
    """
    Plot the extraction validation experiment
    """

    group = "id"
    
    # Read the taxon profiles
    df = pandas.read_csv(taxa, sep=sep, header=0)

    # Filter if requested
    if exclude_tag_substring:
        for substring in exclude_tag_substring.split(','):
            df = df[~df['sample_tag'].str.contains(substring.strip())]
    
    # Filter if requested
    if exclude_id_substring:
        for substring in exclude_id_substring.split(','):
            df = df[~df['sample_id'].str.contains(substring.strip())]


    # Check and get target/offtarget list file paths
    target_paths = get_file_paths(targets=targets)
    
    for i, target_path in enumerate(target_paths):
        
        print(target_path)

        # Read the files into dictionaries
        label_targets = get_taxa_to_include(files=[target_path], substrings=substrings, df=df)

        # Extract the data for targets
        target_classifications = extract_classifications(
            df=df, 
            group=group, 
            label_targets=label_targets,
        )

        plot_comparison_v1(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction")
        plot_comparison_v1(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="machine")
        plot_comparison_v1(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction_machine")



def plot_comparison_v1(target_classifications: List[LibraryResult], panel: str = "ZM1", library_variable: str | None = None):


    if library_variable is not None:
        results = [
            r for r in target_classifications if r.library_data and r.library_data.panel == panel and r.library_data.__dict__[library_variable]
        ]
    else:
        results = [
            r for r in target_classifications if r.library_data and r.library_data.panel == panel
        ]


    lib = LibraryDataFrame(results=results, dataframe=None)
    df =  lib.get_dataframe(qualitative=False)

    fig1, axes = plt.subplots(nrows=4, ncols=2, figsize=(24,48))

    axes = axes.flat
    i = 0

    print(panel, len(target_classifications), len(results), len(df))

    for metric in ("total_rpm", "kmer_rpm", "alignment_rpm", "assembly_contigs_bases"):
        for (nucleic_acid, nucleic_acid_data) in df.groupby("nucleic_acid"):

            nucleic_acid_data = nucleic_acid_data.replace(0, np.nan)

            if library_variable == "extraction_machine":
                hue_order = [
                    "sonication_tanbead",
                    "beads_tanbead",
                    "baseline_tanbead",
                ]
                44
            else:
                hue_order = None
                palette = YESTERDAY_MEDIUM

            p = sns.barplot(nucleic_acid_data, x="label", y=metric, hue=library_variable, hue_order=hue_order, ax=axes[i], palette=palette)
            p1 = sns.stripplot(x="label", y=metric, hue=library_variable, hue_order=hue_order, data=nucleic_acid_data, ax=axes[i], palette=palette, dodge=True if library_variable else False, edgecolor="black", linewidth=2, legend=None)
    
            p.set_yscale('log')
            p1.set_yscale('log')

            p.set_title(f"{panel} {nucleic_acid} ({metric})")
            p.set_ylabel(f"{metric}")
            p.set_xlabel(None)
            legend = p.get_legend()
            legend.set_title(None) 
            sns.move_legend(p, "upper right")


            i += 1

    df.replace(0, np.nan).to_csv(f"{panel}_{library_variable if library_variable else 'all'}.tsv", sep='\t', index=False, header=True)
    fig1.savefig(f"{panel}_{library_variable if library_variable else 'all'}.pdf", dpi=300,  transparent=False)

@app.command()
def plot_qc(
    qc_data: Path = typer.Option(
        ..., help="Quality control table from workflow output"
    ),
    meta_data: Path = typer.Option(
        ..., help="Experiment metadata with variable and treatment columns"
    ),
    min_length_bp: str = typer.Option(
        "75", help="Threshold used for the correlation plots of minimum read length and adapters trimmed"
    ),
    sep: str = typer.Option(
        "\t", help="QC and metadata table delimiter"
    ),
):
    
    """
    Plot quality control data from the extraction experiment
    """

    qc = pandas.read_csv(qc_data, sep=sep, header=0)
    meta = pandas.read_csv(meta_data, sep=sep,header=0)
    
    df = merge_meta_data(meta, qc)

    print(df)

    print(f"Merged quality control and meta data have {len(df)} matching samples")

    fig1, axes = plt.subplots(nrows=5, ncols=2, figsize=(24,60))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="adapter_trimmed_percent", title="Reads with adapters trimmed (%)", ylabel=None, xlabel="Percent of reads trimmed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="qc_min_length_percent", title=f"Read removed < {min_length_bp}bp (%)", ylabel=None, xlabel="Percent removed (%)", dna_axis=axes[2], rna_axis=axes[3])
    plot_comparison_scatter(df=df, x="adapter_trimmed_percent", y="qc_min_length_percent", hue="machine", style="extraction", title="Reads adapter trimmed vs minimum length removed", xlabel="Reads with adapters trimmed (%)", ylabel=f"Reads removed < {min_length_bp}bp (%)", dna_axis=axes[4], rna_axis=axes[5])
    plot_comparison_barplot(df=df, x="mean_length_r1", title="Mean read length (R1) after QC", ylabel=None, xlabel="Base pairs (n)", dna_axis=axes[6], rna_axis=axes[7])
    plot_comparison_barplot(df=df, x="mean_length_r2", title="Mean read length (R2) after QC", ylabel=None, xlabel="Base pairs (n)", dna_axis=axes[8], rna_axis=axes[9])
    fig1.savefig("qc_trim_lengths.pdf", dpi=300,  transparent=False)

    fig1, axes = plt.subplots(nrows=4, ncols=2, figsize=(24,48))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="deduplicated_percent", title="Deduplication (%)", ylabel=None, xlabel="Deduplicated reads (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="low_complexity_percent", title="Low complexity reads (%)", ylabel=None, xlabel="Low complexity reads (n)", dna_axis=axes[2], rna_axis=axes[3])
    plot_comparison_barplot(df=df, x="ercc_constructs", title="ERCC/EDCC synthetic constructs (n)", ylabel=None, xlabel="Synthetic constructs detected (n)", dna_axis=axes[4], rna_axis=axes[5])
    plot_comparison_barplot(df=df, x="ercc_percent", title="ERCC/EDCC reads (%)", ylabel=None, xlabel="Percent of reads (%)", dna_axis=axes[6], rna_axis=axes[7])
    fig1.savefig("qc_dedup_ercc.pdf", dpi=300,  transparent=False)

    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent", title="Reads removed in host depletion (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3])
    plot_comparison_barplot(df=df, x="output_percent", title="Read remaining for classification (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5])
    fig1.savefig("qc_output.pdf", dpi=300,  transparent=False)


    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent", title="Reads removed in host depletion [outliers removed] (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3], max_x_rna=3.5)
    plot_comparison_barplot(df=df, x="output_percent", title="Read remaining for classification [outliers removed] (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5], max_x_dna=1.0)
    fig1.savefig("qc_output_no_outliers.pdf", dpi=300,  transparent=False)

    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent", y="extraction_machine", hue="organism_domain", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent",  y="extraction_machine", hue="organism_domain", title="Reads removed in host depletion (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3])
    plot_comparison_barplot(df=df, x="output_percent",  y="extraction_machine", hue="organism_domain", title="Read remaining for classification (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5])
    fig1.subplots_adjust(wspace=0.25)
    fig1.savefig("qc_output_organism_domain_hue.pdf", dpi=300,  transparent=False)


    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent", y="extraction_machine", hue="organism_domain", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent",  y="extraction_machine", hue="organism_domain", title="Reads removed in host depletion (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3], max_x_rna=3.5)
    plot_comparison_barplot(df=df, x="output_percent",  y="extraction_machine", hue="organism_domain", title="Read remaining for classification (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5], max_x_dna=1.0)
    fig1.subplots_adjust(wspace=0.25)
    fig1.savefig("qc_output_organism_domain_hue_no_outliers.pdf", dpi=300,  transparent=False)


    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent",  y="organism_domain", hue="extraction_machine", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent",  y="organism_domain", hue="extraction_machine", title="Reads removed in host depletion (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3])
    plot_comparison_barplot(df=df, x="output_percent",  y="organism_domain", hue="extraction_machine",title="Read remaining for classification (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5])
    fig1.subplots_adjust(wspace=0.25)
    fig1.savefig("qc_output_extraction_machine_hue.pdf", dpi=300,  transparent=False)


    fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(24,36))
    axes = axes.flat
    plot_comparison_barplot(df=df, x="qc_percent", y="organism_domain", hue="extraction_machine", title="Reads removed in QC (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[0], rna_axis=axes[1])
    plot_comparison_barplot(df=df, x="host_percent", y="organism_domain", hue="extraction_machine", title="Reads removed in host depletion (%)", ylabel=None, xlabel="Reads removed (%)", dna_axis=axes[2], rna_axis=axes[3], max_x_rna=3.5)
    plot_comparison_barplot(df=df, x="output_percent", y="organism_domain", hue="extraction_machine", title="Read remaining for classification (%)", ylabel=None, xlabel="Reads remaining %)", dna_axis=axes[4], rna_axis=axes[5], max_x_dna=1.0)
    fig1.subplots_adjust(wspace=0.25)
    fig1.savefig("qc_output_extraction_machine_hue_no_outliers.pdf", dpi=300,  transparent=False)

def plot_comparison_barplot(df: pandas.DataFrame, x: str, title: str, dna_axis, rna_axis, y: str = "extraction",  hue: str = "machine",  ylabel: str | None = None, xlabel: str | None = None,max_x_dna: float | None = None, max_x_rna: float | None = None):

    if hue == "extraction_machine":
        palette = [
            "#D8511D",
            "#FEB424",
            "#FCC893",
            "#444E7E",
            "#8087AA",
            # "#B7ABBC",
            # "#F9ECE8",
            # "#FD8700",
        ]
    else:
        palette = YESTERDAY_MEDIUM

    if max_x_dna:
        df_dna = df[(df["nucleic_acid"] == "DNA") & (df[x] < max_x_dna)]
    else:
        df_dna = df[df["nucleic_acid"] == "DNA"]


    if max_x_rna:
        df_rna = df[(df["nucleic_acid"] == "RNA") & (df[x] < max_x_rna)]
    else:
        df_rna = df[df["nucleic_acid"] == "RNA"]

    p = sns.barplot(data=df_dna, x=x, y=y, hue=hue, ax=dna_axis, palette=palette, legend=True)
    p = sns.stripplot(data=df_dna, x=x, y=y, hue=hue, ax=dna_axis, palette=palette, dodge=True, alpha=0.9, edgecolor="black", linewidth=1, legend=None)
    
    p.set_title(f"{title} (DNA)")
    p.set_ylabel(f"{ylabel}" if ylabel else None)
    p.set_xlabel(f"{xlabel}" if xlabel else None)
    legend = p.get_legend()
    legend.set_title(None) 
    sns.move_legend(p, "upper right")

    p = sns.barplot(data=df_rna, x=x, y=y, hue=hue, ax=rna_axis, palette=palette, legend=True)
    p = sns.stripplot(data=df_rna, x=x, y=y, hue=hue, ax=rna_axis, palette=palette, dodge=True, alpha=0.9, edgecolor="black", linewidth=1, legend=None)
    
    p.set_title(f"{title} (RNA)")
    p.set_ylabel(f"{ylabel}" if ylabel else None)
    p.set_xlabel(f"{xlabel}" if xlabel else None)
    legend = p.get_legend()
    legend.set_title(None) 
    sns.move_legend(p, "upper right")

    


def plot_comparison_scatter(df: pandas.DataFrame, x: str, y: str, hue: str | None, style: str | None, title: str, ylabel: str | None, xlabel: str | None, dna_axis, rna_axis):

    p = sns.scatterplot(data=df[df["nucleic_acid"] == "DNA"], x=x, y=y, hue=hue, style=style, ax=dna_axis, palette=YESTERDAY_MEDIUM, s=70, alpha=0.9, edgecolor="black", linewidth=1, legend=True)
    
    p.set_title(f"{title} (DNA)")
    p.set_ylabel(f"{ylabel}" if ylabel else None)
    p.set_xlabel(f"{xlabel}" if xlabel else None)
    
    handles, labels = p.get_legend_handles_labels()
    p.legend(handles=handles[1:3]+handles[4:], labels=labels[1:3]+labels[4:])
    sns.move_legend(p, "upper right")

    p = sns.scatterplot(data=df[df["nucleic_acid"] == "RNA"], x=x, y=y, hue=hue, style=style, ax=rna_axis, palette=YESTERDAY_MEDIUM, s=70, alpha=0.9, edgecolor="black", linewidth=1, legend=True)
    
    p.set_title(f"{title} (RNA)")
    p.set_ylabel(f"{ylabel}" if ylabel else None)
    p.set_xlabel(f"{xlabel}" if xlabel else None)
    
    handles, labels = p.get_legend_handles_labels()
    p.legend(handles=handles[1:3]+handles[4:], labels=labels[1:3]+labels[4:])
    sns.move_legend(p, "upper right")


def plot_comparison(target_classifications: List[LibraryResult], organism: str, dna_only: bool = False, rna_only: bool = False, output: Path = "organims.pdf", remap_rpm: bool = False, zoom_range: List[str] = ["10.0", "1.0", "0.1"], simulation: bool = False):

    if dna_only and rna_only:
        raise ValueError("Cannot specify `dna_only` and `rna_only` at the same time")
    
    # if library_variable is not None:
    #     results = [
    #         r for r in target_classifications if r.library_data and r.library_data.panel == panel and r.library_data.__dict__[library_variable]
    #     ]
    # else:

    results = [
        r for r in target_classifications if r.library_data and r.library_data.organism == organism
    ]

    lib = LibraryDataFrame(results=results, dataframe=None)
    df =  lib.get_dataframe(qualitative=False)

    # Libraries with organism spike-in
    df = df[df["label"] == organism]

    # Replace 0 with NAN to show dropouts more clearly
    df = df.replace(0, np.nan)

    # Make the dilution range a categorical data type
    df["copies_ul"] = df["copies_ul"].astype(str)

    df_dna = df[df["nucleic_acid"] == "DNA"]
    df_rna = df[df["nucleic_acid"] == "RNA"]

    # df_dna = df_dna[df_dna["extraction_machine"] != "Baseline - EZ1"]
    
    if dna_only or rna_only:
        ncols = 3
        nrows = 2
        size = (36,24)
    # else:
    #     ncols = 2
    #     nrows = 3
    #     size = (24,36)

    fig1, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=size)

    for i, metric in enumerate(("kmer_rpm", "remap_rpm" if remap_rpm else "alignment_rpm", "assembly_contigs_bases")):

        if dna_only or rna_only:
            dna_ax = axes[0][i]
            rna_ax = axes[0][i]

            dna_low_ax = axes[1][i]
            rna_low_ax = axes[1][i]
        # else:
        #     dna_ax = axes[i][0]
        #     rna_ax = axes[i][1]

        if metric == "kmer_rpm":
            ylabel = "K-mer RPM\n"
        elif metric == "alignment_rpm":
            ylabel = "Alignment RPM\n"
        else:
            ylabel = "Assembled Bases\n"

        if dna_only:
            p = sns.lineplot(data=df_dna, x="copies_ul", y=metric, hue="extraction_machine", hue_order=HUE_ORDER_SIMULATION if simulation else HUE_ORDER_EXPERIMENT, ax=dna_ax, linewidth=4, palette=EXTRACTION_MACHINE_PALETTE, marker='o', markersize=12)

            p.set_title(f"{organism} (DNA)")
            p.set_xlabel("Genome copies per microliter (copies/ul)")
            p.set_ylabel(ylabel)
            
            legend = p.get_legend()
            if legend:
                legend.set_title(None) 
                sns.move_legend(p, "upper right")

            df_dna_low = df_dna[df_dna["copies_ul"].isin(zoom_range)]

            p = sns.stripplot(data=df_dna_low, x="copies_ul", y=metric, hue="extraction_machine", hue_order=HUE_ORDER_SIMULATION if simulation else HUE_ORDER_EXPERIMENT, ax=dna_low_ax, palette=EXTRACTION_MACHINE_PALETTE, size=12)

            if metric != "assembly_contigs_bases":
                p.axhline(y=10, color='black', linestyle='--')



            p.set_title(f"{organism} (DNA)")
            p.set_xlabel("Genome copies per microliter (copies/ul)")
            p.set_ylabel(ylabel)
            p.set_ylim(0)
            
            legend = p.get_legend()
            if legend:
                legend.set_title(None) 
                sns.move_legend(p, "upper right")

        if rna_only:
            p = sns.lineplot(data=df_rna, x="copies_ul", y=metric, hue="extraction_machine", hue_order=HUE_ORDER_SIMULATION if simulation else HUE_ORDER_EXPERIMENT, ax=rna_ax, linewidth=4, palette=EXTRACTION_MACHINE_PALETTE, marker='o', markersize=12)

            p.set_title(f"{organism} (RNA)")
            p.set_xlabel("Genome copies per microliter (copies/ul)")
            p.set_ylabel(ylabel)
            
            legend = p.get_legend()
            if legend:
                legend.set_title(None) 
                sns.move_legend(p, "upper right")
            
            df_rna_low = df_rna[df_rna["copies_ul"].isin(zoom_range)]

            p = sns.stripplot(data=df_rna_low, x="copies_ul", y=metric, hue="extraction_machine", hue_order=HUE_ORDER_SIMULATION if simulation else HUE_ORDER_EXPERIMENT, ax=rna_low_ax, palette=EXTRACTION_MACHINE_PALETTE, size=12)
            
            if metric != "assembly_contigs_bases":
                p.axhline(y=10, color='black', linestyle='--')

            p.set_title(f"{organism} (RNA)")
            p.set_xlabel("Genome copies per microliter (copies/ul)")
            p.set_ylabel(ylabel)
            p.set_ylim(0)
            
            legend = p.get_legend()
            if legend:
                legend.set_title(None) 
                sns.move_legend(p, "upper right")
    
    fig1.savefig(output, dpi=300,  transparent=False)
    df.to_csv(Path(output).with_suffix(".csv"), index=False, header=True)


    # fig1, axes = plt.subplots(nrows=4, ncols=2, figsize=(24,48))

    # axes = axes.flat
    # i = 0
    # for metric in ("total_rpm", "kmer_rpm", "alignment_rpm", "assembly_contigs_bases"):
    #     for (nucleic_acid, nucleic_acid_data) in df.groupby("nucleic_acid"):

    #         nucleic_acid_data = nucleic_acid_data.replace(0, np.nan)
    #         nucleic_acid_data[metric] = np.log10(nucleic_acid_data[metric])

    #         print(nucleic_acid_data)
            
    #         if library_variable == "extraction_machine_primer":
    #             hue_order = [
    #                 "sonication-ez1-NEB",
    #                 "sonication-tanbead-NEB",
    #                 "beads-ez1-NEB",
    #                 "beads-tanbead-NEB",
    #                 "sonication-ez1-BIO",
    #                 "sonication-tanbead-BIO",
    #                 "beads-ez1-BIO",
    #                 "beads-tanbead-BIO",
    #             ]
    #             palette = [
    #                 "#444E7E",
    #                 "#8087AA",
    #                 "#B7ABBC",
    #                 "#F9ECE8",
    #                 "#D8511D",
    #                 "#FD8700",
    #                 "#FEB424",
    #                 "#FCC893"
    #             ]
    #         else:
    #             hue_order = None
    #             palette = YESTERDAY_MEDIUM

    #         p = sns.barplot(nucleic_acid_data, x="label", y=metric, hue=library_variable, hue_order=hue_order, ax=axes[i], palette=palette)
    #         sns.stripplot(x="label", y=metric, hue=library_variable, hue_order=hue_order, data=nucleic_acid_data, ax=axes[i], palette=palette, dodge=True if library_variable else False, edgecolor="black", linewidth=2, legend=None)
    
    #         p.set_title(f"{panel} {nucleic_acid} ({metric})")
    #         p.set_ylabel(f"Log10 ({metric})")
    #         p.set_xlabel(None)
    #         legend = p.get_legend()
    #         legend.set_title(None) 
    #         sns.move_legend(p, "upper right")


    #         i += 1

    # fig1.savefig(f"{panel}_{library_variable if library_variable else 'all'}.pdf", dpi=300,  transparent=False)


def merge_meta_data(meta: pandas.DataFrame, other: pandas.DataFrame, other_column: str = "id") -> pandas.DataFrame:


    # Initialize an empty list to store filtered DataFrames
    filtered_dfs = []

    # Iterate over each unique sample_id in meta
    for sample_id in meta['sample_id'].unique():
        # Filter qc where 'id' starts with the current 'sample_id'
        mask = other[other_column].str.startswith(sample_id)
        filtered_df = other[mask].copy()
        filtered_df['sample_id'] = sample_id  # Add 'sample_id' for merging
        filtered_dfs.append(filtered_df)

    # Concatenate all filtered DataFrames
    concatenated_df = pandas.concat(filtered_dfs, ignore_index=True)

    # Merge with meta on 'sample_id' to bring in the descriptions
    result = pandas.merge(concatenated_df, meta, on='sample_id')

    return result



@app.command()
def plot_reference_alignments(
    workflows: List[Path] = typer.Argument(
        ..., help="Reference alignment validation workflow output folder"
    ),
    meta_data: Path = typer.Option(
        ..., help="Extraction experiment meta data table (.tsv)"
    ),     
    quality_control: bool = typer.Option(
        False, help="Read quality control conducted before reference mapping"
    ),   
    simulation: bool = typer.Option(
        False, help="Experiment was conducted with Cipher simulated reads"
    ),
    log_scale: bool = typer.Option(
        False, help="Log scale for metric to plot (y-axis)"
    ),
    metric: str = typer.Option(
        "rpm", help="Metric to plot, must be one of: rpm, alignments, coverage"
    ),
    title: str = typer.Option(
        "Single reference alignments", help="Plot title"
    ),
    title_size: int = typer.Option(
        48, help="Plot title fontsize"
    ),
    output: Path = typer.Option(
        "reference_alignments.pdf", help="Output plot file"
    ),

):
    
    """
    Plot the reference alignment validation results
    """

    meta = pandas.read_csv(meta_data, sep="\t", header=0)

    workflow_data = []
    for workflow in workflows:
        total_reads_dataframe = get_total_reads(workflow=workflow, quality_control=quality_control)
        alignment_summary_dataframe = get_alignment_summaries(workflow=workflow)

        data = pandas.merge(total_reads_dataframe, alignment_summary_dataframe, on="sample_id", how="left")

        df = pandas.merge(meta, data, on="sample_id", how="left")

        df["rpm"] = (df["alignments"] / df["total_reads"])*1e06
        df["workflow"] = [workflow.name for _ in df.iterrows()]

        workflow_data.append(df)

    df = pandas.concat(workflow_data)

    dataframes = {}
    for ((organism, nucleic_acid), data) in df.groupby(["organism", "nucleic_acid"]):
        if not simulation:
            if organism in RNA_ORGANISMS and nucleic_acid == "RNA":
                dataframes[organism] = data
            elif organism not in RNA_ORGANISMS and nucleic_acid == "DNA":
                dataframes[organism] = data
        else:
            dataframes[organism] = data
        
    ncols = 3
    nrows = 4
    size = (36,48)
    
    fig1, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=size)
    axes = axes.flat

    for (i, (organism, df)) in enumerate(dataframes.items()):
        if organism == "Mycobacterium tuberculosis" and not simulation:
            df = df[df["extraction_machine"] != "Baseline - EZ1"]
            
        plot_reference_alignment_data(df=df, ax=axes[i], metric=metric, title=organism, simulation=simulation, log_scale=log_scale)

    fig1.suptitle(title, fontweight="bold", fontsize=title_size)
    fig1.savefig(output, dpi=300, transparent=False)


def plot_reference_alignment_data(df: pandas.DataFrame, ax: Axes, metric: str = "rpm", title: str = "Organism", simulation: bool = False, log_scale: bool = False):
    
    # Replace 0 with NAN to show dropouts more clearly
    df = df.replace(0, np.nan)

    # Make the dilution range a categorical data type
    df["copies_ul"] = df["copies_ul"].astype(str)

    if log_scale:
        df[metric] = np.log10(df[metric])

    print(df)

    p = sns.lineplot(data=df, x="copies_ul", y=metric, hue="workflow", ax=ax, linewidth=4, palette=EXTRACTION_MACHINE_PALETTE, marker='o', markersize=12)

    p = sns.stripplot(data=df, x="copies_ul", y=metric, hue="workflow",  ax=ax, palette=EXTRACTION_MACHINE_PALETTE, size=12)

    p.set_title(title)
    p.set_xlabel("Genome copies per microliter (copies/ul)")
    p.set_ylabel(f"{metric.capitalize()}{'(log10)' if log_scale else ''}")
    p.set_ylim(0)
    
    legend = p.get_legend()
    if legend:
        legend.set_title(None) 
        sns.move_legend(p, "upper right")

    return p

def get_total_reads(workflow: Path, quality_control: bool) -> pandas.DataFrame:

    if quality_control:
        fastp_scan_outputs = workflow/"workflow"/"quality_control"/"read_qc"/"scan"
    else:
        fastp_scan_outputs = workflow/"workflow"/"validation"/"reference_alignment"
    
    data = []
    for file in fastp_scan_outputs.glob("*.json"):
        sample_id = file.stem
        total_reads = read_fastp_total_reads(file=file)
        data.append([sample_id, total_reads])

    return pandas.DataFrame(data, columns=["sample_id", "total_reads"])

def get_alignment_summaries(workflow: Path) -> pandas.DataFrame:

    vircov_outputs = workflow/"workflow"/"validation"/"reference_alignment"
    
    data = []
    for file in vircov_outputs.glob("*.tsv"):
        sample_id = file.stem
        genome_size, alignments, coverage = read_vircov_summary(file=file)
        data.append([sample_id, genome_size, alignments, coverage])

    return pandas.DataFrame(data, columns=["sample_id", "genome_size", "alignments", "coverage"])


def read_fastp_total_reads(file: Path):

    with file.open() as fin:
        data = json.load(fin)

    return data["summary"]["before_filtering"]["total_reads"]

def read_vircov_summary(file: Path):

    df = pandas.read_csv(file, sep="\t", header=None, names=[
        "accession", "regions", "reads", "alignments", "bp", "refsize", "coverage", "header", "other"
    ], index_col=None)

    genome_size, alignments, coverage = sum(df["refsize"]), sum(df["alignments"]), sum(df["coverage"])

    return genome_size, alignments, coverage

