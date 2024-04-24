import typer
import warnings

import typer
import pandas
import warnings
import numpy as np
import seaborn as sns

from typing import Optional, List
from pathlib import Path
from dataclasses import dataclass
from matplotlib import pyplot as plt

from ..utils import read_qc_table, YESTERDAY_MEDIUM
from ..taxa.terminal import get_file_paths, get_taxa_to_include, extract_classifications, LibraryResult, LibraryDataFrame

pandas.set_option('display.max_columns', None)
warnings.filterwarnings("ignore")

def stdin_callback(value: Optional[Path]) -> Path:
    return value if value else Path('/dev/stdin')

app = typer.Typer(add_completion=False)

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

@app.command()
def plot_experiment(
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
):
    
    """
    Plot the extraction validation experiment
    """

    group = "cerebro_id"
    

    # Read the taxon profiles
    df = pandas.read_csv(taxa, sep=",", header=0)

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
            label_targets=label_targets 
        )

        print(target_classifications)

        plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="primer")
        plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="pellet")
        plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction")
        plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="machine")
        plot_comparison(target_classifications=target_classifications, panel=f"ZM{i+1}", library_variable="extraction_machine_primer")




def plot_comparison(target_classifications: List[LibraryResult], panel: str = "ZM1", library_variable: str | None = None):

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

    print(df)

    fig1, axes = plt.subplots(nrows=4, ncols=2, figsize=(24,48))

    axes = axes.flat
    i = 0
    for metric in ("total_rpm", "kmer_rpm", "alignment_rpm", "assembly_contigs_bases"):
        for (nucleic_acid, nucleic_acid_data) in df.groupby("nucleic_acid"):

            nucleic_acid_data = nucleic_acid_data.replace(0, np.nan)
            nucleic_acid_data[metric] = np.log10(nucleic_acid_data[metric])

            print(nucleic_acid_data)
            
            if library_variable == "extraction_machine_primer":
                hue_order = [
                    "sonication-ez1-NEB",
                    "sonication-tanbead-NEB",
                    "beads-ez1-NEB",
                    "beads-tanbead-NEB",
                    "sonication-ez1-BIO",
                    "sonication-tanbead-BIO",
                    "beads-ez1-BIO",
                    "beads-tanbead-BIO",
                ]
                palette = [
                    "#444E7E",
                    "#8087AA",
                    "#B7ABBC",
                    "#F9ECE8",
                    "#D8511D",
                    "#FD8700",
                    "#FEB424",
                    "#FCC893"
                ]
            else:
                hue_order = None
                palette = YESTERDAY_MEDIUM

            p = sns.barplot(nucleic_acid_data, x="label", y=metric, hue=library_variable, hue_order=hue_order, ax=axes[i], palette=palette)
            sns.stripplot(x="label", y=metric, hue=library_variable, hue_order=hue_order, data=nucleic_acid_data, ax=axes[i], palette=palette, dodge=True if library_variable else False, edgecolor="black", linewidth=2, legend=None)
    
            p.set_title(f"{panel} {nucleic_acid} ({metric})")
            p.set_ylabel(f"Log10 ({metric})")
            p.set_xlabel(None)
            legend = p.get_legend()
            legend.set_title(None) 
            sns.move_legend(p, "upper right")


            i += 1

    fig1.savefig(f"{panel}_{library_variable if library_variable else 'all'}.pdf", dpi=300,  transparent=False)
