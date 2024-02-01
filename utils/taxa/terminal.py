import typer
import pandas
import warnings
import numpy as np
import seaborn as sns
from pathlib import Path
from dataclasses import dataclass
from matplotlib import pyplot as plt
from typing import List, Optional, Tuple, Dict
from ..utils import read_qc_table, YESTERDAY_MEDIUM

from matplotlib.colors import LinearSegmentedColormap, to_rgb

pandas.set_option('display.max_columns', None)
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
def plot_heatmap_top(
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
def plot_heatmap_target(
    taxa: Path = typer.Option(
        ..., help="Taxa overview summary table"
    ),
    targets: str = typer.Option(
         "me_control_1.txt,me_control_2.txt", help="Target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) "
    ),
    offtargets: str = typer.Option(
         "me_control_misc.txt", help="Off target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) - checked for presence in samples"
    ),
    target_labels: Optional[str] = typer.Option(
        None, help="Off target files with list of taxonomic identifiers to extract (comma-delimited string if multiple files or taxids per line) - checked for presence in samples"
    ),
):
    """
    Plot target and offtarget taxa heatmap for a set of sample taxonomic profiles returned from Cerebro API
    """

    # Set grouping variable
    group = "cerebro_id"

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

    # Check and get target/offtarget list file paths
    target_paths, offtarget_paths = get_file_paths(targets=targets, offtargets=offtargets)
    
    # Read the files into dictionaries
    target_taxids = get_taxa_to_include(files=target_paths)
    offtarget_taxids = get_taxa_to_include(files=offtarget_paths)

    # Create the data for targets
    target_classifications = extract_classifications(
        df=df, 
        group=group, 
        label_taxids=target_taxids, 
    )

    # Create the data for offtargets
    offtarget_classifications = extract_classifications(
        df=df, 
        group=group, 
        label_taxids=offtarget_taxids, 
    )

    # Shape the target data into the final heatmap format
    create_heatmap(library_results=target_classifications, qualitative=False, label_map=label_map)
    create_heatmap(library_results=target_classifications, qualitative=True, label_map=label_map)

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
    assembly_contigs: Optional[int] = None
    assembly_contigs_bases: Optional[int] = None

@dataclass
class LibraryResult:
    group: str
    sample_id: str
    sample_tag: str
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

def extract_classifications(
    df: pandas.DataFrame, 
    group: str,
    label_taxids: Dict[str, List[int]], 
) -> List[LibraryResult]:
    
    classifications = []
    for grp, library in df.groupby(group):
        for label, taxids in label_taxids.items():

            # Is this target in the library profile?
            target_profile = library.loc[library['taxid'].isin(taxids), :]
            

            # Libraries only ever have a unqiue sample identifier
            # but for now we want to group by model and extract it
            # may change on testing
            sample_id = library.sample_id.unique()[0]
            sample_tag = library.sample_tag.unique()[0]

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
                quantitative=quantitative,
                qualitative=qualititative
            ))

    return classifications

def get_file_paths(targets: str, offtargets: str) -> Tuple[List[Path], List[Path]]:

    target_paths = []
    for target_file in targets.strip().split(","):
        if Path(target_file).exists():
            target_paths.append(Path(target_file))
        else:
            print(f"Target file {target_file} could not be found")

    offtarget_paths = []
    for offtarget_file in offtargets.strip().split(","):
        if Path(offtarget_file).exists():
            offtarget_paths.append(Path(offtarget_file))
        else:
            print(f"Offtarget file {offtarget_file} could not be found")

    return target_paths, offtarget_paths


def get_taxa_to_include(files: List[Path]) -> Dict[str, List[int]]:

    taxids = dict()
    for file in files:
        with file.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue

                content = line.strip().split("\t")

                for c in content:
                    if not c:
                        continue
                
                tids = [int(tid.strip()) for tid in content[0].strip().split(",")]
                
                label = content[1].strip()

                taxids[label] = tids

    return taxids
                

def create_heatmap(
    library_results: List[LibraryResult], 
    qualitative: bool = False, 
    label_map: Optional[Dict[str, str]] = None
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
        plot_qualitative_heatmap(data=data)
    else:
        plot_quantitative_heatmap(data=data)

def create_qualitative_dataframe(
    data: Dict[str, List[QualitativeResult]],
    field: str = "present"
):
    column_names = ["Library"]
    rows = []
    for sample_name, quantitative_results in data.items():
        row = [sample_name]

        for result in quantitative_results:
            # Column headers if not already present
            if result.label not in column_names:
                column_names.append(result.label)
            row.append(int(getattr(result, field)))
        rows.append(row)

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

        for result in quantitative_results:
            # Column headers if not already present
            if result.label not in column_names:
                column_names.append(result.label)
            row.append(getattr(result, field))
        rows.append(row)

    df = pandas.DataFrame(rows, columns=column_names).set_index("Library").sort_index()
    
    return df


def plot_qualitative_heatmap(data: Dict[str, List[QualitativeResult]], outfile: Path = Path("qualitative.png")):


    dfs = [
        create_qualitative_dataframe(data=data, field="present"),
        create_qualitative_dataframe(data=data, field="kmer"),
        create_qualitative_dataframe(data=data, field="alignment"),
        create_qualitative_dataframe(data=data, field="assembly")
    ]

    fig1, axes = plt.subplots(nrows=2, ncols=2, figsize=(36,30))
    
    axes = axes.flat

    colors = [
        ["#9DA7BF", "#97AD3D"], ["#9DA7BF", "#55CFD8"], ["#9DA7BF", "#B46DB3"], ["#9DA7BF", "#F9B40E"]
    ]

    classifications = [
        "Combined", "K-mer", "Alignment", "Assembly"
    ]

    for i, df in enumerate(dfs):
        dataframe = df.transpose()
        print(dataframe)
        x_labels = dataframe.columns.tolist()
        ax = axes[i]

        p = sns.heatmap(dataframe, ax=ax, square=False, cmap=colors[i], linewidths=2, robust=True, xticklabels=1, yticklabels=1, cbar=False)
        
        ax.xaxis.tick_top()
        ax.set_xticklabels(x_labels, rotation=45, ha="center")
        ax.set_xlabel(f"{classifications[i]}", fontsize=24)
        ax.tick_params(axis='x', labelsize="large")  
        ax.tick_params(axis='y', labelsize="large") 

    fig1.suptitle('ZeptoMetrix Target Detection (Qualitative)', fontsize=36)

    fig1.savefig(outfile, dpi=300,  transparent=False)


def plot_quantitative_heatmap(data: Dict[str, List[QuantitativeResult]], outfile: Path = Path("quantitative.png")):

    vmax = 1
    vmax_bp = 1000

    dfs = [
        create_quantitative_dataframe(data=data, field="total_rpm"),
        create_quantitative_dataframe(data=data, field="kmer_rpm"),
        create_quantitative_dataframe(data=data, field="alignment_rpm"),
        create_quantitative_dataframe(data=data, field="assembly_contigs_bases")
    ]

    fig1, axes = plt.subplots(nrows=2, ncols=2, figsize=(36,30))
    
    axes = axes.flat

    colors = [
        ["#f9f9f9", "#58a787"], ["#f9f9f9", "#4c849a"], ["#f9f9f9", "#c1bd38"], ["#f9f9f9", "#5d5686"]
    ]

    classifications = [
        "Combined RPM", "K-mer RPM", "Alignment RPM", "Assembly BP"
    ]

    for i, df in enumerate(dfs):
        dataframe = df.transpose()
        dataframe.replace(0, np.nan, inplace=True)

        print(dataframe)
        x_labels = dataframe.columns.tolist()
        ax = axes[i]

        rgb_colors = [to_rgb(color) for color in colors[i]]

        # Create a custom colormap
        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", rgb_colors)

        p = sns.heatmap(dataframe, ax=ax, square=False, cmap=custom_cmap, linewidths=2, robust=True, xticklabels=1, yticklabels=1, cbar=True, vmax=vmax if i < 3 else vmax_bp)
        
        ax.xaxis.tick_top()
        ax.set_xticklabels(x_labels, rotation=90, ha="center")
        ax.set_xlabel(f"{classifications[i]}", fontsize=24)
        ax.tick_params(axis='x', labelsize="large")  
        ax.tick_params(axis='y', labelsize="large") 

        cbar = ax.collections[0].colorbar

        # Get the current tick labels
        labels = [item.get_text() for item in cbar.ax.get_yticklabels()]

        # Modify the label of the top value
        if labels:
            labels[-1] = f'{vmax}+ rpm' if i < 3 else f'{vmax_bp}+ bp'
            labels[0] = f'0 rpm' if i < 3 else f'0 bp'

            cbar.set_ticks(cbar.get_ticks())
            cbar.set_ticklabels(labels)

    plt.subplots_adjust(hspace=0.4)

    fig1.suptitle('ZeptoMetrix Target Detection (Quantitative)', fontsize=36)

    fig1.savefig(outfile, dpi=300,  transparent=False)

    for i, df in enumerate(dfs):
        df.replace(0, np.nan, inplace=True)
        df["Metric"] = [classifications[i] for _ in df.iterrows()]
    
    pandas.concat(dfs).to_csv("quantitative.tsv", sep="\t")