import typer
import pandas
import numpy as np
import seaborn as sns

from pathlib import Path
from matplotlib import pyplot as plt

app = typer.Typer(add_completion=False)

@app.command()
def plot_experiment(
    qc: Path = typer.Option(
        ..., help="Quality control summary table (.tsv)"
    ), 
    alignment: Path = typer.Option(
        ..., help="Vircov alignment summary table (.tsv)"
    ),
    dilutions: Path = typer.Option(
        ..., help="Dilution table for samples (.csv)"
    ),
    output: str = typer.Option(
        "lod", help="Output plots basename"
    ),
    aligner: str = typer.Option(
        "bowtie2", help="Aligner used in alignment detection"
    ),
):
    
    """
    LOD experiment plots (September 2024)
    """

    pandas.set_option('display.width', 350)
    pandas.set_option('display.max_columns', 8) 
    pandas.set_option('display.max_rows', None)  

    taxid_reference = [
        'Aspergillus niger', 
        'Cryptococcus neoformans', 
        "Haemophilus influenzae", 
        "Mycobacterium tuberculosis", 
        "Streptococcus pneumoniae",
        "Toxoplasma gondii",
        'HSV-1', 
        'MVEV',
    ]

    qc_df = pandas.read_csv(qc, sep="\t", header=0)
    alignment_df = pandas.read_csv(alignment, sep="\t", header=0)
    dilutions_df = pandas.read_csv(dilutions, sep=",", header=0)


    # Fix dataframe identifiers
    qc_df["id"] = qc_df["id"].str.replace(r'_S\d+$', '', regex=True)
    alignment_df["id"] = alignment_df["id"].str.replace(r'_S\d+$', '', regex=True)
    dilutions_df["id"] = dilutions_df["id"].str.replace(" ", "")
    dilutions_df["dilution"] = dilutions_df["dilution"].str.replace(" ", "")
    

    # Get reference data and combined alignment dataframes
    data_df = pandas.merge(dilutions_df, qc_df, on='id', how='outer')
    alignment_df = pandas.merge(data_df, alignment_df, on='id', how='outer')

    # Fix repeat annotations and identifiers
    alignment_df['repeat'] = alignment_df['id'].str.extract(r'RPT(\d+)$') 
    alignment_df['id'] = alignment_df['id'].str.replace(r'__RPT\d+$', '', regex=True)
    data_df['repeat'] = data_df['id'].str.extract(r'RPT(\d+)$')
    data_df['id'] = data_df['id'].str.replace(r'__RPT\d+$', '', regex=True)

    # Expand with expected taxid from experiment
    reference_df = data_df.loc[data_df.index.repeat(len(taxid_reference))].copy().reset_index(drop=True)
    reference_df['taxid'] = taxid_reference * len(data_df)

    # Get dilution dataframes and concat
    data_plot_dilutions = []
    data_full_dilutions = []
    for dilution in data_df["dilution"].unique():    
        print(f"Getting data for dilution: {dilution}")

        (df_plot, df_full) = get_dilution_data(alignment_df, reference_df, dilution=dilution)
        data_plot_dilutions.append(df_plot)
        data_full_dilutions.append(df_full)
    
    df_plot_dilutions = pandas.concat(data_plot_dilutions, axis=0, ignore_index=True)
    df_full_dilutions = pandas.concat(data_full_dilutions, axis=0, ignore_index=True)

    # Setup plots
    fig, axes = plt.subplots(nrows=len(taxid_reference), ncols=2, figsize=(24,48))
    fig_log, axes_log = plt.subplots(nrows=len(taxid_reference), ncols=2, figsize=(24,48))

    dilution_order = data_df["dilution"].unique().tolist()
    
    row = 0 
    for (taxid, taxid_data) in df_plot_dilutions.groupby('taxid'):
        print(f"Creating plot for taxid: {taxid}")

        dna = taxid_data[taxid_data["id"].str.contains("__DNA__")]
        rna = taxid_data[taxid_data["id"].str.contains("__RNA__")]

        dna_ax = axes[row][0]
        rna_ax = axes[row][1]
        dna_ax_log = axes_log[row][0]
        rna_ax_log = axes_log[row][1]
        
        create_dilution_plot(dna, dna_ax, f"{taxid} (DNA)", aligner, dilution_order, log_scale=False)
        create_dilution_plot(rna, rna_ax, f"{taxid} (RNA)", aligner, dilution_order, log_scale=False)

        create_dilution_plot(dna, dna_ax_log, f"{taxid} (DNA)", aligner, dilution_order, log_scale=True)
        create_dilution_plot(rna, rna_ax_log, f"{taxid} (RNA)", aligner, dilution_order, log_scale=True)

        row += 1
    
    fig.savefig(f"{output}.pdf", dpi=300, transparent=False)
    fig_log.savefig(f"{output}_log.pdf", dpi=300, transparent=False)

    reorder_columns(reference_df, ["id", "taxid", "repeat"]).to_csv(f"{output}_reference.csv", index=False)    

    reorder_columns(df_full_dilutions, ["id", "taxid", "dilution", "repeat", "scan_alignments_total", "remap_alignments_total", "rpm_scan_total", "rpm_remap_total", "input_reads"]).to_csv(f"{output}_data_full.csv", index=False)
    reorder_columns(df_plot_dilutions, ["id", "taxid", "dilution", "repeat", "scan_alignments_total", "remap_alignments_total", "rpm_scan_total", "rpm_remap_total", "input_reads"]).to_csv(f"{output}_data_plot.csv", index=False)

    
def reorder_columns(df, priority_columns):
    
    # Create a list of the remaining columns, excluding the priority columns
    remaining_columns = [col for col in df.columns if col not in priority_columns]

    # Reorder the DataFrame columns
    df = df[priority_columns + remaining_columns]

    return df


def create_dilution_plot(taxid_data, ax, title, aligner, dilution_order, log_scale):

    # Drop data points without detection
    taxid_data = taxid_data.dropna(subset=['rpm_scan_total'])

    print(taxid_data)

    reversed_palette = sns.color_palette("BuGn", 11)[::-1]

    p = sns.barplot(taxid_data, x="dilution", y="rpm_scan_total", ax=ax, palette=reversed_palette, order=dilution_order)
    p1 = sns.stripplot(taxid_data, x="dilution", y="rpm_scan_total", ax=ax, palette=reversed_palette, edgecolor="black", linewidth=2, legend=None, order=dilution_order)

    if log_scale:
        p.set_yscale('log')
        p1.set_yscale('log')

    p.set_title(f"{title}")
    p.set_ylabel(f"Log10 RPM ({aligner})" if log_scale else f"RPM ({aligner})")
    p.set_xlabel(None)

    legend = p.get_legend()
    if legend:
        legend.set_title(None) 
        sns.move_legend(p, "upper right")

def get_dilution_data(alignment_df: pandas.DataFrame, reference_df: pandas.DataFrame, dilution: str):

    dilution_df = alignment_df[alignment_df["dilution"] == dilution]
    dilution_ref_df = reference_df[reference_df["dilution"] == dilution]
    
    rpm_data = []
    # Loop over id and taxid
    for ((id, taxid), data_df) in dilution_df.groupby(['id', 'taxid'], as_index=False):

        # Loop over repeat
        for (repeat, repeat_df) in data_df.groupby("repeat"):

            # Step 1: Compute RPM for scan and remap alignments
            repeat_df['rpm_scan'] = (repeat_df['scan_alignments'] / repeat_df['input_reads']) * 1e6
            repeat_df['rpm_remap'] = (repeat_df['remap_alignments'] / repeat_df['input_reads']) * 1e6

            # Step 2: Sum the 'rpm_scan' and 'rpm_remap' columns
            rpm_scan_total = repeat_df['rpm_scan'].sum()
            rpm_remap_total = repeat_df['rpm_remap'].sum()

            scan_align_total = repeat_df['scan_alignments'].sum()
            remap_align_total = repeat_df['remap_alignments'].sum()

            # Step 3: Keep one representative row for other columns
            # Since other values are the same, we can use .iloc[0] to get a representative value
            repeat_data_row = repeat_df.iloc[0][["id", "repeat", "taxid", "dilution", "input_reads"]].copy()

            # Step 4: Add the summed values
            repeat_data_row["rpm_scan_total"] = rpm_scan_total
            repeat_data_row["rpm_remap_total"] = rpm_remap_total
            repeat_data_row["scan_alignments_total"] = scan_align_total
            repeat_data_row["remap_alignments_total"] = remap_align_total

            # Convert the Series into a DataFrame (one row DataFrame)
            repeat_data_df = pandas.DataFrame([repeat_data_row])

            rpm_data.append(repeat_data_df)
    
    if len(rpm_data) == 0:
        dilution_rpm_df = dilution_ref_df[["id", "repeat", "taxid", "dilution", "input_reads"]].copy()
        dilution_rpm_df["rpm_scan_total"] = np.nan
        dilution_rpm_df["rpm_remap_total"] = np.nan
        dilution_rpm_df["scan_alignments_total"] = np.nan
        dilution_rpm_df["remap_alignments_total"] = np.nan

        dilution_data_plot = dilution_rpm_df[['id', 'repeat', 'dilution', 'taxid', 'input_reads', 'scan_alignments_total', 'remap_alignments_total', 'rpm_scan_total', 'rpm_remap_total']]
        dilution_data_full = dilution_rpm_df
    else:
        dilution_rpm_df = pandas.concat(rpm_data, axis=0, ignore_index=True)

        # Merge the two dataframes on 'id', 'repeat', 'taxid', and 'dilution'
        # We use 'outer' to retain all rows, even if some taxids are missing in alignment_df
        merged_df = pandas.merge(
            dilution_ref_df, 
            dilution_rpm_df, 
            on=['id', 'repeat', 'taxid', 'dilution'], 
            how='outer', 
            suffixes=('_ref', '_align')
        )

        # Step 2: Fill missing values in 'input_reads' from the reference dataframe for rows where taxid was missing in alignment_df
        # If 'input_reads_align' is NaN, use 'input_reads_ref'
        merged_df['input_reads'] = merged_df['input_reads_align'].combine_first(merged_df['input_reads_ref'])


        # Step 3: Fill any other missing columns with default values (e.g., RPM columns, alignment-specific data)
        # merged_df['rpm_scan_total'].fillna(0, inplace=True)
        # merged_df['rpm_remap_total'].fillna(0, inplace=True)
        
        # Step 4: Drop the extra columns from reference dataframe (if they are not needed)
        merged_df.drop(['input_reads_ref', 'input_reads_align'], axis=1, inplace=True)

        # Step 5: Optional sorting based on id, repeat, and taxid (for better readability)
        merged_df.sort_values(by=['id', 'taxid', 'repeat'], inplace=True)

        # Step 6: Optional sorting based on id, repeat, and taxid (for better readability)
        dilution_data_plot = merged_df[['id', 'repeat', 'dilution', 'taxid', 'input_reads', 'scan_alignments_total', 'remap_alignments_total', 'rpm_scan_total', 'rpm_remap_total']]
        dilution_data_full = merged_df

    return dilution_data_plot, dilution_data_full


@app.command()
def plot_lod2(
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