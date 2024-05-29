
import typer
import pandas
import warnings

from pathlib import Path

app = typer.Typer(add_completion=False)

@app.command()
def create_validation_sheet(
    fastq_dir: Path = typer.Option(
        ..., help="Directory containing Illumina FASTQ files (NextSeq)"
    ),
    ref_dir: Path = typer.Option(
        ..., help="Directory containing FASTA reference genomes for alignment"
    ),
    meta_data: Path = typer.Option(
        ..., help="Meta data table from extraction experiment"
    ),
    output: Path = typer.Option(
        "sample_sheet.csv", help="Meta data table from extraction experiment"
    ),
):
    """ Create a sample sheet for the reference alignment validation workflow """

    df = pandas.read_csv(meta_data, sep="\t", header=0)

    read_files = list(fastq_dir.glob("*.fastq.gz"))

    data = []
    for i, row in df.iterrows():
        reads = [file for file in read_files if file.name.startswith(row["sample_id"])]

        if len(reads) != 2:
            print(reads)
            print(f"Could not detect paired end reads for sample identifier: {row['sample_id']}")
            continue

        ref = ref_dir / row["reference"]
        if not (ref_dir / ref).exists():
            print(f"Failed to detect reference file {row['reference']} at: {ref}")

        data.append([row["sample_id"], reads[0], reads[1], ref])

    pandas.DataFrame(data, columns=["sample_id", "forward_path", "reverse_path", "reference_path"]).to_csv(output, sep=",", index=False, header=True)
