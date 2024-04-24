import typer

from .ercc.terminal import app as ercc_app
from .phage.terminal import app as phage_app
from .taxa.terminal import app as taxa_app
from .extraction.terminal import app as extraction_app

app = typer.Typer(add_completion=False)

app.add_typer(ercc_app, name="ercc")
app.add_typer(phage_app, name="phage")
app.add_typer(taxa_app, name="taxa")
app.add_typer(extraction_app, name="extraction")


