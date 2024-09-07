import typer

from .ercc.terminal import app as ercc_app
from .phage.terminal import app as phage_app
from .taxa.terminal import app as taxa_app
from .extraction.terminal import app as extraction_app
from .sheets.terminal import app as sheets_app
from .controls.terminal import app as controls_app


app = typer.Typer(add_completion=False)

app.add_typer(ercc_app, name="ercc")
app.add_typer(phage_app, name="phage")
app.add_typer(taxa_app, name="taxa")
app.add_typer(extraction_app, name="extraction")
app.add_typer(sheets_app, name="sheets")
app.add_typer(controls_app, name="controls")

