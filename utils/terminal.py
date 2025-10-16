import typer

from .lod.terminal import app as lod_app
from .ercc.terminal import app as ercc_app
from .phage.terminal import app as phage_app
from .taxa.terminal import app as taxa_app
from .sheets.terminal import app as sheets_app
from .manuscript.terminal import app as manuscript_app


app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)

app.add_typer(lod_app, name="lod")
app.add_typer(ercc_app, name="ercc")
app.add_typer(phage_app, name="phage")
app.add_typer(taxa_app, name="taxa")
app.add_typer(sheets_app, name="sheets")
app.add_typer(manuscript_app, name="manuscript")

