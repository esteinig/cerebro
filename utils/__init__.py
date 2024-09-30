import rich.traceback

# Disable the enhanced console traceback
rich.traceback.install(show_locals=False, suppress=[__name__])
