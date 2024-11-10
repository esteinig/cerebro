from setuptools import setup, find_packages

setup(
    name="cerebro-utils",
    url="https://github.com/esteinig/cerebro",
    author="Eike J. Steinig",
    author_email="eike.steinig@unimelb.edu.au",
    packages=find_packages(),
    include_package_data=True, 
    package_data={
        'utils.assets': ['ercc.tsv']
    },
    install_requires=[
        "typer",
        "pandas",
        "seaborn",
        "scipy"
    ],
    entry_points="""
        [console_scripts]
        cerebro-utils=utils.terminal:app
    """,
    version="0.7.0",
    license="MIT",
    description="Python utilities for plotting Cerebro experiment results",
)