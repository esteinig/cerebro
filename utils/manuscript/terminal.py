#!/usr/bin/env python3
import sys
from pathlib import Path
from typing import Optional, List, Union

import typer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scikit_posthocs as sp
import scipy.stats as st

app = typer.Typer(add_completion=False)

ParamsOrder = ["32b", "14b", "8b", "4b"]
QuantRows = ["q8", "q4"]
Pairs = {
    "sens_spec": ("sensitivity", "specificity"),
    "ppv_npv": ("ppv", "npv"),
}

def _coerce_and_clean(df: pd.DataFrame) -> pd.DataFrame:
    req = {"params", "quant", "clinical", "group", "metric", "value"}
    missing = req - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    df["params"] = df["params"].astype(str)
    df["quant"] = df["quant"].astype(str).str.lower()
    df["clinical"] = df["clinical"].astype(str).str.lower()
    df["group"] = df["group"].astype(str).str.lower()
    df["metric"] = df["metric"].astype(str).str.lower()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])
    return df

def _order_present(cats: List[str], order: List[str]) -> List[str]:
    present = set(cats)
    ordered = [c for c in order if c in present]
    ordered += [c for c in cats if c not in ordered]
    return ordered

def _parse_palette(palette: str, clinical_levels: List[str]) -> Union[str, List[str], dict]:
    """
    Accepts seaborn/matplotlib palette name, or comma-separated colors.
    Returns something valid for seaborn 'palette='.
    If two colors provided, map to ['clinical','noclinical'] by level order.
    """
    pal = palette.strip()
    if "," in pal:
        cols = [c.strip() for c in pal.split(",") if c.strip()]
        if len(cols) == 1:
            return cols
        # map by order of levels
        if len(cols) >= len(clinical_levels):
            return {lvl: cols[i] for i, lvl in enumerate(clinical_levels)}
        return cols
    return pal

def _plot_panel(
    ax,
    data: pd.DataFrame,
    x_order: List[str],
    overlay: str,
    title: str,
    pal_rep,
    pal_con,
    include_consensus: bool,
    hue_order,
    legend: bool,
):
    d_rep = data[data["group"] == "replicate"]
    d_con  = data[data["group"] == "consensus"]

    # Replicate stripplot (positions established by hue+dodge)
    sns.stripplot(
        data=d_rep, x="params", y="value",
        hue="clinical", hue_order=hue_order, order=x_order, dodge=True,
        alpha=0.9, edgecolor="black", linewidth=1.3, size=6, ax=ax,
        palette=pal_rep,
    )

    # Overlays from replicates only
    if overlay == "bar":
        sns.barplot(
            data=d_rep, x="params", y="value",
            hue="clinical", hue_order=hue_order, order=x_order, dodge=True,
            errorbar=None, alpha=0.35, ax=ax, palette=pal_rep,
        )
    elif overlay == "violin":
        sns.violinplot(
            data=d_rep, x="params", y="value",
            hue="clinical", hue_order=hue_order, order=x_order, dodge=True,
            cut=0, inner=None, alpha=0.25, ax=ax, palette=pal_rep,
        )

    # Consensus dots: SAME hue+dodge → land in each condition’s lane,
    # but colored by a separate palette and slightly larger
    if include_consensus and not d_con.empty:
        sns.stripplot(
            data=d_con, x="params", y="value",
            hue="clinical", hue_order=hue_order, order=x_order, dodge=True,
            alpha=0.9, edgecolor="black", linewidth=1.3, size=6, zorder=10,
            ax=ax, palette=pal_con,
        )
        
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(None)
    ax.set_ylabel("Performance (%)", fontsize=12)
    ax.tick_params(axis="both", which="major", labelsize=12, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2, width=0.6)
    ax.set_ylim(50, 105)

    if not legend:
        ax.legend().remove()
    else:
        handles, labels = ax.get_legend_handles_labels()
        uniq = {}
        for h, l in zip(handles, labels):
            if l in hue_order and l not in uniq:
                uniq[l] = h
        if uniq:
            ax.legend(uniq.values(), uniq.keys(), title="Condition", loc="best")
        else:
            ax.legend().remove()

    # Thinner axis lines
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_linewidth(0.8)
    
    for spine in ["top", "right"]:
        ax.spines[spine].set_linewidth(0)


def _plot_metric_grouped(
    ax,
    data: pd.DataFrame,
    x_order: List[str],
    overlay: str,
    title: str,
    pal_rep,
    pal_con,
    include_consensus: bool,
    legend: bool,
    clinical: bool
):
    """Grouped bars across params with hue=metric (sensitivity,specificity)."""
    
    hue_order = ["sensitivity", "specificity"]
    
    d_rep = data[data["group"] == "replicate"]
    d_con = data[data["group"] == "consensus"]


def _plot_metric_grouped(
    ax, data, x_order, overlay, title, pal_rep, pal_con, include_consensus, legend, clinical
):
    hue_order = ["sensitivity", "specificity"]
    d_rep = data[data["group"] == "replicate"]
    d_con = data[data["group"] == "consensus"]

    # choose which clinical slice to show
    target = "clinical" if clinical else "noclinical"
    d_rep = d_rep[d_rep["clinical"] == target]
    d_con  = d_con[d_con["clinical"]  == target]

    # points
    sns.stripplot(
        data=d_rep, x="params", y="value",
        hue="metric", hue_order=hue_order, order=x_order, dodge=True,
        alpha=0.9, edgecolor="black", linewidth=1.3, size=6, ax=ax,
        palette=pal_rep,
    )

    # overlays
    if overlay == "bar":
        sns.barplot(
            data=d_rep, x="params", y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            errorbar=None, alpha=0.35, ax=ax, palette=pal_rep,
        )
    elif overlay == "violin":
        sns.violinplot(
            data=d_rep, x="params", y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            cut=0, inner=None, alpha=0.25, ax=ax, palette=pal_rep,
        )

    if include_consensus and not d_con.empty:
        sns.stripplot(
            data=d_con, x="params", y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            alpha=0.9, edgecolor="black", linewidth=1.3, size=6, zorder=10,
            ax=ax, palette=pal_con,
        )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel(None)
    ax.set_ylabel("Performance (%)", fontsize=12)
    ax.tick_params(axis="both", which="major", labelsize=12, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2, width=0.6)
    ax.set_ylim(50, 105)

    if not legend:
        ax.legend().remove()
    else:
        ax.legend(title="Metric", loc="best")

    for spine in ["bottom", "left"]:
        ax.spines[spine].set_linewidth(0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_linewidth(0)
    

@app.command()
def plot_gpt(
    input_tsv: Path = typer.Argument(..., exists=True, readable=True, help="Replicate summary TSV"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Path to save figure (png/pdf/svg)"),
    overlay: str = typer.Option("none", "--overlay", "-O", case_sensitive=False, help="none | bar | violin"),
    pair: str = typer.Option("sens_spec", "--pair", "-p", case_sensitive=False, help="sens_spec | ppv_npv"),
    include_consensus: bool = typer.Option(False, "--include-consensus", help="Overlay consensus dots"),
    consensus_palette: str = typer.Option(
        "dark", "--consensus-palette",
        help="Palette for consensus dots. Name or comma colors for clinical,noclinical."
    ),
    palette: str = typer.Option(
        "deep", "--palette", "-P",
        help="Seaborn palette name or comma-separated colors for clinical,noclinical (e.g., 'C0,C1' or '#1f77b4,#ff7f0e')"
    ),
    height: float = typer.Option(4.0, "--height", help="Row height in inches"),
    width: float = typer.Option(10.0, "--width", help="Figure width in inches"),
    tight: bool = typer.Option(True, "--tight/--no-tight", help="Use tight_layout"),
    dpi: int = typer.Option(200, "--dpi", help="DPI when saving"),
    show: bool = typer.Option(False, "--show", help="Display window instead of just saving"),
    legend: bool = typer.Option(False, "--legend", help="Display legends"),
    clinical_comparison: bool = typer.Option(False, "--clinical-comparison", help="If set, use original clinical/noclinical layout. If omitted, switch to 1x2 quant layout with metric-grouped bars."),
    clinical: bool = typer.Option(False, "--clinical", help="If set, use clinical group to plot when not using --clinical-comparison"),
):
    """
    If --clinical-comparison: 2x2 grid (rows=quant q8,q4; cols=metric pair). Hue=clinical.
    If omitted: 1x2 grid (cols=q8,q4). Hue=metric (sensitivity,specificity) as grouped bars across params.
    """


    overlay = overlay.lower()
    if overlay not in {"none", "bar", "violin"}:
        typer.echo("overlay must be one of: none, bar, violin", err=True)
        raise typer.Exit(code=2)

    pair = pair.lower()
    if pair not in Pairs:
        typer.echo("pair must be one of: sens_spec, ppv_npv", err=True)
        raise typer.Exit(code=2)

    df = pd.read_csv(input_tsv, sep="\t")
    df = _coerce_and_clean(df)

    keep = df["group"].eq("replicate") | (include_consensus & df["group"].eq("consensus"))
    df = df[keep]

    df = df[df["quant"].isin(QuantRows)]
    if df.empty:
        typer.echo("No data for q8 or q4 after filtering.", err=True)
        raise typer.Exit(code=1)

    # Categorical ordering
    present_params = sorted(df["params"].unique(), key=lambda s: (ParamsOrder + [s]).index(s) if s in ParamsOrder else 999)
    x_order = _order_present(present_params, ParamsOrder)

    sns.set_context("talk")

    if clinical_comparison:
        # Original behavior
        m1, m2 = Pairs[pair]
        d_pair = df[df["metric"].isin([m1, m2])]
        if d_pair.empty:
            typer.echo(f"No rows for metrics {m1} or {m2}.", err=True)
            raise typer.Exit(code=1)

        clinical_levels = [lvl for lvl in ["clinical", "noclinical"]
                           if lvl in d_pair["clinical"].unique().tolist()]
        hue_order = clinical_levels
        pal_rep = _parse_palette(palette, clinical_levels)
        pal_con = _parse_palette(consensus_palette, clinical_levels)

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(width, height * 2), sharex=True, sharey=True)
        for r, q in enumerate(QuantRows):
            sub = d_pair[d_pair["quant"] == q]
            for c, metric in enumerate((m1, m2)):
                ax = axes[r, c]
                panel = sub[sub["metric"] == metric]
                if panel.empty:
                    ax.set_title(f"{metric.capitalize()} [no data]", fontsize=12)
                    ax.axis("off")
                    continue
                _plot_panel(
                    ax=ax,
                    data=panel,
                    x_order=x_order,
                    overlay=overlay,
                    title=f"{metric.capitalize()} [{q}]",
                    pal_rep=pal_rep,
                    pal_con=pal_con,
                    include_consensus=include_consensus,
                    hue_order=hue_order,
                    legend=legend
                )
    else:
        m1, m2 = Pairs["sens_spec"]
        d_pair = df[df["metric"].isin([m1, m2])]
        if d_pair.empty:
            typer.echo("No rows for metrics sensitivity or specificity.", err=True)
            raise typer.Exit(code=1)

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(width, height), sharey=True)
                
        # normalize axes to a flat array of Axes

        # Use provided palette name to get two colors
        pal_rep = _parse_palette(palette, 2)
        pal_con = _parse_palette(consensus_palette, 2)

        for i, q in enumerate(QuantRows):
            ax = axes[i]
            sub = d_pair[d_pair["quant"] == q]
            if sub.empty:
                ax.set_title(f"{q} [no data]", fontsize=12)
                ax.axis("off")
                continue

            _plot_metric_grouped(
                ax=ax,
                data=sub,
                x_order=x_order,
                overlay=overlay,
                title=f"{q}",
                pal_rep=pal_rep,
                pal_con=pal_con,
                include_consensus=include_consensus,
                legend=legend,
                clinical=clinical
            )

    sns.despine(top=False, right=False)

    if tight:
        plt.tight_layout()

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi)
        typer.echo(f"Saved: {output}")
    if show or output is None:
        plt.show()


@app.command()
def friedman_test(
    input_tsv: Path = typer.Argument(..., exists=True, readable=True, help="Replicate summary TSV"),
    quant: str = typer.Option("q8", "--quant", "-q", help="Quantization level to test (e.g. q8)"),
    metric: str = typer.Option("sensitivity", "--metric", "-m", help="Metric to test (sensitivity/specificity/ppv/npv)"),
    condition: str = typer.Option("clinical", "--condition", "-c", help="Condition to test (clinical/noclinical)"),
    posthoc: Optional[str] = typer.Option(
        None, "--posthoc", "-p",
        help="Posthoc: 'conover' for  Conover or 'wilcoxon' for pairwise Wilcoxon."
    ),
    p_adjust: str = typer.Option("holm", "--p-adjust",
                                help="P-adjust method for posthoc (e.g. holm, bonferroni, fdr_bh)"),
    group_col: str = typer.Option("params", "--group-col", help="Column defining parameter groups (e.g. params)"),
    replicate_col: str = typer.Option("replicate_id", "--replicate-col",
                                      help="Column defining replicates (if absent, auto-index used)"),
):
    """
    Perform Friedman test for differences in diagnostic performance across parameters
    at a fixed quantization and metric. Optionally run pairwise posthoc tests with correction.
    """
    df = pd.read_csv(input_tsv, sep="\t")

    df = df[df["group"] != "consensus"]
    df = df[df["clinical"] == condition]

    if "replicate_id" not in df.columns:
        # auto-assign replicate index within each (quant,clinical,group,metric)
        df = df.copy()
        df["replicate_id"] = df.groupby(["params", "quant", "clinical", "group", "metric"]).cumcount() + 1

    print(df)

    # Filter to target quant and metric
    sub = df[(df["quant"] == quant) & (df["metric"] == metric) & (df["group"] == "replicate")]
    if sub.empty:
        typer.echo(f"No replicate data for quant={quant}, metric={metric}", err=True)
        raise typer.Exit(code=1)

    sub.to_csv(input_tsv.with_name(f"posthoc_{metric}_{quant}.data.tsv"), sep="\t")
    
    sub["replicate_id"] = [f"replicate_{i}" for i in sub["replicate_id"]]

    # Wide pivot: each column = parameter (e.g. 32b, 14b, 8b, 4b)
    pivot = sub.pivot_table(index=replicate_col, columns=group_col, values="value")
    pivot = pivot.dropna(axis=0, how="any")  # remove incomplete replicate rows

    if pivot.shape[1] < 3:
        typer.echo("Need ≥3 parameter groups for Friedman test", err=True)
        raise typer.Exit(code=1)

    typer.echo(f"Testing {metric} across parameters for quant={quant}")
    stat, p = st.friedmanchisquare(*[pivot[c] for c in pivot.columns])
    typer.echo(f"Friedman χ² = {stat:.3f}, p = {p:.4g}")

    print(pivot)

    if posthoc:
        if posthoc.lower() == "wilcoxon":
            res = sp.posthoc_wilcoxon(pivot.to_numpy().T, zero_method="pratt", correction=True, p_adjust=p_adjust)
            res.columns = pivot.columns.copy()
            res.index = pivot.columns.copy()
            tag = f"wilcoxon_{p_adjust}"
        else:
            res = sp.posthoc_conover_friedman(
                pivot.to_numpy(), p_adjust=p_adjust,
            )
            tag = f"conover_friedman_{p_adjust}"

        out = input_tsv.with_name(f"posthoc_{metric}_{quant}_{tag}.tsv")
        res.to_csv(out, sep="\t")
        pd.set_option("display.float_format", lambda x: f"{x:0.4g}")
        print(res)
        typer.echo(f"Saved: {out}")
    
