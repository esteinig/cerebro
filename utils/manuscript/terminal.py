#!/usr/bin/env python3
import sys
from pathlib import Path
from typing import Optional, List, Union
from matplotlib.lines import Line2D

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
    "sens-spec": ("sensitivity", "specificity"),
    "ppv-npv": ("ppv", "npv"),
}


def _legend_right(ax, entries, consensus_color=None, title: Optional[str] = None, title_fontsize: int = 8):
    handles = [
        Line2D([0], [0], marker='o', linestyle='None',
               markersize=5, markerfacecolor=c,
               markeredgecolor='black', markeredgewidth=1.0)
        for _, c in entries
    ]
    labels = [lbl for lbl, _ in entries]

    if consensus_color is not None:
        handles.append(
            Line2D([0], [0], marker='o', linestyle='None',
                   markersize=5, markerfacecolor=consensus_color,
                   markeredgecolor='black', markeredgewidth=1.0)
        )
        labels.append("Consensus")

    leg = ax.legend(handles, labels,
              loc='center left', bbox_to_anchor=(1.02, 0.5),
              frameon=False, borderaxespad=0.,
              handletextpad=0.4, labelspacing=0.3,
              fontsize=8, title=title, title_fontsize=title_fontsize)
    
    if leg.get_title() is not None:
        leg.get_title().set_fontweight("bold")


def get_title(q: str, subtitle: str) -> str:

    if q == "q8":
        return f"8-bit {subtitle}"
    elif q == "q4":
        return f"4-bit {subtitle}"
    elif q == "q2":
        return f"2-bit {subtitle}"
    else:
        return f"{q} {subtitle}"

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

    if "reasoning" in df.columns:
        df["reasoning"] = df["reasoning"].astype(str).str.lower()
        
    return df

def _order_present(cats: List[str], order: List[str]) -> List[str]:
    present = set(cats)
    ordered = [c for c in order if c in present]
    ordered += [c for c in cats if c not in ordered]
    return ordered

def _parse_palette(palette: str) -> Union[str, List[str], dict]:
    """
    Accepts seaborn/matplotlib palette name, or comma-separated colors.
    Returns something valid for seaborn 'palette='.
    If two colors provided, map to ['clinical','noclinical'] by level order.
    """
    pal = palette.strip()
    if "," in pal:
        cols = [c.strip() for c in pal.split(",") if c.strip()]
        return cols
    return pal

# -------------------- p-adjust helpers --------------------
def _padjust(pvals: List[float], method: str) -> List[float]:
    m = len(pvals)
    method = (method or "none").lower()
    if method == "none" or m == 0:
        return pvals
    p = np.asarray(pvals, dtype=float)

    if method == "bonferroni":
        return np.minimum(p * m, 1.0).tolist()

    if method == "holm":
        order = np.argsort(p)
        adj = np.empty_like(p)
        rank = np.empty_like(order)
        rank[order] = np.arange(1, m + 1)
        adj_vals = (m - rank + 1) * p
        # step-down monotonicity
        adj_sorted = np.maximum.accumulate(adj_vals[order][::-1])[::-1]
        adj[order] = np.minimum(adj_sorted, 1.0)
        return adj.tolist()

    if method in {"fdr_bh", "bh", "benjamini-hochberg"}:
        order = np.argsort(p)
        ranked = np.empty_like(order)
        ranked[order] = np.arange(1, m + 1)
        q = p * m / ranked
        # enforce monotonicity from largest to smallest
        q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
        out = np.empty_like(p)
        out[order] = np.minimum(q_sorted, 1.0)
        return out.tolist()

    # fallback: none
    return pvals

def _p_to_stars(p: float, thresholds=(0.05, 0.01, 0.001)) -> str:
    if p < thresholds[2]:
        return "***"
    if p < thresholds[1]:
        return "**"
    if p < thresholds[0]:
        return "*"
    return "ns"



# --- simple bracket drawing for two-level hue within a category --------------

def _annotate_bracket(ax, x_index: int, center: float, dx: float, y: float, h: float, text: str, fontsize: int = 10, fontweight: str = "bold"):
    """Draw a bracket over two points at x positions [center-dx, center+dx]."""
    x1 = center - dx
    x2 = center + dx
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], linewidth=1.2, color="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + h + 0.2, text, ha="center", va="bottom", fontsize=fontsize, fontweight=fontweight)

def _stars_from_row(row: pd.Series, use_adj: bool, thresholds=(0.05, 0.01, 0.001)) -> str:
    p = float(row["p_adj"] if use_adj else row["p_raw"])
    return _p_to_stars(p, thresholds)

def _plot_panel(
    ax,
    data: pd.DataFrame,
    x_col: str,
    x_order: List[str],
    overlay: str,
    title: str,
    pal_rep,
    pal_con,
    include_consensus: bool,
    hue_order,
    legend: bool,
    hline_y: Optional[float] = None,
    legend_title: Optional[str] = None,
    reasoning: bool = False
):
    d_rep = data[data["group"] == "replicate"]
    d_con  = data[data["group"] == "consensus"]

    # Replicate stripplot (positions established by hue+dodge)
    sns.stripplot(
        data=d_rep, x=x_col, y="value",
        hue="reasoning" if reasoning else "clinical", hue_order=hue_order, order=x_order, dodge=True,
        alpha=0.7, edgecolor="black", linewidth=1.3, size=6, ax=ax,
        palette=pal_rep,
    )

    # Overlays from replicates only
    if overlay == "bar":
        sns.barplot(
            data=d_rep, x=x_col, y="value",
            hue="reasoning" if reasoning else "clinical", hue_order=hue_order, order=x_order, dodge=True,
            errorbar=None, alpha=0.6, ax=ax, palette=pal_rep, edgecolor="black",
            linewidth=1.5, gap=0.1
        )
    elif overlay == "violin":
        sns.violinplot(
            data=d_rep, x=x_col, y="value",
            hue="reasoning" if reasoning else "clinical", hue_order=hue_order, order=x_order, dodge=True,
            cut=0, inner=None, alpha=0.6, ax=ax, palette=pal_rep,
        )

    # Consensus dots: SAME hue+dodge → land in each condition’s lane,
    # but colored by a separate palette and slightly larger
    if include_consensus and not d_con.empty:
        sns.stripplot(
            data=d_con, x=x_col, y="value",
            hue="reasoning" if reasoning else "clinical", hue_order=hue_order, order=x_order, dodge=True,
            alpha=0.9, edgecolor="black", linewidth=1.3, size=6, zorder=10,
            ax=ax, palette=pal_con,
        )
    
    if legend_title is None and title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    else:
        ax.set_title("")
    
    ax.set_xlabel(None)
    ax.set_ylabel("Performance (%)", fontsize=12)
    ax.tick_params(axis="both", which="major", labelsize=12, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2, width=0.6)
    ax.set_ylim(50, 105)

    if hline_y is not None:
        ax.axhline(hline_y, color="darkgray", linestyle=":", linewidth=1.2, zorder=1)

    
    if legend:
        handles_, labels_ = ax.get_legend_handles_labels()
        level_to_color = {}
        hue_levels = list(hue_order)
        for h, l in zip(handles_, labels_):
            if l in hue_levels and l not in level_to_color:
                try:
                    c = h.get_facecolor()
                    if hasattr(c, '__len__') and len(c) and hasattr(c[0], '__len__'):
                        c = c[0]
                except Exception:
                    c = getattr(h, 'get_color', lambda: 'C0')()
                level_to_color[l] = c

        entries = [(lvl.capitalize(), level_to_color[lvl]) for lvl in hue_levels if lvl in level_to_color]

        consensus_color = None
        if include_consensus and not data[data["group"] == "consensus"].empty:
            consensus_color = pal_con[0] if isinstance(pal_con, (list, tuple)) else pal_con

        _legend_right(ax, entries, consensus_color, title=legend_title, title_fontsize=10)
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
    x_col: str,
    x_order: List[str],
    overlay: str,
    title: str,
    pal_rep,
    pal_con,
    include_consensus: bool,
    legend: bool,
    clinical: bool,
    hline_y: Optional[float] = None,
    legend_title: Optional[str] = None,
):
    """Grouped bars across params with hue=metric (sensitivity,specificity)."""
    
    hue_order = ["sensitivity", "specificity"]
    d_rep = data[data["group"] == "replicate"]
    d_con = data[data["group"] == "consensus"]

    # choose which clinical slice to show
    target = "clinical" if clinical else "noclinical"
    d_rep = d_rep[d_rep["clinical"] == target]
    d_con  = d_con[d_con["clinical"]  == target]

    # points
    sns.stripplot(
        data=d_rep, x=x_col, y="value",
        hue="metric", hue_order=hue_order, order=x_order, dodge=True,
        alpha=0.7, edgecolor="black", linewidth=1.3, size=6, ax=ax,
        palette=pal_rep,
    )

    # overlays
    if overlay == "bar":
        sns.barplot(
            data=d_rep, x=x_col, y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            errorbar=None, alpha=0.6, ax=ax, palette=pal_rep, edgecolor="black",
            linewidth=1.5, gap=0.1
        )
    elif overlay == "violin":
        sns.violinplot(
            data=d_rep, x=x_col, y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            cut=0, inner=None, alpha=0.6, ax=ax, palette=pal_rep,
        )

    if include_consensus and not d_con.empty:
        sns.stripplot(
            data=d_con, x=x_col, y="value",
            hue="metric", hue_order=hue_order, order=x_order, dodge=True,
            alpha=0.9, edgecolor="black", linewidth=1.3, size=6, zorder=10,
            ax=ax, palette=pal_con,
        )

    # Title control mirrors _plot_panel
    if legend_title is None and title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    else:
        ax.set_title("")

    ax.set_xlabel(None)
    ax.set_ylabel("Performance (%)", fontsize=12)
    ax.tick_params(axis="both", which="major", labelsize=12, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2, width=0.6)
    ax.set_ylim(50, 105)

    if legend:
        handles_, labels_ = ax.get_legend_handles_labels()
        level_to_color = {}
        for h, l in zip(handles_, labels_):
            if l in hue_order and l not in level_to_color:
                try:
                    c = h.get_facecolor()
                    if hasattr(c, '__len__') and len(c) and hasattr(c[0], '__len__'):
                        c = c[0]
                except Exception:
                    c = getattr(h, 'get_color', lambda: 'C0')()
                level_to_color[l] = c

        entries = [(lvl.capitalize(), level_to_color[lvl]) for lvl in hue_order if lvl in level_to_color]
        consensus_color = None
        if include_consensus and not data[data["group"] == "consensus"].empty:
            consensus_color = pal_con[0] if isinstance(pal_con, (list, tuple)) else pal_con

        _legend_right(ax, entries, consensus_color, title=legend_title, title_fontsize=10)
    else:
        ax.legend().remove()

    for spine in ["bottom", "left"]:
        ax.spines[spine].set_linewidth(0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_linewidth(0)
    
    if hline_y is not None:
        ax.axhline(hline_y, color="darkgray", linestyle=":", linewidth=1.2, zorder=1)


@app.command()
def plot_gpt(
    input_tsv: Path = typer.Argument(..., exists=True, readable=True, help="Replicate summary TSV"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Path to save figure (png/pdf/svg)"),
    overlay: str = typer.Option("none", "--overlay", "-O", case_sensitive=False, help="none | bar | violin"),
    pair: str = typer.Option("sens-spec", "--pair", "-p", case_sensitive=False, help="sens-spec | ppv-npv"),
    include_consensus: bool = typer.Option(False, "--include-consensus", help="Overlay consensus dots"),
    consensus_palette: str = typer.Option("#f78462,#f78462", "--consensus-palette", help="Palette for consensus dots. Name or comma colors for clinical,noclinical."),
    palette: str = typer.Option("#768f84,#dfbf88ff", "--palette", "-P",help="Seaborn palette name or comma-separated colors for clinical,noclinical (e.g., 'C0,C1' or '#1f77b4,#ff7f0e')"),
    height: float = typer.Option(4.0, "--height", help="Row height in inches"),
    width: float = typer.Option(10.0, "--width", help="Figure width in inches"),
    tight: bool = typer.Option(True, "--tight/--no-tight", help="Use tight_layout"),
    dpi: int = typer.Option(200, "--dpi", help="DPI when saving"),
    show: bool = typer.Option(False, "--show", help="Display window instead of just saving"),
    legend: bool = typer.Option(False, "--legend", help="Display legends"),
    subtitle: str = typer.Option("", "--subtitle", "-s", case_sensitive=False, help="Additional subtitle for each panel"),
    comparison: bool = typer.Option(False, "--comparison", help="If set, use original clinical/noclinical layout. If omitted, switch to 1x2 quant layout with metric-grouped bars."),
    combined: bool = typer.Option(False, "--combined", help="With --comparison, plot a single row where x-axis is '{params}-{q4,q8}'. Two subplots: Sensitivity | Specificity."),
    clinical: bool = typer.Option(False, "--clinical", help="If set, use clinical group to plot when not using --comparison (default: noclinical)"),
    reasoning: bool = typer.Option(False, "--reasoning", help="If set, use reasoning column if present for --comparison (default: uses clinical column)"),
    hline: Optional[float] = typer.Option(None, "--hline", help="Draw a horizontal dotted dark-gray line at this y-value on every subplot"),
    stats_tsv: Optional[Path] = typer.Option(None, "--stats-tsv",help="TSV from condition_test for adding bracket annotations in --clinical-comparison --combined"),
    stats_alpha: float = typer.Option(0.05, "--stats-alpha", help="Alpha threshold for significance stars using adjusted p-values"),
    stats_use_adjusted: bool = typer.Option(True, "--stats-use-adjusted/--stats-use-raw", help="Use adjusted p-values for stars"),
    labels_tsv: Optional[Path] = typer.Option(None, "--labels-tsv", help="TSV with columns 'id' and 'label'. If provided, non-comparison mode uses label values on the x-axis"),
):
    """
    If --comparison: original clinical/noclinical layout. With --combined, x is '{params}-{quant}'.
    If omitted (non-comparison): 1x2 grid (q8,q4). Hue=metric. If --labels-tsv is given and input has 'id',
    x-axis uses 'label' values instead of 'params'.
    """


    overlay = overlay.lower()
    if overlay not in {"none", "bar", "violin"}:
        typer.echo("overlay must be one of: none, bar, violin", err=True)
        raise typer.Exit(code=2)

    pair = pair.lower()
    if pair not in Pairs:
        typer.echo("pair must be one of: sens-spec, ppv-npv", err=True)
        raise typer.Exit(code=2)

    df = pd.read_csv(input_tsv, sep="\t")
    df = _coerce_and_clean(df)

    # Optional labels mapping
    labels_df = None
    if labels_tsv is not None:
        labels_df = pd.read_csv(labels_tsv, sep="\t")
        need = {"id", "label"}
        missing = need - set(labels_df.columns)
        if missing:
            typer.echo(f"labels-tsv missing columns: {sorted(missing)}", err=True)
            raise typer.Exit(code=2)
        if "id" not in df.columns:
            typer.echo("Input TSV has no 'id' column but --labels-tsv was provided.", err=True)
            raise typer.Exit(code=2)
        # preserve labels file order
        labels_df = labels_df.copy()
        labels_df["id"] = labels_df["id"].astype(str)
        labels_df["label"] = labels_df["label"].astype(str)
        labels_df["__label_order__"] = np.arange(len(labels_df))
        df["id"] = df["id"].astype(str)
        df = df.merge(labels_df[["id", "label", "__label_order__"]], on="id", how="left")

    keep = df["group"].eq("replicate") | (include_consensus & df["group"].eq("consensus"))
    df = df[keep]

    df = df[df["quant"].isin(QuantRows)]
    if df.empty:
        typer.echo("No data for q8 or q4 after filtering.", err=True)
        raise typer.Exit(code=1)

    # Categorical ordering
    present_params = sorted(df["params"].unique(), key=lambda s: (ParamsOrder + [s]).index(s) if s in ParamsOrder else 999)
    x_order = _order_present(present_params, ParamsOrder)


    # Categorical ordering for labels (if any)
    use_labels = labels_df is not None and df["label"].notna().any()
    if use_labels:
        # Prefer labels file order, then fall back to order of appearance
        lbls = df[["label", "__label_order__"]].dropna().drop_duplicates().sort_values("__label_order__")["label"].tolist()
        if not lbls:
            lbls = df["label"].dropna().unique().tolist()
        x_order_labels = lbls

    # Optional stats table
    stats_df = None
    if stats_tsv is not None and stats_tsv.exists():
        stats_df = pd.read_csv(stats_tsv, sep="\t")
        # expected columns: metric, paramq, n, stat, p_raw, p_adj, method, adjust
        # normalize
        for col in ("metric", "paramq"):
            if col in stats_df.columns:
                stats_df[col] = stats_df[col].astype(str).str.lower()

    sns.set_context("talk")

    if comparison:
        
        m1, m2 = Pairs[pair]
        d_pair = df[df["metric"].isin([m1, m2])]

        if d_pair.empty:
            typer.echo(f"No rows for metrics {m1} or {m2}.", err=True)
            raise typer.Exit(code=1)

        comparison_condition = ["think", "nothink"] if reasoning else ["clinical", "noclinical"] 
        comparison_column =  d_pair["reasoning"]if reasoning else d_pair["clinical"] 

        levels = [
            lvl for lvl in comparison_condition
            if lvl in comparison_column.unique().tolist()
        ]

        hue_order = levels

        print(levels)

        pal_rep = _parse_palette(palette)
        pal_con = _parse_palette(consensus_palette)

        if combined:
            
            # Build combined x-axis "{params}-{quant}" over both q8 and q4
            d_comb = d_pair.copy()
            d_comb["paramq"] = d_comb["params"].astype(str) + "-" + d_comb["quant"].astype(str)

            # Order: for each params in ParamsOrder, list param-q8 then param-q4 if present
            present_paramq = []
            for p in x_order:
                for q in QuantRows:
                    key = f"{p}-{q}"
                    if key in set(d_comb["paramq"].unique()):
                        present_paramq.append(key)

            if not present_paramq:
                typer.echo("No combined param-quant levels found.", err=True)
                raise typer.Exit(code=1)

            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(width, height), sharex=True, sharey=True)
            fig.subplots_adjust(right=0.86)

            for c, metric in enumerate((m1, m2)):
                
                ax = axes[c]
                panel = d_comb[d_comb["metric"] == metric]
                
                if panel.empty:
                    ax.set_title(f"{metric.capitalize()} [no data]", fontsize=12)
                    ax.axis("off")
                    continue

                print(panel)

                _plot_panel(
                    ax=ax,
                    data=panel,
                    x_col="paramq",
                    x_order=present_paramq,
                    overlay=overlay,
                    title=f"{metric.capitalize()} {subtitle}",
                    pal_rep=pal_rep,
                    pal_con=pal_con,
                    include_consensus=include_consensus,
                    hue_order=hue_order,
                    legend=legend,
                    hline_y=hline,
                    legend_title=metric.capitalize(),  # legend title overrides subplot title,
                    reasoning=reasoning
        
                )
                
                # --- significance brackets per category using stats_tsv ---

                if stats_df is not None:
                    stats_sub = stats_df[stats_df["metric"] == metric]
                    if not stats_sub.empty:
                        # parameters for bracket placement
                        # approximate hue dodge half-width
                        dx = 0.18
                        # compute per-category top y
                        for i, cat in enumerate(present_paramq):
                            srow = stats_sub[stats_sub["paramq"] == cat]
                            if srow.empty:
                                continue
                            pval = float(srow.iloc[0]["p_adj" if stats_use_adjusted else "p_raw"])
                            stars = _p_to_stars(pval)
                            if stats_use_adjusted and pval >= stats_alpha:
                                # still annotate with "ns" to be explicit
                                stars = "ns"
                            y_top = panel.loc[panel["paramq"] == cat, "value"].max()
                            if pd.isna(y_top):
                                continue
                            y = float(y_top) + 2.0  # padding above max
                            _annotate_bracket(ax, i, i, dx, y, h=1.0, text=stars, fontsize=10, fontweight="bold" if "*" in stars else None)


            for ax in axes:
                ax.set_xlabel(None)
        else:
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(width, height * 2), sharex=True, sharey=True)        
            fig.subplots_adjust(right=0.86)
            
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
                        x_col="params",
                        x_order=x_order,
                        overlay=overlay,
                        title=f"{metric.capitalize()} [{get_title(q, subtitle)}]",
                        pal_rep=pal_rep,
                        pal_con=pal_con,
                        include_consensus=include_consensus,
                        hue_order=hue_order,
                        legend=legend,
                        hline_y=hline,
                        legend_title=None, 
                        reasoning=False
                    )
    else:
        m1, m2 = Pairs["sens-spec"]
        d_pair = df[df["metric"].isin([m1, m2])]
        if d_pair.empty:
            typer.echo("No rows for metrics sensitivity or specificity.", err=True)
            raise typer.Exit(code=1)

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(width, height), sharey=True)
        fig.subplots_adjust(right=0.86)

        pal_rep = _parse_palette(palette)
        pal_con = _parse_palette(consensus_palette)

        # choose x mapping
        x_col = "label" if use_labels else "params"
        x_order = x_order_labels if use_labels else x_order_params

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
                x_col=x_col,
                x_order=x_order,
                overlay=overlay,
                title=get_title(q, subtitle),
                pal_rep=pal_rep,
                pal_con=pal_con,
                include_consensus=include_consensus,
                legend=legend,
                clinical=clinical,
                hline_y=hline,
                legend_title=None,
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
    
@app.command()
def condition_test(
    input_tsv: Path = typer.Argument(..., exists=True, readable=True, help="Replicate summary TSV"),
    metrics: str = typer.Option("sens-spec", "--metrics", "-m", help="Metric to test (sens-spec | ppv-npv)"),
    condition_col: str = typer.Option("clinical", "--condition-col", help="Column containing the two conditions"),
    level_a: str = typer.Option("clinical", "--level-a", help="First condition level"),
    level_b: str = typer.Option("noclinical", "--level-b", help="Second condition level"),
    p_adjust: str = typer.Option("holm", "--p-adjust", help="none | bonferroni | holm | fdr_bh"),
    replicate_col: str = typer.Option("replicate_id", "--replicate-col",
                                      help="Replicate ID column; if absent it will be generated within each cell"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Where to write the TSV"),
):
    """
    Non-parametric paired comparison between two conditions per {params}-{quant} cell.
    Uses Wilcoxon signed-rank on matched replicates.
    Produces: metric, paramq, n, stat, p_raw, p_adj, method, adjust
    """
    df_data = pd.read_csv(input_tsv, sep="\t")
    df_clean = _coerce_and_clean(df_data)
    
    metrics = metrics.lower()

    level_a = str(level_a).lower()
    level_b = str(level_b).lower()

    if metrics == "sens-spec" or metrics == "sensitivity-specificity":
        metrics_single = ["sensitivity", "specificity"]
    elif metrics == "npv-ppv":
        metrics_single = ["npv", "ppv"]
    else:
        typer.echo("Metrics allows only: sens-spec or npv-ppv.", err=True)
        raise typer.Exit(code=1)

    data = []
    for metric in metrics_single:
        df = df_clean[(df_clean["group"] == "replicate") & (df_clean["metric"] == metric)]
        if df.empty:
            typer.echo(f"No replicate rows for requested metric: {metric}", err=True)
            raise typer.Exit(code=1)

        # ensure replicate IDs exist within each (params, quant, condition, metric)
        if replicate_col not in df.columns:
            df = df.copy()
            df["replicate_id"] = df.groupby(["params", "quant", condition_col, "metric"]).cumcount() + 1
            replicate_col = "replicate_id"

        df["paramq"] = df["params"].astype(str) + "-" + df["quant"].astype(str)

        out_rows = []
        for paramq, sub in df.groupby("paramq"):
            a = sub[sub[condition_col] == level_a].set_index(replicate_col)["value"]
            b = sub[sub[condition_col] == level_b].set_index(replicate_col)["value"]
            # align paired replicates
            pairs = a.index.intersection(b.index)
            if len(pairs) < 1:
                continue
            a_vals = a.loc[pairs].to_numpy()
            b_vals = b.loc[pairs].to_numpy()
            # If all differences are zero, wilcoxon raises. Handle gracefully.
            if np.allclose(a_vals, b_vals):
                stat, p = 0.0, 1.0
            else:
                stat, p = st.wilcoxon(a_vals, b_vals, zero_method="pratt", alternative="two-sided", correction=False)
            out_rows.append({"metric": metric, "paramq": paramq, "n": len(pairs), "stat": stat, "p_raw": p})

        if not out_rows:
            typer.echo("No comparable pairs found.", err=True)
            raise typer.Exit(code=1)

        res = pd.DataFrame(out_rows)
        res["p_adj"] = _padjust(res["p_raw"].tolist(), p_adjust)
        res["method"] = "wilcoxon_signed_rank"
        res["adjust"] = p_adjust

        data.append(res)


    if output is None:
        output = input_tsv.with_name(f"condition_test_{metrics}_{p_adjust}.tsv")

    if data:
        data_out = pd.concat(data)
        data_out.to_csv(output, sep="\t", index=False)
        typer.echo(f"Saved: {output}")
    else:
        typer.echo(f"No data generated for metrics: {metrics}", err=True)
        raise typer.Exit(code=1)
