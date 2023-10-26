import typer
import warnings

from pathlib import Path
from ..utils import read_qc_table, YESTERDAY_MEDIUM
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

app = typer.Typer(add_completion=False)

@app.command()
def plot_spikes(
    qc_table: Path = typer.Option(
        ..., help="Quality control table"
    ),
):
    
    """
    Plot the Detroit cell spike-in experiment quality control results
    """

    df = read_qc_table(qc_table, remove_ntc=False, min_ercc_constructs=None)

    fig1, axes1 = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
    
    # Accuracy of experiment
        
    ax1 = axes1[0]
    sns.barplot(x="sample_group", y="phage_reads", hue="phage_spike", data=df, ax=ax1, palette=YESTERDAY_MEDIUM, log=True)
    ax1.set_title(f"\nT4 Reads")
    ax1.set_xlabel("\n")
    ax1.set_ylabel("Reads aligned")
    ax1.set_ylim(10e-1, None)
    legend = ax1.get_legend()
    legend.set_title(None) 
    sns.move_legend(ax1, "upper left")
    sns.despine()
    ax1.grid(False)

    ax2 = axes1[1]
    sns.barplot(x="sample_group", y="phage_coverage_percent", hue="phage_spike", data=df, ax=ax2, palette=YESTERDAY_MEDIUM)
    ax2.set_title(f"\nT4 Coverage")
    ax2.set_xlabel("\n")
    ax2.set_ylabel("Coverage (%)")
    legend = ax2.get_legend()
    legend.set_title(None) 
    sns.move_legend(ax2, "upper left")
    sns.despine()
    ax2.grid(False)


    fig1_suptitle = "Phage spike-in: experiment recovery "
    fig1.suptitle(fig1_suptitle, fontsize=16, fontweight="bold")
    fig1.savefig("phage_recovery_experiment.png", dpi=300,  transparent=False)



@app.command()
def plot_subsample_spikes(
    qc_table: Path = typer.Option(
        ..., help="Quality control table"
    ),
):

    df = read_qc_table(qc_table, remove_ntc=False, min_ercc_constructs=None)

    fig1, axes1 = plt.subplots(nrows=2, ncols=3, figsize=(18,12))


    for (i, (sample_group, group_df)) in enumerate(df.groupby("sample_group")):
        
        ax1 = axes1[0][i]
        sns.stripplot(x="sampled", y="phage_reads", hue="phage_spike", data=group_df, ax=ax1, palette=YESTERDAY_MEDIUM, size=8, alpha=0.7)
        ax1.invert_xaxis()
        ax1.set_yscale('log')
        ax1.set_ylim(10e-1, 10e5)
        ax1.set_title(f"\n{sample_group}")
        ax1.set_xlabel("\n")
        ax1.set_ylabel("Reads aligned")
        legend = ax1.get_legend()
        legend.set_title(None) 
        # sns.move_legend(ax1, "upper left")
        sns.despine()
        ax1.grid(False)
        
        ax2 = axes1[1][i]
        sns.stripplot(x="sampled", y="phage_coverage_percent", hue="phage_spike", data=group_df, ax=ax2, palette=YESTERDAY_MEDIUM, size=8, alpha=0.7)
        ax2.invert_xaxis()
        ax2.set_title(f"\n{sample_group}")
        ax2.set_xlabel("\n")
        ax2.set_ylabel("Coverage (%)")
        legend = ax2.get_legend()
        legend.set_title(None) 
        # sns.move_legend(ax1, "upper left")
        sns.despine()
        ax2.grid(False)

    # ax3 = axes3[panel_indices[i][0], panel_indices[i][1]]
    # sns.barplot(x="input_mass", y="biomass_error", hue="nucleic_acid", data=group_df, ax=ax3, palette=YESTERDAY_MEDIUM)
    # ax3.set_title(f"\nInput: {ercc_input_mass} pg")
    # ax3.set_xlabel("Detroit cell input biomass (pg)")
    # ax3.set_ylabel("Delta estimated - cell input biomass (pg)\n")
    # legend = ax3.get_legend()
    # legend.set_title(None) 
    # sns.despine()
    # ax3.grid(False)

    fig1_suptitle = "Phage spike-in: subsample recovery "
    fig1.suptitle(fig1_suptitle, fontsize=16, fontweight="bold")
    fig1.savefig("phage_recovery_subsample.png", dpi=300,  transparent=False)

    # fig2_suptitle = biomass_label + f"\n{pipeline_variant}"
    # fig2.suptitle(fig2_suptitle, fontsize=16, fontweight="bold")

    # fig3_suptitle = biomass_label+" ERROR" + f"\n{pipeline_variant}"
    # fig3.suptitle(fig3_suptitle, fontsize=16, fontweight="bold")


    # fig2.savefig("biomass_estimated.png", dpi=300,  transparent=False)
    # fig3.savefig("biomass_error.png", dpi=300,  transparent=False)

    # fig4, axes4 = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

    # sns.stripplot(x="input_mass", y="ercc_percent", hue="ercc_input_mass", data=df, palette=YESTERDAY_MEDIUM, ax=axes4, size=16)
    # axes4.set_xlabel("\nDetroit cell input biomass (pg)")
    # axes4.set_ylabel("Proportion of total reads (%)\n")
    # legend = axes4.get_legend()
    # legend.set_title("Input (pg)") 
    # sns.despine()
    # axes4.grid(False)

    # fig4_suptitle = f"ERCC/EDCC Proportions\n{pipeline_variant}"
    # fig4.suptitle(fig4_suptitle, fontsize=16, fontweight="bold")
    # fig4.savefig("constructs_proportions.png", dpi=300,  transparent=False)

    plt.close()




@app.command()
def plot_spike_host(
    qc_table: Path = typer.Option(
        ..., help="Quality control table"
    ),
    plot: Path = typer.Option(
        "phage_host_biomass_recovery_estimates.png", help="Output plot path"
    ),
    sep: str = typer.Option(
        ",", help="Quality control table delimiter"
    ),
):
    
    """
    Plot the T4 Phage x Detroit cell spike-in experiment quality control results
    """

    df = read_qc_table(qc_table, remove_ntc=False, min_ercc_constructs=50, sep=sep)

    num_host_spikes =  len(df["host_spike"].unique())
    
    fig1, axes = plt.subplots(nrows=num_host_spikes, ncols=2, figsize=(12,6*num_host_spikes))

    for (i, (host_spike, host_spike_df)) in enumerate(df.groupby("host_spike")):
        
        host_spike_df["phage_spike"] = host_spike_df["phage_spike"].astype(str)
        
        # Accuracy of experiment
            
        ax1 = axes[i][0]
        sns.barplot(x="sample_group", y="phage_reads", hue="phage_spike", data=host_spike_df, ax=ax1, palette=YESTERDAY_MEDIUM, log=host_spike != 0)
        sns.stripplot(x="sample_group", y="phage_reads", hue="phage_spike", data=host_spike_df, ax=ax1, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)
        ax1.set_title(f"\nT4 reads aligned ({host_spike} ng)" if host_spike != 0 else "\nT4 reads aligned (controls, no phage, no host)")
        ax1.set_xlabel("\n")
        ax1.set_ylabel("Reads aligned")
        ax1.set_ylim(10e-1 if host_spike != 0 else 0, None if host_spike != 0 else 10)
        legend = ax1.get_legend().set_visible(False)
        # legend.set_title(None) 
        # sns.move_legend(ax1, "upper left")
        sns.despine()
        ax1.grid(False)

        ax2 = axes[i][1]
        sns.barplot(x="sample_group", y="phage_coverage_percent", hue="phage_spike", data=host_spike_df, ax=ax2, palette=YESTERDAY_MEDIUM)

        sns.stripplot(x="sample_group", y="phage_coverage_percent", hue="phage_spike", data=host_spike_df, ax=ax2, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)

        ax2.set_title(f"\nT4 genome coverage ({host_spike} ng)" if host_spike != 0 else "\nT4 genome coverage (controls, no phage, no host)")
        ax2.set_xlabel("\n")
        ax2.set_ylabel("Coverage (%)")
        
        ax2.set_ylim(0, 110)
        legend = ax2.get_legend().set_visible(False)

        # legend.set_title(None) 
        # sns.move_legend(ax2, "upper left")
        sns.despine()
        ax2.grid(False)
    
    handles, labels = axes[1, 0].get_legend_handles_labels()

    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, title="Phage spike-in")
    plt.subplots_adjust(bottom=0.2)

    fig1.savefig(plot, dpi=300,  transparent=False)


@app.command()
def plot_host_biomass(
    qc_table: Path = typer.Option(
        ..., help="Quality control table"
    ),
    plot: Path = typer.Option(
        "phage_host_estimates.png", help="Output plot path"
    ),
    sep: str = typer.Option(
        ",", help="Quality control table delimiter"
    ),
    ylim: bool = typer.Option(
        False, help="Force limit of y-axis on biomass plots to 2000"
    ),
    adjust_csf: bool = typer.Option(
        False, help="Adjust the CSF estimates using negative controls"
    ),
):
    
    """
    Plot the T4 Phage x Detroit cell spike-in experiment spike-in host biomass estimates (EDCC)
    """

    df = read_qc_table(qc_table, remove_ntc=False, min_ercc_constructs=50, sep=sep)

    num_phage_spike =  len(df["phage_spike"].dropna().unique())

    # For CSF-1 and CSF-1 there is host biomass already present 
    # this is evident in the 208 and 209 samples CSF-1 CSF-2 no phage no host
    # and the replicate average should be substracted from any other biomass
    # estimate to assess the spike-in host biomass

    if adjust_csf:
        df_csf1_control = df.loc[(df["sample_group"] == "CSF-1") & (df["phage_spike"].isna()) & (df["host_spike"] == 0)]
        csf1_total_biomass_average = df_csf1_control["total_biomass_deduplicated"].mean()
        csf1_host_biomass_average = df_csf1_control["host_biomass"].mean()
        
        df_csf2_control = df.loc[(df["sample_group"] == "CSF-2") & (df["phage_spike"].isna()) & (df["host_spike"] == 0)]
        csf2_total_biomass_average = df_csf2_control["total_biomass_deduplicated"].mean()
        csf2_host_biomass_average = df_csf2_control["host_biomass"].mean()

    fig1, axes = plt.subplots(nrows=num_phage_spike, ncols=2, figsize=(12,6*num_phage_spike))

    for (i, (phage_spike, phage_spike_df)) in enumerate(df.groupby("phage_spike")):
        
        phage_spike_df["host_spike"] = (phage_spike_df["host_spike"]*1000).astype(str)
        phage_spike_df["host_spike"] = [f"{conc} pg" for conc in phage_spike_df["host_spike"]]
        
        if adjust_csf:
            phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-1", "total_biomass_deduplicated"] = phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-1", "total_biomass_deduplicated"]-csf1_total_biomass_average
            phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-2", "total_biomass_deduplicated"] = phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-2", "total_biomass_deduplicated"]-csf2_total_biomass_average

            phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-1", "host_biomass"] = phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-1", "host_biomass"]-csf1_host_biomass_average
            phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-2", "host_biomass"] = phage_spike_df.loc[phage_spike_df["sample_group"] == "CSF-2", "host_biomass"]-csf2_host_biomass_average

        # Accuracy of experiment

        ax2 = axes[i][0]
        sns.barplot(x="sample_group", y="total_biomass_deduplicated", hue="host_spike", data=phage_spike_df, ax=ax2, palette=YESTERDAY_MEDIUM, log=False)
        sns.stripplot(x="sample_group", y="total_biomass_deduplicated", hue="host_spike", data=phage_spike_df, ax=ax2, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)
        ax2.set_title(f"\nTotal biomass (phage: {phage_spike})" if phage_spike != 0 else "\nTotal biomass (controls, no phage, no host)")
        ax2.set_xlabel("\n")
        ax2.set_ylabel("Total biomass (pg)")
        ax2.set_ylim(0, 2000 if ylim else None)
        legend = ax2.get_legend().set_visible(False)
        # legend.set_title(None) 
        # sns.move_

        ax3 = axes[i][1]
        sns.barplot(x="sample_group", y="host_biomass", hue="host_spike", data=phage_spike_df, ax=ax3, palette=YESTERDAY_MEDIUM, log=False)
        sns.stripplot(x="sample_group", y="host_biomass", hue="host_spike", data=phage_spike_df, ax=ax3, palette=YESTERDAY_MEDIUM, dodge=True, edgecolor="black", linewidth=2)
        ax3.set_title(f"\nHost biomass (phage: {phage_spike})" if phage_spike != 0 else "\nHost biomass (controls, no phage, no host)")
        ax3.set_xlabel("\n")
        ax3.set_ylabel("Host biomass (pg)")
        ax3.set_ylim(0, 2000 if ylim else None)
        legend = ax3.get_legend().set_visible(False)

        # legend.set_title(None) 
        # sns.move_legend(ax1, "upper left")
        sns.despine()
        ax3.grid(False)
    
    handles, labels = axes[1, 0].get_legend_handles_labels()

    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, title="Host spike-in")
    plt.subplots_adjust(bottom=0.2)
    fig1.savefig(plot, dpi=300,  transparent=False)
