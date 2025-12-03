"""Run all simulation experiments."""
import shutil
from pathlib import Path
from typing import List

from pymetadata.console import console

from pkdb_models.models import losartan
from pkdb_models.models.losartan.helpers import run_experiments
from pkdb_models.models.losartan.experiments.studies import *
from pkdb_models.models.losartan.experiments.misc import *
from pkdb_models.models.losartan.experiments.scans.scan_parameters import LosartanParameterScan

from sbmlutils import log

from sbmlsim.plot import Figure
Figure.legend_fontsize = 8
Figure.fig_dpi = 300


logger = log.get_logger(__name__)

EXPERIMENTS = {
    "studies": [
        Azizi1999,
        Bae2011,
        Doig1993,
        Donzelli2014,
        FDA1995S60,
        FDA1995S67,
        Fischer2002,
        Goldberg1995,
        Goldberg1995a,
        Han2009a,
        Huang2021,
        Kim2016,
        Kobayashi2008,
        Lee2003b,
        Li2009,
        Lo1995,
        Munafo1992,
        Oh2012,
        Ohtawa1993,
        Puris2019,
        Sekino2003,
        Shin2020,
        Sica1995,
        Tanaka2014,
        Yasar2002a,
        # Christen1991a,  # excluded: ang1 and ang2 protocols
    ],
    "abcb1": [
        Shin2020,
    ],
    "cyp2c9": [
        Bae2011,
        Han2009a,
        Huang2021,
        Lee2003b,
        Li2009,
        Sekino2003,  # PK & PD
        Yasar2002a
    ],
    "dose_dependency": [
        Doig1993,  # PD
        Goldberg1995a,  # PK & PD
        Munafo1992,  # PK & PD
        Ohtawa1993,  # PK & PD
    ],
    "hepatic_impairment": [
        FDA1995S67,
    ],
    "renal_impairment": [
        Sica1995,
    ],
    "pharmacodynamic": [
        Azizi1999,
        # Christen1991a,  # FIXME: ang1 and ang2 protocols
        Doig1993,
        Goldberg1995,
        Goldberg1995a,
        Ohtawa1993,
        Munafo1992,
        Sekino2003,
    ],
    "misc": [
        DoseDependencyExperiment,
        GeneticVariantsExperiment,
    ],
    "scan": [
        LosartanParameterScan,
    ]

}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scan"]

def run_simulation_experiments(
    selected: str = None,
    experiment_classes: List = None,
    output_dir: Path = None
) -> None:
    """Run losartan simulation experiments."""

    Figure.fig_dpi = 600
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = losartan.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = losartan.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """Run experiments."""

    # selected = "all"
    # selected = "studies"
    # selected = "pharmacodynamic"
    # selected = "cyp2c9"
    # selected = "abcb1"

    run_simulation_experiments(selected="studies")