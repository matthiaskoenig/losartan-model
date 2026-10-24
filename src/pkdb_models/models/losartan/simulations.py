"""The simulation experiments of the losartan model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.misc import (
    DoseDependencyExperiment,
    GeneticVariantsExperiment,
    RaasExperiment,
)
from pkdb_models.models.losartan.experiments.scans.scan_parameters import (
    LosartanParameterScan,
)
from pkdb_models.models.losartan.experiments.studies import (
    FDA1995S60,
    FDA1995S67,
    Azizi1999,
    Bae2011,
    Christen1991a,
    Doig1993,
    Donzelli2014,
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
)

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [
        Azizi1999,
        Bae2011,
        # the angiotensin I and angiotensin II infusion protocols are not
        # implemented, the study is not fitted
        Christen1991a,
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
        Yasar2002a,
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
        Christen1991a,
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
        # the RAAS model alone with constant E3174 concentrations
        RaasExperiment,
    ],
    "scans": [
        LosartanParameterScan,
    ],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scans"]


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
