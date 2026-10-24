"""Base class of the losartan simulation experiments."""

from collections.abc import Mapping, Sequence
from typing import ClassVar, NamedTuple

import pandas as pd
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.impairment import RENAL_COLORS_END_STAGE, RENAL_FUNCTION_GFR
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.losartan_pk import (
    calculate_losartan_pd,
    calculate_losartan_pk,
)


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    losp: Quantity  # losartan potassium, the salt of the tablets
    los: Quantity
    e3174: Quantity
    l158: Quantity
    ren: Quantity
    anggen: Quantity
    ang1: Quantity
    ang2: Quantity
    ald: Quantity


class LosartanSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the losartan simulation experiments.

    The model doses losartan as the free acid (`PODOSE_los`, `IVDOSE_los`,
    `Ri_los`): a dose of losartan potassium (all tablets and solutions) is
    converted with `losartan_dose`. The plasma concentrations are the
    observables `[Cveplasma_<sid>]`, the venous blood concentrations
    `[Cve_<sid>]` (blood to plasma ratio 0.53 losartan, 0.58 E3174). The
    groups of a study are set by changes of the model: the CYP2C9 genotype
    (`cyp2c9_activity`, `LI__f_cyp2c9`), interacting drugs
    (`cyp2c9_interactions`, `LI__f_cyp2c9`), the ABCB1 genotype
    (`abcb1_activity`, `GU__f_abcb1`), the renal function (`renal_map`,
    `KI__f_renal_function`, a creatinine clearance of X ml/min is X/101), the
    cirrhosis (`cirrhosis_map`, `f_cirrhosis`), the body weight (`BW`) and
    the baselines of the RAAS and the blood pressure (`*_ref`).
    """

    pkpd_model = MODEL
    selections = (
        "time",
        # dosing
        "IVDOSE_los",
        "PODOSE_los",
        "IVDOSE_e3174",
        "Ri_los",
        "Ri_e3174",
        # venous plasma
        "[Cveplasma_los]",
        "[Cveplasma_e3174]",
        "[Cveplasma_l158]",
        "[Cveplasma_los_e3174]",
        "[Cveplasma_total]",
        "mr_e3174_los_plasma",
        # venous blood
        "[Cve_los]",
        "[Cve_e3174]",
        # urine
        "Aurine_los",
        "Aurine_e3174",
        "Aurine_l158",
        "Aurine_los_e3174",
        "Aurine_total",
        "mr_e3174_los_urine",
        # feces
        "Afeces_los",
        "Afeces_e3174",
        "Afeces_l158",
        "Afeces_total",
        "mr_e3174_los_feces",
        # cases
        "BW",
        "KI__f_renal_function",
        "KI__GFR",
        "f_cirrhosis",
        "LI__f_cyp2c9",
        "GU__f_abcb1",
        "[e3174]",
        "[anggen]",
        "[ang1]",
        "ang1_change",
        "ang1_ratio",
        "[ang2]",
        "ang2_change",
        "ang2_ratio",
        "[ren]",
        "ren_change",
        "ren_ratio",
        "[ald]",
        "ald_change",
        "ald_ratio",
        "fe_e3174",
        "f_at1",
        "ANGGEN2ANG1",
        "ANG1ANG2",
        "ANG2DEG",
        "RENSEC",
        "RENDEG",
        "ALDSEC",
        "ALDDEG",
        "SBP",
        "SBP_change",
        "SBP_ratio",
        "DBP",
        "DBP_change",
        "DBP_ratio",
        "MAP",
    )

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "Time"
    label_los = "Losartan"
    label_e3174 = "E3174"
    label_l158 = "other metabolites"
    label_total = "LOS+E3174+metabolites"
    label_los_e3174 = f"{label_los}+{label_e3174}\n"
    label_mr_los_e3174 = f"{label_los}/{label_e3174}\n"
    label_mr_e3174_los = f"{label_e3174}/{label_los}\n"

    label_los_urine = label_los + " urine"
    label_e3174_urine = label_e3174 + " urine"
    label_l158_urine = label_l158 + " urine"
    label_total_urine = label_total + " urine"
    label_los_e3174_urine = label_los_e3174 + "urine"
    label_mr_los_e3174_urine = label_mr_los_e3174 + " urine"
    label_mr_e3174_los_urine = label_mr_e3174_los + " urine"

    label_los_feces = label_los + " feces"
    label_e3174_feces = label_e3174 + " feces"
    label_l158_feces = label_l158 + " feces"
    label_total_feces = label_total + " feces"
    label_mr_los_e3174_feces = label_mr_los_e3174 + " feces"
    label_mr_e3174_los_feces = label_mr_e3174_los + " feces"

    labels: ClassVar[Mapping[str, str]] = {
        "time": label_time,
        "[Cveplasma_los]": label_los,
        "[Cveplasma_e3174]": label_e3174,
        "[Cveplasma_l158]": label_l158,
        "[Cveplasma_los_e3174]": label_los_e3174,
        "[Cveplasma_total]": label_total,
        "[Cve_los]": label_los + " blood",
        "[Cve_e3174]": label_e3174 + " blood",
        "mr_e3174_los_plasma": label_mr_e3174_los,
        "mr_los_e3174_plasma": label_mr_los_e3174,
        "Aurine_los": label_los_urine,
        "Aurine_e3174": label_e3174_urine,
        "Aurine_l158": label_l158_urine,
        "Aurine_los_e3174": label_los_e3174_urine,
        "Aurine_total": label_total_urine,
        "mr_e3174_los_urine": label_mr_e3174_los_urine,
        "mr_los_e3174_urine": label_mr_los_e3174_urine,
        "Afeces_los": label_los_feces,
        "Afeces_e3174": label_e3174_feces,
        "Afeces_l158": label_l158_feces,
        "Afeces_total": label_total_feces,
        "mr_e3174_los_feces": label_mr_e3174_los_feces,
        "mr_los_e3174_feces": label_mr_los_e3174_feces,
        "[e3174]": label_e3174 + " (RAAS)",  # e3174 in Vplasma
        "[anggen]": "angiotensinogen",  # angiotensinogen in Vplasma
        "[ang1]": "angiotensin I",  # angiotensin I in Vplasma
        "ang1_change": "angiotensin I change\n",
        "ang1_ratio": "angiotensin I ratio\n",
        "[ang2]": "angiotensin II",  # angiotensin II in Vplasma
        "ang2_change": "angiotensin II change\n",
        "ang2_ratio": "angiotensin II ratio\n",
        "[ren]": "renin",  # renin in Vplasma
        "ren_change": "renin change\n",
        "ren_ratio": "renin ratio\n",
        "[ald]": "aldosterone",  # aldosterone in Vplasma
        "ald_change": "aldosterone change\n",
        "ald_ratio": "aldosterone ratio\n",
        "fe_e3174": "AT1 blockade",  # [-] AT1 receptors blocked by E3174
        "f_at1": "AT1 signal",  # [-] relative to the baseline
        "ANGGEN2ANG1": "ANGGEN2ANG1",  # angiotensinogen to angiotensin I (renin)
        "ANG1ANG2": "ANG1ANG2",  # angiotensin I to angiotensin II (ANG1ANG2)
        "ANG2DEG": "ANG2DEG",  # angiotensin II degradation (ANG2DEG)
        "RENSEC": "renin secretion\n(RENSEC)",  # renin secretion (RENSEC)
        "RENDEG": "renin degradation\n(RENDEG)",  # renin degradation (RENDEG)
        "ALDSEC": "aldosterone secretion\n(ALDSEC)",  # aldosterone secretion (ALDSEC)
        "ALDDEG": "aldosterone degradation\n(ALDDEG)",  # aldosterone degradation (ALDDEG)
        "SBP": "SBP",
        "SBP_change": "SBP change\n",
        "SBP_ratio": "SBP ratio\n",
        "DBP": "DBP",
        "DBP_change": "DBP change\n",
        "DBP_ratio": "DBP ratio\n",
        "MAP": "MAP",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "nM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_los = unit_metabolite
    unit_e3174 = unit_metabolite
    unit_l158 = unit_metabolite
    unit_los_e3174 = unit_metabolite
    unit_total = unit_metabolite

    unit_los_urine = unit_metabolite_urine
    unit_e3174_urine = unit_metabolite_urine
    unit_l158_urine = unit_metabolite_urine
    unit_los_e3174_urine = unit_metabolite_urine
    unit_total_urine = unit_metabolite_urine

    unit_los_feces = unit_metabolite_feces
    unit_e3174_feces = unit_metabolite_feces
    unit_l158_feces = unit_metabolite_feces
    unit_total_feces = unit_metabolite_feces

    unit_mr = "dimensionless"

    units: ClassVar[Mapping[str, str]] = {
        "time": unit_time,
        "[Cveplasma_los]": unit_los,
        "[Cveplasma_e3174]": unit_e3174,
        "[Cveplasma_l158]": unit_l158,
        "[Cveplasma_total]": unit_total,
        "[Cveplasma_los_e3174]": unit_los_e3174,
        "[Cve_los]": unit_los,
        "[Cve_e3174]": unit_e3174,
        "Aurine_los": unit_los_urine,
        "Aurine_e3174": unit_e3174_urine,
        "Aurine_l158": unit_l158_urine,
        "Aurine_los_e3174": unit_los_e3174_urine,
        "Aurine_total": unit_total_urine,
        "Afeces_los": unit_los_feces,
        "Afeces_e3174": unit_e3174_feces,
        "Afeces_l158": unit_l158_feces,
        "Afeces_total": unit_total_feces,
        "mr_e3174_los_plasma": unit_mr,
        "mr_los_e3174_plasma": unit_mr,
        "mr_e3174_los_urine": unit_mr,
        "mr_los_e3174_urine": unit_mr,
        "mr_e3174_los_feces": unit_mr,
        "mr_los_e3174_feces": unit_mr,
        "[e3174]": "nM",  # e3174 in Vplasma [0-1000 nM]
        "[anggen]": "pM",  # angiotensinogen in Vplasma
        "[ang1]": "pM",  # angiotensin I in Vplasma
        "ang1_change": "pM",
        "ang1_ratio": "dimensionless",
        "[ang2]": "pM",  # angiotensin II in Vplasma
        "ang2_change": "pM",
        "ang2_ratio": "dimensionless",
        "[ren]": "pM",  # renin in Vplasma
        "ren_change": "pM",
        "ren_ratio": "dimensionless",
        "[ald]": "pM",  # aldosterone in Vplasma
        "ald_change": "pM",
        "ald_ratio": "dimensionless",
        "fe_e3174": "dimensionless",  # [-] AT1 receptors blocked by E3174
        "f_at1": "dimensionless",  # [-] AT1 signal
        "ANGGEN2ANG1": "pmol/min",  # angiotensinogen to angiotensin I (renin)
        "ANG1ANG2": "pmol/min",  # angiotensin I to angiotensin II (ACE)
        "ANG2DEG": "pmol/min",  # angiotensin II degradation (ANG2DEG)
        "RENSEC": "pmol/min",  # renin secretion (RENSEC)
        "RENDEG": "pmol/min",  # renin degradation (RENDEG)
        "ALDSEC": "pmol/min",  # aldosterone secretion (ALDSEC)
        "ALDDEG": "pmol/min",  # aldosterone degradation (ALDDEG)
        "SBP": "mmHg",
        "SBP_change": "mmHg",
        "SBP_ratio": "dimensionless",
        "DBP": "mmHg",
        "DBP_change": "mmHg",
        "DBP_ratio": "dimensionless",
        "MAP": "mmHg",
    }

    # ----------- Renal function --------------
    #: the renal function of the classes of renal impairment relative to the
    #: healthy GFR of 101 ml/min, a creatinine clearance of X ml/min is X/101
    renal_map = RENAL_FUNCTION_GFR
    renal_colors = RENAL_COLORS_END_STAGE

    # ----------- Genotypes CYP2C9 --------------
    # see https://www.pharmgkb.org/page/cyp2c9RefMaterials
    # LI__f_cyp2c9, the activity of a diplotype is the mean of its alleles
    cyp2c9_allele_activity: ClassVar[dict[str, float]] = {
        "*1": 1.0,
        "*2": 0.60,  # [Kusama2009] 0.65; [Wang2014] 0.4858
        # in vitro 0.14-0.23 [Kusama2009, Maekawa2009, Wang2014]; in vivo the
        # metabolic ratio E3174/losartan of *3/*3 is 0.05 of *1/*1 (losartan
        # AUC x2.4, E3174 AUC x0.12 [Bae2011]; E3174 AUC x0.07 [Yasar2002a]),
        # of *1/*3 0.40-0.50 [Han2009a, Sekino2003]: the model's *3/*3 with
        # 0.05 gives losartan AUC x2.2, E3174 AUC x0.11
        "*3": 0.05,
        "*13": 0.05,  # [Wang2014] 0.082; [Maekawa2009] 17.6/704 = 0.025
    }
    cyp2c9_activity: ClassVar[dict[str, float]] = {
        "*1/*1": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*1"]) / 2,
        "*1/*2": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*2"]) / 2,
        "*1/*3": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*3"]) / 2,
        "*1/*13": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*13"]) / 2,
        "*2/*2": (cyp2c9_allele_activity["*2"] + cyp2c9_allele_activity["*2"]) / 2,
        "*2/*3": (cyp2c9_allele_activity["*2"] + cyp2c9_allele_activity["*3"]) / 2,
        "*3/*3": (cyp2c9_allele_activity["*3"] + cyp2c9_allele_activity["*3"]) / 2,
    }
    cyp2c9_colors: ClassVar[dict[str, str]] = {
        "*1/*1": "black",
        "*1/*2": "tab:blue",
        "*1/*3": "tab:red",
        "*1/*13": "tab:green",
        "*2/*2": "tab:purple",
        "*2/*3": "tab:cyan",
        "*3/*3": "tab:orange",
    }

    # ----------- Interactions on CYP2C9 --------------
    #: the activity of CYP2C9 (`LI__f_cyp2c9`) with an interacting drug, from
    #: the observed change of the metabolic ratio E3174/losartan (AUC), which
    #: the model reproduces with the activity; documented, not fitted, the
    #: arms are validation data
    cyp2c9_interactions: ClassVar[dict[str, float]] = {
        # bucolome 300 mg/day: losartan AUC x1.67, E3174 AUC x0.165
        # (ratio 0.10) [Kobayashi2008]
        "bucolome": 0.10,
        # phenytoin 4 mg/kg/day: E3174 AUC -63 %, losartan AUC +17 %
        # (ratio 0.32) [Fischer2002]
        "phenytoin": 0.32,
        # silymarin 140 mg tid: metabolic ratio 2.63/5.79 = 0.45 in *1/*1
        # [Han2009a]
        "silymarin": 0.45,
        # fermented red ginseng: metabolic ratio x0.77 [Kim2016]
        "ginseng": 0.77,
    }

    # ----------- Genotypes ABCB1/MDR1/P-GP --------------
    # see https://www.pharmgkb.org/page/cyp2c9RefMaterials
    # GU__f_abcb1
    abcb1_allele_activity: ClassVar[dict[str, float]] = {
        # exon 21 (changes in protein expression, no changes in mRNA)
        "c.2677 G": 1.0,  # 13.9 [0.20-56.3] / 13.9 [0.20-56.3] = 1.00 [Siegmund2002]
        "c.2677 T": 0.122,  # 1.69 [0.00-57.3] / 13.9 [0.20-56.3] = 0.122 [Siegmund2002]
        # exon 26 (changes in protein expression, no changes in mRNA)
        "c.3435 C": 1.0,  # 1275/1275 [Hoffmeyer2000]; 10.1 [0.20-56.3] / 10.1 [0.20-56.3] = 1.0 [Siegmund2002]
        "c.3435 T": 0.49,  # 627/1275 = 0.49 [Hoffmeyer2000]; 6.37 [0.00-57.3] / 10.1 [0.20-56.3] = 0.63 [Siegmund2002]
    }
    abcb1_activity: ClassVar[dict[str, float]] = {
        "GG/CC": (
            2 * abcb1_allele_activity["c.2677 G"]
            + 2 * abcb1_allele_activity["c.3435 C"]
        )
        / 4,
        "GT/CT": (
            abcb1_allele_activity["c.2677 G"]
            + abcb1_allele_activity["c.2677 T"]
            + abcb1_allele_activity["c.3435 C"]
            + abcb1_allele_activity["c.3435 T"]
        )
        / 4,
        "TT/TT": (
            2 * abcb1_allele_activity["c.2677 T"]
            + 2 * abcb1_allele_activity["c.3435 T"]
        )
        / 4,
    }
    abcb1_colors: ClassVar[dict[str, str]] = {
        "GG/CC": "black",
        "GT/CT": "#d6272855",
        "TT/TT": "tab:red",
    }

    # molecular weights of losartan (potassium salt), its metabolites and the RAAS
    Mr: ClassVar[MolecularWeights] = MolecularWeights(
        losp=Q(461.0, "g/mole"),
        los=Q(422.911, "g/mole"),
        e3174=Q(436.894, "g/mole"),
        # the hydroxylated losartan L-158,795/796 (C22H23ClN6O2)
        l158=Q(438.9, "g/mole"),
        ren=Q(45057, "g/mole"),
        anggen=Q(52670, "g/mole"),
        ang1=Q(1296.499, "g/mole"),
        ang2=Q(1046.179, "g/mole"),
        ald=Q(360.444, "g/mole"),
    )

    def losartan_dose(self, dose: Quantity) -> Quantity:
        """Get the dose of losartan (free acid) of a dose of losartan potassium.

        The tablets and solutions (Cozaar, DuP 753, MK-954) contain losartan
        potassium (461.0 g/mole), the model doses the free acid (422.9
        g/mole): 50 mg losartan potassium are 45.9 mg losartan.

        Args:
            dose: dose of losartan potassium (a mass).

        Returns:
            The dose of losartan as free acid in mg.
        """
        return (dose * self.Mr.los / self.Mr.losp).to("mg")

    def calculate_losartan_pk(
        self, scans: Sequence[str] | None = None
    ) -> dict[str, pd.DataFrame]:
        """Calculate the pharmacokinetic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.

        Returns:
            The parameters of every scan by simulation key, per draw of the
            fitted parameters in a run over them (`uncertainty`).
        """
        across = None if self.uncertainty is None else self.uncertainty.id
        keys = scans if scans else self._simulations
        return {
            key: calculate_losartan_pk(
                experiment=self, res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }

    def calculate_losartan_pd(
        self, scans: Sequence[str] | None = None
    ) -> dict[str, dict[str, pd.DataFrame]]:
        """Calculate the pharmacodynamic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.

        Returns:
            The parameters of every scan by simulation key and readout, per
            draw of the fitted parameters in a run over them (`uncertainty`).
        """
        across = None if self.uncertainty is None else self.uncertainty.id
        keys = scans if scans else self._simulations
        return {
            key: calculate_losartan_pd(res=self.results[f"task_{key}"], across=across)
            for key in keys
        }
