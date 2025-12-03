"""
Reusable functionality for multiple simulation experiments.
"""
from collections import namedtuple
from typing import Dict

import pandas as pd

from pkdb_models.models.losartan import MODEL_PATH
from pkdb_models.models.losartan.losartan_pk import calculate_losartan_pk, calculate_losartan_pd
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task


# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "losp los e3174 l158 ren anggen ang1 ang2 ald")

class LosartanSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "Time"
    label_los = "Losartan"
    label_e3174 = "E3174"
    label_l158 = "L158"
    label_total = "LOS+E3174+L158"
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

    labels: Dict[str, str] = {
        "time": "time",
        "[Cve_los]": label_los,
        "[Cve_e3174]": label_e3174,
        "[Cve_l158]": label_l158,
        "[Cve_los_e3174]": label_los_e3174,
        "[Cve_total]": label_total,

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

        "[e3174]": "e3174",  # e3174 in Vplasma

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

        "fe_e3174": "E3174 effect",  # [-] effect via exp3174

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
        "DBP": "DBP\n",
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

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_los]": unit_los,
        "[Cve_e3174]": unit_e3174,
        "[Cve_l158]": unit_l158,
        "[Cve_total]": unit_los,
        "[Cve_los_e3174]": unit_los_e3174,
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

        "fe_e3174": "dimensionless",  # [-] effect via exp3174

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

    # ----------- Genotypes CYP2C9 --------------
    # see https://www.pharmgkb.org/page/cyp2c9RefMaterials
    # LI__f_cyp2c9
    cyp2c9_allele_activity = {
        "*1": 1.0,
        "*2": 0.60,  # [Kusama2009] 0.65; [Wang2014] 0.4858
        "*3": 0.17,  # [Kusama2009] 0.14; [Maekawa2009] 161/704 = 0.23; [Wang2014] 0.1989
        "*13": 0.05,  # [Wang2014] 0.082; [Maekawa2009] 17.6/704 = 0.025
    }
    cyp2c9_activity = {
        "*1/*1": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*1"]) / 2,
        "*1/*2": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*2"]) / 2,
        "*1/*3": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*3"]) / 2,
        "*1/*13": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*13"]) / 2,
        "*2/*2": (cyp2c9_allele_activity["*2"] + cyp2c9_allele_activity["*2"]) / 2,
        "*2/*3": (cyp2c9_allele_activity["*2"] + cyp2c9_allele_activity["*3"]) / 2,
        "*3/*3": (cyp2c9_allele_activity["*3"] + cyp2c9_allele_activity["*3"]) / 2,
    }
    cyp2c9_colors = {
        "*1/*1": "black",
        "*1/*2": "tab:blue",
        "*1/*3": "tab:red",
        "*1/*13": "tab:green",
        "*2/*2": "tab:purple",
        "*2/*3": "tab:cyan",
        "*3/*3": "tab:orange",
    }

    # ----------- Genotypes ABCB1/MDR1/P-GP --------------
    # see https://www.pharmgkb.org/page/cyp2c9RefMaterials
    # GU__f_abcb1
    abcb1_allele_activity = {
        # exon 21 (changes in protein expression, no changes in mRNA)
        "c.2677 G": 1.0,  # 13.9 [0.20-56.3] / 13.9 [0.20-56.3] = 1.00 [Siegmund2002]
        "c.2677 T": 0.122,  # 1.69 [0.00-57.3] / 13.9 [0.20-56.3] = 0.122 [Siegmund2002]

        # exon 26 (changes in protein expression, no changes in mRNA)
        "c.3435 C": 1.0,  # 1275/1275 [Hoffmeyer2000]; 10.1 [0.20-56.3] / 10.1 [0.20-56.3] = 1.0 [Siegmund2002]
        "c.3435 T": 0.49,  # 627/1275 = 0.49 [Hoffmeyer2000]; 6.37 [0.00-57.3] / 10.1 [0.20-56.3] = 0.63 [Siegmund2002]
    }
    abcb1_activity = {
        "GG/CC": (2 * abcb1_allele_activity["c.2677 G"] + 2 * abcb1_allele_activity["c.3435 C"]) / 4,
        "GT/CT": (abcb1_allele_activity["c.2677 G"] + abcb1_allele_activity["c.2677 T"] + abcb1_allele_activity["c.3435 C"] + abcb1_allele_activity["c.3435 T"]) / 4,
        "TT/TT": (2 * abcb1_allele_activity["c.2677 T"] + 2 * abcb1_allele_activity["c.3435 T"]) / 4,
    }
    abcb1_colors = {
        "GG/CC": "black",
        "GT/CT": "#d6272855",
        "TT/TT": "tab:red",
    }

    # ----------- Renal map --------------
    # glomerular filtration rate (GFR); estimated (eGFR); creatinine clearance; 100 ml/min
    # KI__f_renal_function
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 69.5 / 101.0,  # 0.69
        "Moderate renal impairment": 32.5 / 101.0,  # 0.32
        "Severe renal impairment": 19.5 / 101.0,  # 0.19
        # "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        # "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    # f_cirrhosis
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.399,  # CPT A
        "Moderate cirrhosis": 0.698,  # CPT B
        "Severe cirrhosis": 0.813,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",  # CPT C
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""
        changes = {
            # Pharmacokinetic
            # 20250708_183921__4fba0/LOSARTAN_LSQ_PK
            # 'ftissue_los': Q_(0.14494266425307856, 'l/min'),  # [0.01 - 10]
            # 'Kp_los': Q_(3.262148956929227, 'dimensionless'),  # [1 - 200]

            # 'GU__LOSABS_k': Q_(0.03617722689911161, '1/min'),  # [0.0001 - 1]
            # 'GU__f_LOSEFL_k': Q_(1.020035540752619, 'dimensionless'),  # [0.1 - 10]
            # 'GU__METEXC_k': Q_(2.6964160819935065e-05, '1/min'),  # [1e-05 - 0.1]

            # 'LI__E3174EX_k': Q_(0.011200118227296746, '1/min'),  # [0.001 - 10]
            # 'LI__LOS2E3174_Vmax': Q_(0.0007251463173167296, 'mmol/min/l'),  # [1e-05 - 100]
            # 'LI__E3174L158_k': Q_(0.0011340090993695589, '1/min'),  # [1e-05 - 10]
            # 'LI__MBIEX_k': Q_(0.066315021777011, '1/min'),  # [1e-05 - 1]

            # 'KI__LOSEX_k': Q_(0.0773896174258953, '1/min'),  # [0.0001 - 1]
            # 'KI__E3174EX_k': Q_(0.028927173263331843, '1/min'),  # [0.0001 - 1]
            # 'KI__L158EX_k': Q_(0.28890834506773044, '1/min'),  # [0.0001 - 1]

            # Pharmacodynamic
            # 20250711_231400__8d0b3/LOSARTAN_LSQ_PD
            # >>> !Optimal parameter 'ALDSEC_k' within 5% of lower bound! <<<
            # 'ANGGEN2ANG1_k': Q_(0.10030014354408065, 'l/min'),  # [0.01 - 100000.0]
            # 'E50_e3174': Q_(0.00029107274011198355, 'mM'),  # [5e-07 - 0.05]
            # 'ALDSEC_k': Q_(1.0105027952685962e-06, 'mmole/min'),  # [1e-06 - 1]
            # 'BP_ald_fe': Q_(0.3119992626947359, 'dimensionless'),  # [0.1 - 0.6]
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return LosartanSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",

                # dosing
                "IVDOSE_los",
                "PODOSE_los",
                "IVDOSE_e3174",

                # venous plasma
                "[Cve_los]",
                "[Cve_e3174]",
                "[Cve_l158]",
                "[Cve_los_e3174]",
                "[Cve_total]",
                "mr_e3174_los_plasma",
                # "mr_los_e3174_plasma",

                # urine
                "Aurine_los",
                "Aurine_e3174",
                "Aurine_l158",
                "Aurine_los_e3174",
                "Aurine_total",
                "mr_e3174_los_urine",
                # "mr_los_e3174_urine",

                # feces
                "Afeces_los",
                "Afeces_e3174",
                "Afeces_l158",
                "Afeces_total",
                "mr_e3174_los_feces",
                # "mr_los_e3174_feces",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',

                'LI__f_cyp2c9',
                'GU__f_abcb1',

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
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            losp=self.Q_(461.0, "g/mole"),
            los=self.Q_(422.911, "g/mole"),
            e3174=self.Q_(436.894, "g/mole"),
            l158=self.Q_(436.894, "g/mole"),
            ren=self.Q_(45057, "g/mole"),
            anggen=self.Q_(52670, "g/mole"),
            ang1=self.Q_(1296.499, "g/mole"),
            ang2=self.Q_(1046.179, "g/mole"),
            ald=self.Q_(360.444, "g/mole"),
        )
    # 461.0/422.911

    # --- Pharmacokinetic parameters ---
    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        "cl_renal": "Renal clearance",
        "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
        "Aurine_eat": "Enalaprilat urine",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "µmole/l*hr",
        "cl": "ml/min",
        "cl_renal": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
        "Aurine_eat": "µmole",
    }

    def calculate_losartan_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pk parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_losartan_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_losartan_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs

    def calculate_losartan_pd(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pd parameters for simulations (scans)"""
       pd_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_losartan_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_losartan_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       return pd_dfs