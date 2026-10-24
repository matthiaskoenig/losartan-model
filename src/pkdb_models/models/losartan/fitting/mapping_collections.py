"""Parameter fit problems for losartan.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the training data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The training data are the single doses, or the exactly encoded repeated doses,
of oral or intravenous losartan and E3174 in healthy or hypertensive subjects
with the CYP2C9 wildtype or a genotype which was not reported, without an
interacting drug. The validation data has four parts, each a different
question: the interactions (`INTERACTION_VALIDATION`, the arms with a
coadministered drug, simulated with CYP2C9 factors of the literature which are
not fitted), the doses (`DOSE_VALIDATION`, the repeated doses and the doses of
other conditions which the single doses of the training predict), the
populations (`POPULATION_VALIDATION`, the CYP2C9 and ABCB1 genotypes, obesity,
cirrhosis and renal impairment, simulated with factors which are not fitted)
and the limitations (`LIMITATION_VALIDATION`, data of a process the model does
not describe). The summaries of data which is fitted one by one are excluded
(`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import mapping_metadata, observable_filter
from pkpdsim.metadata import Coadministration, Genotype, Health
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.losartan import DATA_PATHS, LOSARTAN_PATH
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

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share, so a named mapping is named with
#: its study
type NamedMappings = dict[str, set[str]]

EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
    Azizi1999,
    Bae2011,
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
    Kim2016,  # the losartan groups of the dextromethorphan study
    Kobayashi2008,
    Lee2003b,
    Li2009,
    Lo1995,
    Munafo1992,
    Oh2012,  # the losartan groups of the caffeine study
    Ohtawa1993,
    Puris2019,  # the losartan groups of the chlorzoxazone study
    Sekino2003,
    Shin2020,
    Sica1995,
    Tanaka2014,  # the losartan groups of the caffeine study
    Yasar2002a,
]

#: health of the training data
CONTROL_HEALTH: frozenset[Health] = frozenset({Health.HEALTHY, Health.HYPERTENSION})

#: genotypes of the training data, the wildtype or not reported
CONTROL_GENOTYPES: frozenset[Genotype] = frozenset(
    {
        Genotype.CYP2C9_1_1,
        Genotype.ABCB1_2677GG_3435CC,
        Genotype.NR,
    }
)

#: coadministered drugs without an interaction with losartan: their arms are
#: training data like the arms without a coadministration. Furosemide on the
#: evening before the dose depletes the sodium of Azizi1999 and Doig1993, the
#: baselines of the RAAS and the blood pressure are those of each arm
NO_INTERACTION: frozenset[Coadministration] = frozenset(
    {Coadministration.NONE, Coadministration.FUROSEMIDE}
)

#: observables of the pharmacokinetics
PK_OBSERVABLES: set[str] = {
    "Cveplasma_los",
    "Cveplasma_e3174",
    "Cveplasma_l158",
    "Cveplasma_los_e3174",
    "Cveplasma_total",
    "mr_e3174_los_plasma",
    "Aurine_los",
    "Aurine_e3174",
    "Aurine_l158",
    "Aurine_los_e3174",
    "Aurine_total",
    "mr_e3174_los_urine",
    "mr_los_e3174_urine",
    "Afeces_los",
    "Afeces_e3174",
    "Afeces_l158",
    "Afeces_total",
    "mr_e3174_los_feces",
}

#: observables of the pharmacodynamics, the RAAS and the blood pressure
PD_OBSERVABLES: set[str] = {
    "ren",
    "ren_change",
    "ren_ratio",
    "ang1",
    "ang1_change",
    "ang1_ratio",
    "ang2",
    "ang2_change",
    "ang2_ratio",
    "ald",
    "ald_change",
    "ald_ratio",
    "SBP",
    "SBP_change",
    "SBP_ratio",
    "DBP",
    "DBP_change",
    "DBP_ratio",
    "MAP",
}

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {
    # the change of the aldosterone at 40 mg has no error of its own (Fig 6
    # has the SEM of the placebo and of 120 mg only), the placebo corrected
    # change rises to +17 pg/ml above the baseline at 24-33 hr, where 80 and
    # 120 mg are at -3 to +3; it was 115 of the cost of 242 of the first PD
    # fit (PD_20261019_133445__26dc7), weighted with the SD of the placebo
    "Munafo1992": {"fm_po_aldosterone_corrected_LOS40"},
}

#: the validation of the interactions: the arms with a coadministered drug,
#: simulated with CYP2C9 factors of the literature which are not fitted (see
#: `cyp2c9_interactions` of the experiments), and the arms of phenotyping
#: cocktails. `filter_training` lets the validation data through for this
#: reason, every other interaction stays out of the selection altogether
INTERACTION_VALIDATION: NamedMappings = {
    # phenytoin, induction and inhibition of CYP2C9, 10 daily doses
    "Fischer2002": {
        "fm_po_los50_losartan_LOS50_PHE4",
        "fm_po_los50_exp3174_LOS50_PHE4",
        "fm_po_los50_losartan_urine_LOS50_PHE4",
        "fm_po_los50_exp3174_urine_LOS50_PHE4",
    },
    # cimetidine, no factor (CYP3A4: losartan AUC x1.18, E3174 x1.00)
    "Goldberg1995": {"fm_po_LOS100_CIM400_losartan", "fm_po_LOS100_CIM400_exp3174"},
    # silymarin, inhibition of CYP2C9
    "Han2009a": {
        "fm_po_los50_losartan_1_SIL140",
        "fm_po_los50_exp3174_1_SIL140",
        "fm_po_los50_losartan_3_SIL140",
        "fm_po_los50_exp3174_3_SIL140",
    },
    # breviscapine, no effect in the data, simulated without one
    "Huang2021": {
        "fm_po_los50_los_breviscapine_1",
        "fm_po_los50_e3174_breviscapine_1",
        "fm_po_los50_los_breviscapine_3",
        "fm_po_los50_e3174_breviscapine_3",
    },
    # a phenotyping cocktail (control) and red ginseng (CYP2C9)
    "Kim2016": {
        "fm_po_los50_los_control",
        "fm_po_los50_e3174_control",
        "fm_po_los50_los_ginseng",
        "fm_po_los50_e3174_ginseng",
    },
    # bucolome, inhibition of CYP2C9
    "Kobayashi2008": {"fm_po_los25_los_LOS_BUC", "fm_po_los25_e3174_LOS_BUC"},
    # the healthy controls of a phenotyping cocktail (a compounded capsule)
    "Puris2019": {
        "fm_po_los12p5_losartan_healthy",
        "fm_po_los12p5_exp3174_healthy",
    },
    # the Inje phenotyping cocktail
    "Tanaka2014": {
        "fm_po_los50_losartan",
        "fm_po_los50_exp3174",
        "fm_po_los50_losartan_urine",
        "fm_po_los50_exp3174_urine",
    },
}

#: the validation of the doses: the repeated doses and the doses of other
#: conditions, which the single doses of the training predict
DOSE_VALIDATION: NamedMappings = {
    # 12.5 mg in the Basel cocktail period with omeprazole and midazolam (no
    # interaction, AUC ratio of the full cocktail 0.99, 0.90-1.10): the
    # absorption is slower than in every other study (losartan Cmax 2.5-fold
    # below the model at 1 hr, E3174 tmax 6-8 hr instead of 3-4 hr). As
    # training data it was 23 % of the cost of the PK fit and lengthened the
    # half-life of losartan to 3.2 hr (Vss 73 l, Lo1995: 2.1 hr, 34 l)
    "Donzelli2014": {
        "fm_po_los12p5_los_los_ome_mid",
        "fm_po_los12p5_e3174_los_ome_mid",
    },
    # the single doses of salt replete subjects (PD) and eight daily doses
    "Christen1991a": {
        *(
            f"fm_po_{sid}_LOS{dose}"
            for sid in ["ang2", "aldosterone"]
            for dose in ["10", "20", "40"]
        ),
        *(
            f"fm_po_{sid}_LOS{dose}M"
            for sid in ["pra", "ang2"]
            for dose in ["5", "10", "20", "40"]
        ),
    },
    # seven daily doses of 100 mg, days 2-6 after breakfast
    "Ohtawa1993": {
        "fm_po_losartan_LOS100M_multi",
        "fm_po_exp3174_LOS100M_multi",
        "fm_po_losartan_urine_LOS100M_multi",
        "fm_po_exp3174_urine_LOS100M_multi",
        "fm_po_pra_ratio_LOS100M_multi",
        "fm_po_ang2_LOS100M_multi",
        "fm_po_sbp_LOS100M_multi",
        "fm_po_dbp_LOS100M_multi",
    },
    # 100 mg once daily for 7 days, hypertensive patients with a normal renal
    # function on their diuretics
    "Sica1995": {
        "fm_po_los100_losartan_group1",
        "fm_po_los100_exp3174_group1",
        "fm_po_los100_losartan_urine_group1",
        "fm_po_los100_exp3174_urine_group1",
    },
}

#: the validation of the populations: the genotypes, obesity, cirrhosis and
#: renal impairment, simulated with factors which are not fitted (the CYP2C9
#: activity of the genotypes, the body weight, `cirrhosis_map`, `renal_map`)
POPULATION_VALIDATION: NamedMappings = {
    # sodium depleted with furosemide, a light breakfast: the exposure is lower
    # than the one of the other 50 mg studies
    "Azizi1999": {"fm_po_LOS50_losartan", "fm_po_LOS50_exp3174"},
    # CYP2C9 *3/*3, one subject
    "Bae2011": {"fm_po_los50_los_3", "fm_po_los50_e3174_3"},
    # cirrhosis, Child-Pugh A and B
    "FDA1995S67": {
        f"fm_{sid}_{arm}"
        for arm in ["LOS50PO", "LOS10IV"]
        for sid in ["losartan", "losartan_urine", "exp3174", "exp3174_urine"]
    }
    | {"fm_exp3174_E10IV", "fm_exp3174_urine_E10IV"},
    # CYP2C9 *1/*3
    "Han2009a": {"fm_po_los50_losartan_3", "fm_po_los50_exp3174_3"},
    "Huang2021": {"fm_po_los50_los_placebo_3", "fm_po_los50_e3174_placebo_3"},
    # CYP2C9 *1/*2 and *1/*3
    "Lee2003b": {
        f"fm_po_los50_{sid}_{g}" for sid in ["los", "e3174"] for g in ["2", "3"]
    },
    # CYP2C9 *1/*3 and *1/*13
    "Li2009": {
        f"fm_po_los50_{sid}_{g}" for sid in ["los", "e3174"] for g in ["1_3", "1_13"]
    },
    # obesity before and one year after a gastric bypass, phenotyping cocktail
    "Puris2019": {
        f"fm_po_los12p5_{sid}_{arm}"
        for sid in ["losartan", "exp3174", "losartan_urine", "exp3174_urine"]
        for arm in ["obese", "one_year_after"]
    },
    # CYP2C9 *1/*3
    "Sekino2003": {
        "fm_po_los25_exp3174_losartan_plasma_3",
        "fm_po_los25_exp3174_losartan_urine_3",
    },
    # ABCB1 2677/3435 GT/CT and TT/TT
    "Shin2020": {
        f"fm_po_los50_{sid}_{g}"
        for sid in ["losartan", "exp3174", "losartan_exp3174", "losartan_exp3174_urine"]
        for g in ["GT_CT", "TT_TT"]
    },
    # renal impairment (CLcr 59 and 17 ml/min) on 7 daily doses and diuretics
    "Sica1995": {
        f"fm_po_los100_{sid}_{group}"
        for sid in ["losartan", "exp3174", "losartan_urine", "exp3174_urine"]
        for group in ["group2", "group3"]
    },
    # CYP2C9 genotypes of the Swedish (50 mg) and the Spanish (25 mg) subjects
    "Yasar2002a": {
        *(
            f"fm_po_los50_{sid}_{g}"
            for sid in ["losartan", "e3174"]
            for g in ["2_2", "3_3"]
        ),
        *(
            f"fm_po_los50_{sid}_{g}"
            for sid in ["losartan_urine", "e3174_urine", "losartan_exp3174_urine"]
            for g in ["1_2", "1_3", "2_2", "2_3", "3_3"]
        ),
        *(
            f"fm_po_los25_losartan_exp3174_urine_{g}_spanish"
            for g in ["1_2", "1_3", "3_3"]
        ),
    },
}

#: the validation of the limitations of the model: valid data of a process the
#: model does not describe, which would distort the fitted parameters of the
#: processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    # the placebo arms: meals, posture and the circadian rhythm change the RAAS
    # and the blood pressure, the model without a dose is flat
    "Azizi1999": {
        "fm_po_placebo_renin",
        "fm_po_placebo_ang1",
        "fm_po_placebo_ang2",
        "fm_po_placebo_map",
    },
    "Christen1991a": {
        "fm_po_ang2_placebo",
        "fm_po_aldosterone_placebo",
        "fm_po_pra_placeboM",
        "fm_po_ang2_placeboM",
    },
    "Doig1993": {
        "fm_po_placebo_pra",
        "fm_po_placebo_aldosterone",
        "fm_po_placebo_supine_sbp",
        "fm_po_placebo_supine_dbp",
    },
    "Goldberg1995": {"fm_po_control_pra_ratio"},
    # the placebo and the diastolic pressure of the weeks of daily doses: the
    # data are flat from trough to peak, the direct effect of the model returns
    # to the baseline before each dose
    "Goldberg1995a": {
        "fm_po_placebo_ratio_pra",
        "fm_po_placebo_ang2_ratio",
        "fm_po_placebo_ald_change",
        "fm_po_placebo_dbp_change",
        "fm_po_LOS25_dbp_change",
        "fm_po_LOS100_dbp_change",
    },
    # the L-158,641 of the radioactivity is a part of the other metabolites of
    # the model; the fecal losartan and E3174 after the infusion are biliary
    # losartan, which the model does not have (two subjects)
    "FDA1995S60": {
        "fm_l158_urine_LOS100PO",
        "fm_l158_urine_LOS30IV",
        "fm_l158_feces_LOS30IV",
        "fm_losartan_feces_LOS30IV",
        "fm_exp3174_feces_LOS30IV",
    },
    "Munafo1992": {"fm_po_aldosterone_Placebo"},
    # the change of the blood pressure from the predose value without a
    # placebo: SBP -9 to -15 and DBP -11 to -19 mmHg at 6-12 hr after 25 mg in
    # healthy normotensive subjects, more than the trough to peak fall of
    # 100 mg in hypertension (DBP -6 to -10 mmHg [FDA label]); the circadian
    # fall of the pressure is in the data. As training data it pushed BP_emax
    # to its upper bound (0.6, PD_20261019_133445__26dc7)
    "Sekino2003": {
        f"fm_po_los25_{bp}_change_{g}" for bp in ["sbp", "dbp"] for g in ["1", "3"]
    },
    # a microdose in a cocktail capsule, absorbed much slower than a tablet
    "Oh2012": {"fm_po_los2_losartan", "fm_po_los2_exp3174"},
    # the aldosterone falls by as much in the placebo arm (circadian rhythm),
    # the single dose study has no placebo
    "Ohtawa1993": {
        *(f"fm_po_aldosterone_LOS{dose}" for dose in ["25", "50", "100", "200"]),
        "fm_po_aldosterone_LOS100M_multi",
        "fm_po_pra_ratio_placebo_multi",
        "fm_po_ang2_placebo_multi",
        "fm_po_aldosterone_placebo_multi",
        "fm_po_sbp_placebo_multi",
        "fm_po_dbp_placebo_multi",
    },
}

#: the summaries of data which the fits use one by one, excluded from every
#: fit: a mean of subjects would count the subjects twice
AGGREGATE_MAPPINGS: NamedMappings = {
    # arm E, the same 16 subjects as arm A
    "Donzelli2014": {
        "fm_po_los12p5_los_basel_cocktail",
        "fm_po_los12p5_e3174_basel_cocktail",
    },
    # the same six subjects as study 2 of Lo1995, which has the plasma of the
    # infusions and the urine of losartan and E3174
    "FDA1995S60": {
        "fm_losartan_LOS30IV",
        "fm_exp3174_LOS30IV",
        "fm_exp3174_E20IV",
        "fm_losartan_urine_LOS100PO",
        "fm_exp3174_urine_LOS100PO",
        "fm_losartan_urine_LOS30IV",
        "fm_exp3174_urine_LOS30IV",
    },
    # the pooled group, the same 11 subjects as the genotype groups
    "Huang2021": {
        "fm_po_los50_los_placebo",
        "fm_po_los50_e3174_placebo",
        "fm_po_los50_los_breviscapine",
        "fm_po_los50_e3174_breviscapine",
    },
    # the ratios of the two mean curves, derived by the curation
    "Oh2012": {"fm_po_los2_exp3174_losartan"},
    "Tanaka2014": {"fm_po_los50_exp3174_losartan"},
}

#: the validation data by the question it answers
VALIDATION_QUESTIONS: dict[str, NamedMappings] = {
    "interactions": INTERACTION_VALIDATION,
    "doses": DOSE_VALIDATION,
    "populations": POPULATION_VALIDATION,
    "limitations": LIMITATION_VALIDATION,
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(mappings.get(study, set()) for mappings in VALIDATION_QUESTIONS.values())
    )
    for mappings in VALIDATION_QUESTIONS.values()
    for study in mappings
}


def study(fit_mapping: FitMapping) -> str:
    """Get the study of a fit mapping, the id of its experiment."""
    return fit_mapping.reference.experiment.sid


def named(
    mappings: NamedMappings, fit_mapping_key: str, fit_mapping: FitMapping
) -> bool:
    """Check whether a fit mapping is one of the named mappings of its study."""
    return fit_mapping_key in mappings.get(study(fit_mapping), set())


def check_named(mappings: NamedMappings, fit_mappings: FitMappings) -> None:
    """Check that the named mappings are fit mappings of their studies.

    Raises:
        ValueError: for a study which is no experiment, or a key which is no fit
            mapping of its study.
    """
    for sid, keys in mappings.items():
        experiment = fit_mappings.experiments.get(sid)
        if experiment is None:
            raise ValueError(f"'{sid}' is no experiment of the fits.")
        if unknown := sorted(keys - set(experiment._fit_mappings)):
            raise ValueError(f"'{unknown}' are no fit mappings of '{sid}'.")


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No summaries of data which is fitted one by one, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


def filter_training(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only the training data and the validation data.

    The training data are the healthy or hypertensive subjects with the CYP2C9
    wildtype or a genotype which was not reported, without an interacting drug;
    the validation data (`VALIDATION_MAPPINGS`) passes, every other mapping
    stays out of the selection.
    """
    if named(VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping):
        return True
    metadata = mapping_metadata(fit_mapping)
    return (
        metadata.health in CONTROL_HEALTH
        and metadata.genotype in CONTROL_GENOTYPES
        and metadata.coadministration in NO_INTERACTION
    )


#: Only the pharmacokinetic data.
filter_pk = observable_filter(PK_OBSERVABLES)


#: Only the pharmacodynamic data.
filter_pd = observable_filter(PD_OBSERVABLES)


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the data of a fit, the mappings which fail them are
    excluded, and so are the summaries of data which is fitted one by one
    (`AGGREGATE_MAPPINGS`) and the mappings which are neither training nor
    validation data (`filter_training`). The outliers are `OUTLIER_MAPPINGS`,
    the validation data is `VALIDATION_MAPPINGS`. Neither kind enters the cost;
    both are simulated and evaluated with the training data, so a report says
    how the fitted parameters describe data they were not fitted on.

    Args:
        filters: filters of the data of the fit.
        validation: the keys of the validation data by study.

    Returns:
        The fit mapping collections of all kinds by experiment id.

    Raises:
        ValueError: if an outlier is a mapping of no study or of several.
    """
    fit_mappings = FitMappings(
        experiment_classes=EXPERIMENT_CLASSES,
        base_path=LOSARTAN_PATH,
        data_path=DATA_PATHS,
    )
    for mappings in [
        OUTLIER_MAPPINGS,
        AGGREGATE_MAPPINGS,
        validation,
    ]:
        check_named(mappings, fit_mappings)
    # sbmlsim selects the outliers by their keys alone, which must name one
    # mapping of all studies
    outliers = {key for keys in OUTLIER_MAPPINGS.values() for key in keys}
    for key in outliers:
        studies = [
            sid
            for sid, experiment in fit_mappings.experiments.items()
            if key in experiment._fit_mappings
        ]
        if len(studies) != 1:
            raise ValueError(f"The outlier '{key}' is a mapping of {studies}.")
    selection = fit_mappings.select(
        filters=[
            *(filters if isinstance(filters, list) else [filters]),
            filter_aggregates,
            filter_training,
        ],
        outliers=outliers,
        validation=lambda key, fit_mapping: named(validation, key, fit_mapping),
        print_info=False,
    )
    # the counts per experiment and kind, not the single mappings
    selection.print(detail=False)
    # the validation mappings by the question they answer
    selected = {
        (sid, key)
        for sid, kinds in selection.kinds.items()
        for key, kind in kinds.items()
        if kind is MappingKind.VALIDATION
    }

    def pairs(mappings: NamedMappings) -> set[tuple[str, str]]:
        return {(sid, key) for sid, keys in mappings.items() for key in keys}

    console.print(
        "validation: "
        + ", ".join(
            f"{len(selected & pairs(mappings))} {question}"
            for question, mappings in VALIDATION_QUESTIONS.items()
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """Pharmacokinetic data."""
    return mapping_collections([filter_pk])


def f_collections_pd() -> dict[str, list[FitMappingCollection]]:
    """Pharmacodynamic data."""
    return mapping_collections([filter_pd])


if __name__ == "__main__":
    for f in [f_collections_pk, f_collections_pd]:
        console.rule(title=f.__name__, align="left", style="white")
        f()
