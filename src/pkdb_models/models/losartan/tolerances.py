"""The tolerances of the integrator for the losartan model."""

from sbmlsim.model.tolerances import AbsoluteTolerance

#: the states of the renin-angiotensin-aldosterone system. Their amounts are
#: about 1e-13 to 1e-8 mmole (renin 0.1 pM to 20 pM, angiotensin I and II
#: about 10 pM, aldosterone up to 2 nM in 3 l of plasma), far below the amounts
#: of losartan and its metabolites, so they get absolute tolerances of their
#: own scale
RAAS_STATES: tuple[str, ...] = ("ren", "ang1", "ang2", "ald")


def absolute_tolerance(value: float, species: float, raas: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the losartan model.

    The species of losartan and its metabolites are small: the plasma
    concentrations are about 1e-6 to 1e-3 mM (1e-5 mM after the microdose of
    2 mg of the cocktail study Oh2012), the amounts about 1e-6 to 1e-1 mmole.
    The other states are the doses in mg (the targets of the rate rules).

    Args:
        value: tolerance of the other states [mg].
        species: tolerance of the species, of their concentrations [mM] and
            amounts [mmole].
        raas: tolerance of the amounts of the RAAS states [mmole].

    Returns:
        The absolute tolerance per kind of state and of the RAAS states.
    """
    return AbsoluteTolerance(
        amount=species,
        concentration=species,
        other=value,
        ids=dict.fromkeys(RAAS_STATES, raas),
    )


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(1e-10, species=1e-14, raas=1e-18)
