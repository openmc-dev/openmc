from __future__ import annotations

import openmc
from openmc.data import half_life


def _waste_classification(mat: openmc.Material, metal: bool = True) -> str:
    """Classify a material for near-surface waste disposal.

    This method determines a waste classification for a material based on the
    NRC regulations (10 CFR 61.55).

    Parameters
    ----------
    mat : openmc.Material
        The material to classify.
    metal : bool, optional
        Whether or not the material is in metal form. This changes the
        acceptable limits in Tables 1 and 2 for certain nuclides.

    Returns
    -------
    str
        The waste disposal classification, which can be "Class A", "Class B",
        "Class C", or "GTCC" (greater than class C).

    """
    # Determine metrics based on Tables 1 and 2 using sum of fractions rule for
    # mixture of radionuclides from §61.55(a)(7)
    ratio1 = _waste_disposal_rating(mat, 'NRC_long', metal=metal)
    ratio2 = [
        _waste_disposal_rating(mat, 'NRC_short_A', metal=metal),
        _waste_disposal_rating(mat, 'NRC_short_B', metal=metal),
        _waste_disposal_rating(mat, 'NRC_short_C', metal=metal),
    ]

    # Determine which nuclides are present in Table 1 and Table 2
    table1_nuclides_present = (ratio1 > 0.0)
    table2_nuclides_present = any(x > 0.0 for x in ratio2)

    # Helper function for classifying based on Table 2
    def classify_table2(col1, col2, col3):
        if col1 < 1.0:
            return "Class A"
        elif col2 < 1.0:
            return "Class B"
        elif col3 < 1.0:
            return "Class C"
        else:
            return "GTCC"

    if table1_nuclides_present and table2_nuclides_present:
        # Classification based on §61.55(a)(5)
        if ratio1 < 0.1:
            return classify_table2(*ratio2)
        elif ratio1 < 1.0:
            return "Class C" if ratio2[2] < 1.0 else "GTCC"
        else:
            return "GTCC"

    elif table1_nuclides_present:
        # Classification based on §61.55(a)(3)
        if ratio1 < 0.1:
            return "Class A"
        elif ratio1 < 1.0:
            return "Class C"
        else:
            return "GTCC"

    elif table2_nuclides_present:
        # Classification based on §61.55(a)(4)
        return classify_table2(*ratio2)

    else:
        # Classification based on §61.55(a)(6)
        return "Class A"


def _waste_disposal_rating(
        mat: openmc.Material,
        limits: str | dict[str, float] = 'Fetter',
        metal: bool = False,
    ) -> float:
    """Return the waste disposal rating for a material.

    This method returns a waste disposal rating for the material based on a set
    of specific activity limits. The waste disposal rating is a single number
    that represents the sum of the ratios of the specific activity for each
    radionuclide in the material against a nuclide-specific limit. A value less
    than 1.0 indicates that the material "meets" the limits whereas a value
    greater than 1.0 exceeds the limits.

    Parameters
    ----------
    mat : openmc.Material
        The material to classify.
    limits : str or dict, optional
        The name of a predefined set of specific activity limits or a dictionary
        that contains specific activity limits for radionuclides, where keys are
        nuclide names and values are activities in units of [Ci/m3]. The
        predefined options are:

        - 'Fetter': Uses limits from Fetter et al. (1990)
        - 'NRC_long': Uses the 10 CFR 61.55 limits for long-lived radionuclides
        - 'NRC_short_A': Uses the 10 CFR 61.55 class A limits for short-lived
          radionuclides
        - 'NRC_short_B': Uses the 10 CFR 61.55 class B limits for short-lived
          radionuclides
        - 'NRC_short_C': Uses the 10 CFR 61.55 class C limits for short-lived
          radionuclides
    metal : bool, optional
        Whether or not the material is in metal form (only applicable for NRC
        based limits)

    Returns
    -------
    float
        The waste disposal rating for the material.

    """
    if limits == 'Fetter':
        # Specific activity limits for radionuclides with half-lives between 5
        # years and 1e12 years from Table 2 in Fetter
        limits = {
            "Be10": 5.0e3,
            "C14": 6.0e2,
            "Al26": 9.0e-2,
            "Si32": 6.0e2,
            "Cl36": 1.0e1,
            "Ar39": 2.0e4,
            "Ar42": 2.0e4,
            "K40": 2.0e0,
            "Ca41": 1.0e4,
            "Ti44": 2.0e2,
            "Fe60": 1.0e-1,
            "Co60": 3.0e8,
            "Ni59": 9.0e2,
            "Ni63": 7.0e5,
            "Se79": 5.0e1,
            "Kr81": 3.0e1,
            "Sr90": 8.0e5,
            "Nb91": 2.0e2,
            "Nb92": 2.0e-1,
            "Nb94": 2.0e-1,
            "Mo93": 4.0e3,
            "Tc97": 4.0e-1,
            "Tc98": 1.0e-2,
            "Tc99": 6.0e-2,
            "Pb107": 9.0e2,
            "Ag108_m1": 3.0e0,
            "Sn121_m1": 7.0e5,
            "Sn126": 1.0e-1,
            "I129": 2.0e0,
            "Cs137": 5.0e4,
            "Ba133": 2.0e8,
            "La137": 2.0e2,
            "Sm151": 5.0e7,
            "Eu150_m1": 3.0e3,
            "Eu152": 3.0e5,
            "Eu154": 5.0e6,
            "Gd148": 2.0e5,
            "Gd150": 2.0e3,
            "Tb157": 5.0e3,
            "Tb158": 4.0e0,
            "Dy154": 1.0e3,
            "Ho166_m1": 2.0e-1,
            "Hf178_m1": 9.0e3,
            "Hf182": 2.0e-1,
            "Re186_m1": 2.0e1,
            "Ir192_m1": 1.0e0,
            "Pt193": 2.0e8,
            "Hg194": 5.0e-1,
            "Pb202": 6.0e-1,
            "Pb210": 3.0e7,
            "Bi207": 9.0e3,
            "Bi208": 8.0e-2,
            "Bi210_m1": 1.0e0,
            "Po209": 3.0e3,
            "Ra226": 1.0e-1,
            "Ra228": 3.0e7,
            "Ac227": 5.0e5,
            "Th229": 2.0e0,
            "Th230": 3.0e-1,
            "Th232": 1.0e-1,
            "Pa231": 7.0e-1,
            "U232": 3.0e1,
            "U233": 2.0e1,
            "U234": 9.0e1,
            "U235": 2.0e0,
            "Np236": 1.0e0,
            "Np237": 1.0e0,
            "Pu238": 7.0e4,
            "Pu239": 1.0e3,
            "Pu240": 1.0e3,
            "Pu241": 2.0e3,
            "Pu242": 1.0e3,
            "Pu244": 9.0e-1,
            "Am241": 5.0e1,
            "Am242_m1": 3.0e2,
            "Am243": 2.0e0,
            "Cm243": 6.0e2,
            "Cm244": 5.0e5,
            "Cm245": 5.0e0,
            "Cm246": 8.0e2,
            "Cm248": 8.0e2,
        }

    elif limits == 'NRC_long':
        # Specific activity limits for long-lived radionuclides from Table 1 in
        # 10 CFR 61.55 in Ci/m3.
        limits = {
            'C14': 8.0,
            'Tc99': 3.0,
            'I129': 0.08,
        }
        if metal:
            limits['C14'] = 80.0
            limits['Ni59'] = 220.0
            limits['Nb94'] = 0.2

        # Convert values in nCi/g to Ci/m3
        factor = (1e6 * mat.get_mass_density()) / 1e9
        limits.update({
            'Pu241': 3500.0 * factor,
            'Cm242': 20000.0 * factor,
            'Np237': 100.0 * factor,
            'Pu238': 100.0 * factor,
            'Pu239': 100.0 * factor,
            'Pu240': 100.0 * factor,
            'Pu242': 100.0 * factor,
            'Pu244': 100.0 * factor,
            'Am241': 100.0 * factor,
            'Am243': 100.0 * factor,
            'Cm243': 100.0 * factor,
            'Cm244': 100.0 * factor,
            'Cm245': 100.0 * factor,
            'Cm246': 100.0 * factor,
            'Cm247': 100.0 * factor,
            'Cm248': 100.0 * factor,
            'Bk247': 100.0 * factor,
            'Cf249': 100.0 * factor,
            'Cf250': 100.0 * factor,
            'Cf251': 100.0 * factor,
        })

    elif limits == 'NRC_short_A':
        # Get Class A specific activity limits for short-lived radionuclides
        # from Table 2 in 10 CFR 61.55
        limits = {
            'H3': 40.0,
            'Co60': 700.0,
            'Ni63': 35.0 if metal else 3.5,
            'Sr90': 0.04,
            'Cs137': 1.0
        }

        # Add radionuclides with half-lives < 5 years to limits for class A
        five_years = 60.0 * 60.0 * 24.0 * 365.25 * 5.0
        for nuc in mat.get_nuclides():
            if half_life(nuc) is not None and half_life(nuc) < five_years:
                limits[nuc] = 700.0

    elif limits == 'NRC_short_B':
        # Get Class B specific activity limits for short-lived radionuclides
        # from Table 2 in 10 CFR 61.55
        limits = {'Ni63': 700.0 if metal else 70.0, 'Sr90': 150.0, 'Cs137': 44.0}

    elif limits == 'NRC_short_C':
        # Get Class C specific activity limits for short-lived radionuclides
        # from Table 2 in 10 CFR 61.55
        limits = {'Ni63': 7000.0 if metal else 700.0, 'Sr90': 7000.0, 'Cs137': 4600.0}

    # Calculate the sum of the fractions of the activity of each radionuclide
    # compared to the specified limits
    ratio = 0.0
    for nuc, ci_m3 in mat.get_activity(units="Ci/m3", by_nuclide=True).items():
        if nuc in limits:
            ratio += ci_m3 / limits[nuc]
    return ratio
