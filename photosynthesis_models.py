"""
C3 and C4 photosynthesis models plus the coupled stomatal-biochemical model
for diagnosing the differential CO2 fertilization effect (CFE) of soybean
(C3) and maize (C4).

- Unit tests at module bottom (run by `python photosynthesis_models.py`).

═══════════════════════════════════════════════════════════════════════════
1.  WHAT THE MODEL DOES
═══════════════════════════════════════════════════════════════════════════

The model couples a biochemical photosynthesis model (Farquhar–von
Caemmerer–Berry for C3, von Caemmerer 2000 for C4) to a stomatal sub-model.

C3: Farquhar, von Caemmerer, & Berry (1980), with parameter values
    from Bernacchi et al. (2001, 2002) for soybean.
C4: von Caemmerer (2000), with parameter values from Massad et al. (2007),
    Leakey et al. (2006), and von Caemmerer (2013) for maize.


Three independent physiological controls modulate the response of leaf-
level A to atmospheric CO₂ (Ca), atmospheric VPD, and root-zone soil water
content (SWC):

    (I)   STOMATAL APERTURE (gs)        — anchored Medlyn USO × β(SWC)
    (II)  C3 BIOCHEMICAL CAPACITY       — β-ramped Vcmax & Jmax (Flexas)
    (III) C4 CCM CAPACITY (Vpmax)       — β-scaled PEPC under drought

The differential CFE response between species, and its collapse under
drought, emerge from the interaction of these three mechanisms with the
species-specific A-Ci curve geometry. No species-asymmetric stomatal rule
is needed.   

=============================================================================
2.  STOMATAL CONDUCTANCE (gs)
=============================================================================

2.1  Functional form -- anchored Medlyn USO

For both species:

    gs(Ca, VPD, SWC) = beta(SWC) x gs_fc_obs x kernel_ca(Ca, VPD, g1)

where:
    kernel_ca(Ca) = gs_medlyn(Ca, VPD, g1) / gs_medlyn(Ca_ref, VPD, g1)
    gs_medlyn     = 1.6 x (1 + g1/sqrt(VPD)) x A(Ci=Ca*g1/(g1+sqrt(VPD))) / Ca
    Ca_ref        = 370 umol mol^-1

The kernel captures the dimensionless Ca-dependence of stomatal aperture
predicted by Medlyn USO theory. The absolute magnitude is set empirically
by gs_fc_obs (see 2.3). beta(SWC) is a Verhoef-Egea soil-water stress
factor (see section 3).

2.2  Why anchor magnitude empirically rather than predict it

The Medlyn USO derivation is theoretically robust on the SHAPE of the
gs(Ca, VPD) response -- Ci/Ca ~= g1/(g1+sqrt(VPD)) is preserved across Ca
changes -- but pinning down absolute gs requires the full equation
gs = 1.6*(1+g1/sqrt(VPD))*A/Ca, which inherits absolute errors in A. The
no-mesophyll-conductance photosynthesis model used here under-predicts
absolute A, which in turn under-predicts absolute gs by ~15-30%. By
anchoring gs_fc to direct field measurement and using the kernel only as
a Ca-scaling function, we sidestep this magnitude problem while preserving
the theoretical Ca-dependence.

2.3  Empirical anchors (gs at well-watered, Ca_ref=370 ppm, VPD=1.5 kPa)

    gs_fc_soy   = 0.55 mol H2O m^-2 s^-1    Ainsworth & Long (2005)
    gs_fc_maize = 0.30 mol H2O m^-2 s^-1    Leakey et al. (2006)

Maize anchor: peak-season midday gs from Leakey et al. (2006, Plant
Physiology) Table II, averaging the three reproductive-stage measurements
(DOY 190, 201, 215: 0.27, 0.29, 0.36; mean 0.31), excluding cloud-
disrupted DOY 173 and senescence-phase DOY 229.

Soy anchor: Ainsworth & Long (2005) FACE meta-analytic mean. The model
operating point at this gs gives Ci/Ca = 0.77 versus directly observed
~0.73 -- a small offset attributable to omitted mesophyll conductance.
Anchoring to the directly measured gs is preferred to anchoring to the
derived Ci/Ca.

2.4  Selection of the Medlyn slope g1

The Medlyn slope g1 controls the dimensionless Ca-dependence kernel.
In the anchored Medlyn formulation, g1 does not affect the operating
point at Ca_ref (which is pinned by gs_fc); it only affects how strongly
the kernel scales as Ca departs from Ca_ref.

We use the plant-functional-type (PFT) values from Lin et al. (2015,
Nat. Clim. Change 5:459-464, DOI 10.1038/nclimate2550), a synthesis of
314 species across 56 field studies.

    g1_soy   = 5.79 kPa^0.5   Lin C3 cropland (+/- 0.64)
    g1_maize = 1.62 kPa^0.5   Lin C4 grassland (+/- 0.13)

Lin's database does not contain a separate "C4 cropland" PFT; we use
C4 grassland as the closest match for maize.

2.5  Why the same Medlyn rule produces different Ca-responses

The kernel_ca behaves differently for C3 vs C4 because of the
species-specific A-Ci curve geometry. The two species sit in
different regions of their respective A-Ci curves and therefore
respond to Ca-changes differently.

  C3 SOYBEAN (g1=5.79, operating Ci/Ca ~= 0.77):
    Soy operates in the rising portion of its A-Ci curve where A is
    near-linear in Ci. As Ca rises, both Ci and A rise nearly
    proportionally, so gs_medlyn = const*A/Ca is only weakly Ca-dependent
    in the historical range. Across 336->420 ppm at VPD=1.5 kPa, the
    kernel changes by only -4.5%. Across the wider SoyFACE comparison
    range 376->550 ppm, the curvature of the A-Ci curve becomes
    perceptible and the kernel decreases by -20%, falling within the
    range of directly measured midday gs reductions at SoyFACE
    (Bernacchi et al. 2007: -16%; Long et al. 2006: range -10% to -30%
    across measurement protocols; Ainsworth & Long 2005 meta-analytic
    mean for C3 species: -22%).

  C4 MAIZE (g1=1.62, operating Ci/Ca ~= 0.37):
    Maize operates on the saturated plateau of its A-Ci curve where A
    is nearly Ca-independent. As Ca rises, A stays flat, so
    gs_medlyn = const*A/Ca falls inversely with Ca over the entire
    Ca range. The kernel decreases by -20% across 336->420 ppm and
    -32% across 376->550 ppm. The latter agrees with Leakey et al.
    (2006b) Table II directly measured midday gs reduction of -34%.

The key C3 vs C4 contrast: maize gs falls steeply over the full Ca range
(because A is saturated everywhere), while soy gs falls gently over the
historical range and more steeply when Ca approaches the upper Ci where
A-Ci curvature increases.

2.6  Validation -- model predictions vs. published FACE observations

With Lin et al. (2015) PFT-default g1 values, the model is compared to
SoyFACE and maize-FACE midday observations at ambient (376 ppm) ->
elevated (550 ppm), well-watered, VPD=1.5 kPa. The comparison is for
qualitative agreement on direction and order of magnitude -- not parameter
calibration.

  Soybean (vs SoyFACE):
    Dgs:    model -20%   |  obs -16% (Bernacchi 2007), -22% (A&L 2005 meta)
    DCi/Ca: preserved    |  obs preserved
    DA:     model +22%   |  obs +15-25% (SoyFACE meta)

  Maize (vs Leakey et al. 2006b):
    Dgs:    model -32%   |  obs -34% (Tab II midday)
    DCi/Ca: preserved    |  obs preserved
    DA:     model +1%    |  obs ~0% (not significant)

All predicted responses fall within the observed ranges across all
measured quantities for both species. The model is fit-for-purpose for
the diagnostic question: explaining the differential CFE response and
its drought-induced collapse.

=============================================================================
3.  SOIL-WATER STRESS -- Verhoef-Egea linear beta
=============================================================================

    beta(SWC) = clip((SWC - SWC_wp) / (SWC_fc - SWC_wp), 0, 1)

Standard linear beta factor from Verhoef & Egea (2014), Egea et al. (2011);
the form used in CLM, JULES, ORCHIDEE.

Anchors for Corn Belt silt loam (Mollisol):
    SWC_wp = 0.10 m3 m^-3   permanent wilting point
    SWC_fc = 0.31 m3 m^-3   field capacity

Values are from pedotransfer functions for silt loam texture (Saxton &
Rawls 2006; Allen et al. 1998, FAO-56), consistent with field observations
of well-watered VWC at SoyFACE/maize-FACE sites in central Illinois.
Reference values for silt loam at -33 kPa span ~0.28-0.32 m3 m^-3 across
sources; 0.31 reflects empirical data from the SoyFACE site.

beta multiplies BOTH gs (stomatal closure under stress) AND biochemical
capacity (see section 4). The same beta factor for both pathways is a
simplifying assumption that captures the leading-order coupling between
leaf hydraulic status and metabolic capacity without introducing a second
free parameter.

=============================================================================
4.  BIOCHEMICAL DROUGHT RESPONSE -- beta-ramped capacity
=============================================================================

4.1  C3 soybean -- Vcmax and Jmax

    Vcmax_eff = Vcmax x [f_min + (1 - f_min) x beta]
    Jmax_eff  = Jmax  x [f_min + (1 - f_min) x beta]
    f_min = 0.60

The functional form follows the conceptual framework of Flexas & Medrano
(2002, Ann. Bot. 89:183-189), who established that drought response in C3
plants progresses through phases: stomatal limitation dominates at mild-
to-moderate stress (with photosynthetic capacity largely intact), while
metabolic limitations (declines in RuBP regeneration capacity, Rubisco
activity) become dominant only at severe stress. Vcmax and Jmax do not
collapse to zero even under severe drought because Rubisco protein is not
rapidly degraded -- biochemical capacity declines but a substantial residual
fraction is retained until rewatering.

The retained-capacity floor f_min = 0.60 is consistent with the observation
that mild/moderate drought leaves photosynthetic capacity nearly intact
(small Vcmax/Jmax reductions of ~10-20%; Flexas & Medrano 2002 Fig. 1, 2)
and that even at severe drought, Vcmax/Jmax decline to roughly 40-60% of
well-watered values rather than collapsing entirely (Zhou et al. 2014 Tree
Physiol.). f_min = 0.60 represents a deliberate conservative midpoint --
the model purpose is qualitative diagnosis of the DCFE collapse mechanism,
not precise reproduction of Vcmax/Jmax values, and the central conclusion
(maize CFE rises under drought as Ci moves below the C4 saturation knee)
is robust to f_min values in the 0.40-0.70 range.

4.2  C4 maize -- Vpmax (CCM capacity)

    Vpmax_eff = Vpmax x beta

PEPC activity declines under water stress in C4 grasses (Ripley et al.
2007 J. Exp. Bot.; Ghannoum 2009 Ann. Bot. -- see review Fig. 3 and
references therein for the three-phase stomatal/mixed/non-stomatal
progression; Bellasio 2024). The progressive reduction in PEPC
contributes to the non-stomatal component of C4 drought response.
Direct beta scaling (no f_min) so that beta=0 fully shuts down the CCM.
Vcmax (Rubisco) and Jmax for maize are NOT scaled because Rubisco is
in bundle-sheath cells, anatomically protected from cytosolic
dehydration.

=============================================================================
5.  MONTE-CARLO PARAMETER UNCERTAINTY
=============================================================================

Cultivar variation is represented by sampling biochemical parameters from
biologically realistic ranges:

  Soybean:
    Vcmax ~ U[105, 145]                      umol m^-2 s^-1
    J/V   ~ U[1.45, 1.75]                    centred on 1.60
    Rd    ~ U[1.0, 2.0]                      umol m^-2 s^-1
    Jmax  = (J/V) x Vcmax

  Maize:
    Vcmax ~ U[40, 55],   Vpmax ~ U[95, 135]
    Kp    ~ U[60, 100],  Jmax  ~ U[350, 450]

The Vcmax-Jmax coupling for soy via J/V (rather than independent uniforms)
is essential. Vcmax (capacity) and J/V (balance between carboxylation and
electron transport) are conceptually independent biological factors that
both vary across cultivars; sampling them independently with realistic
ranges captures genuine cultivar variability. Sampling Vcmax and Jmax
independently as uniform variates produces unphysical parameter
combinations (e.g. J/V = 1.2, which does not occur in healthy soybean)
that distort the Monte-Carlo distribution.

The chosen J/V range covers documented soybean variation:
    Bernacchi et al. (2005)  J/V ~= 1.66  (SoyFACE in vivo, no-gm fit)
    Rogers et al. (2004)     J/V = 1.6-1.8  (seasonal variation, SoyFACE)
    Bishop et al. (2015)     J/V ~= 1.5-1.7  (across 18 cultivars)

=============================================================================
6.  REFERENCES
=============================================================================

Ainsworth & Long (2005)        New Phytol.        165: 351-372.
Ainsworth & Rogers (2007)      Plant Cell Environ. 30: 258-270.
Ainsworth, Sanz-Saez, Leisner (2025)  Phil. Trans. R. Soc. B  380: 20240230.
Allen et al. (1998)            FAO Irrigation Drainage Paper 56 (FAO, Rome).
Bellasio (2024)                New Phytol.            (CCM drought response).
Bernacchi et al. (2001)        Plant Cell Environ.  24: 253-259.
Bernacchi et al. (2005)        Planta              220: 434-446.
Bernacchi et al. (2007)        Plant Physiol.      143: 134-144.
Bishop et al. (2015)           Plant Cell Environ.  38: 1765-1774.
De Kauwe et al. (2015)         Geosci. Model Dev.    8: 431-452.
Egea et al. (2011)             Agric. For. Meteorol. 151: 1369-1384.
Farquhar, von Caemmerer, Berry (1980)  Planta     149: 78-90.
Flexas & Medrano (2002)        Ann. Bot.            89: 183-189.
Ghannoum (2009)                Ann. Bot.           103: 635-644.
Leakey et al. (2006a)          Plant Cell Environ.  29: 1794-1800.
Leakey et al. (2006b)          Plant Physiol.      140: 779-790.
Leakey et al. (2009)           J. Exp. Bot.         60: 2859-2876.
Lin et al. (2015)              Nat. Clim. Change     5: 459-464.
Long et al. (2006)             Science             312: 1918-1921.
Medlyn et al. (2011)           Glob. Change Biol.   17: 2134-2144.
Ripley et al. (2007)           J. Exp. Bot.         58: 1351-1363.
Rogers et al. (2004)           Plant Cell Environ.  27: 449-458.
Saxton & Rawls (2006)          Soil Sci. Soc. Am. J. 70: 1569-1578.
Verhoef & Egea (2014)          Agric. For. Meteorol. 191: 22-32.
von Caemmerer (2000)           Biochemical Models of Leaf Photosynthesis.
Zhou et al. (2014)             Tree Physiol.        34: 1035-1046.
"""

import numpy as np
from scipy.optimize import brentq


# ===========================================================================
# Constants -- CO2 levels, Ca anchor
# ===========================================================================
CA_LEVELS = [336, 370, 420]
CA_LABELS = ['1979', '2000', '2023']
CA_COLORS = ['#bbbbbb', '#666666', '#1a1a1a']
CA_LOW    = 336
CA_HIGH   = 420
CA_REF    = 370   # reference Ca for gs anchor
VPD_REF   = 1.5   # Anchored to FACE obs

# ===========================================================================
# Stomatal anchors -- well-watered midday gs from FACE observations
# ===========================================================================
# Direct measurements of gs at Ca_ref=370 ppm, VPD~=1.5 kPa, well-watered.
# These anchor the absolute magnitude of gs in the model (see section 2.3).
# gs is the more directly measurable field quantity than the derived Ci/Ca.
GS_FC_OBS_SOY   = 0.55   # Ainsworth & Long (2005) SoyFACE meta-analysis.
                         # Yields Ci/Ca ~= 0.77 in this no-gm photo model
                         # vs observed ~0.73 -- small offset reflects absence
                         # of mesophyll conductance in the photo model.
GS_FC_OBS_MAIZE = 0.30   # Leakey et al. (2006b) Table II midday peak-season
                         # mean (DOY 190, 201, 215). Yields Ci/Ca ~= 0.37.

# Medlyn USO slope, from Lin et al. (2015, Nat. Clim. Change)
# PFT-level fits to a global database of 314 species across 56 field studies.
# Values from De Kauwe et al. (2015, Tab. 1) using Lin's database with the
# standard fitting protocol:
#     C3 cropland   g1 = 5.79 +/- 0.64
#     C4 grassland  g1 = 1.62 +/- 0.13   (no separate C4 cropland PFT in
#                                        Lin's database; C4 grassland is
#                                        the closest match for maize)
G1_SOY   = 5.79  # kPa^0.5  Lin et al. 2015 C3 cropland
G1_MAIZE = 1.62  # kPa^0.5  Lin et al. 2015 C4 grassland
G0       = 0.0   # residual conductance

# ===========================================================================
# Soil water stress (Verhoef-Egea linear beta)
# ===========================================================================
# Soil hydraulic anchors for Corn Belt silt loam (Mollisol). Values from
# pedotransfer functions for silt loam texture (Saxton & Rawls 2006 SSSAJ
# 70:1569; Allen et al. 1998 FAO-56) and consistent with field observations
# of well-watered VWC at SoyFACE/maize-FACE sites.
SWC_WP = 0.10   # permanent wilting point (m3 m^-3)
SWC_FC = 0.31   # field capacity (m3 m^-3)

# ===========================================================================
# Drought biochemistry parameters
# ===========================================================================
F_MIN_C3   = 0.60   # Soy Vcmax,Jmax retained at beta=0. Conceptual basis:
                    # Flexas & Medrano 2002 Ann. Bot. (phased
                    # drought response, severe-stress metabolic limitation);
# Maize Vpmax: direct beta scaling (Ripley et al. 2007 JXB 58:1351;
# Ghannoum 2009 Ann. Bot.; Bellasio 2024), no f_min

# ===========================================================================
# Default species parameters
# ===========================================================================
SOY_PARAMS   = dict(Vcmax=125.0, Jmax=200.0, Rd=1.5)
MAIZE_PARAMS = dict(Vcmax=49.0, Vpmax=120.0, Kp=80.0, Jmax=400.0)

# ===========================================================================
# Species colors (used by plotting scripts)
# ===========================================================================
SOY_COLOR   = '#2e8b57'
MAIZE_COLOR = '#d97706'
DIFF_COLOR  = '#333333'

# ===========================================================================
# Monte Carlo settings
# ===========================================================================
N_MC = 200
RNG  = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# C3 photosynthesis model (Farquhar et al. 1980, Bernacchi et al. 2001/2002)
# ---------------------------------------------------------------------------

def c3_photosynthesis(Ci, Vcmax=125.0, Jmax=200.0, Rd=1.5,
                      Kc=404.9, Ko=278.4, O=210.0, Gamma_star=42.75,
                      light_limited=True):
    """
    Net CO2 assimilation A for a C3 plant given intercellular CO2 (Ci).

    Defaults are calibrated to reproduce the soybean A/Ci curve in
    Ainsworth, Sanz-Saez & Leisner (2025) Phil. Trans. R. Soc. B Figure 2,
    which uses SoyFACE 2023 measurements at 27.5 C. The Ainsworth-reported
    biochemical Vcmax (194.6) and Jmax (260.0) include implicit mesophyll-
    conductance limitation; the effective values for a no-gm Farquhar model
    that reproduce their observed A/Ci shape are Vcmax=125, Jmax=200.

    Kc = 404.9 umol/mol, Ko = 278.4 mmol/mol, Gamma_star = 42.75 umol/mol
    following Bernacchi et al. 2001, 2002 -- Standard 25 C in vivo kinetics

    Rd = 1.5 umol m^-2 s^-1 (Standard dark respiration)
    """
    Ci = np.asarray(Ci, dtype=float)

    Ac = Vcmax * (Ci - Gamma_star) / (Ci + Kc * (1 + O / Ko))

    if not light_limited:
        return Ac - Rd

    J = Jmax  # saturating light assumed (standard A-Ci curve protocol);
              # at field PAR ~1800 umol/m2/s, J ~= Jmax via the
              # non-rectangular hyperbola (Farquhar & Wong 1984)
    Aj = J * (Ci - Gamma_star) / (4 * Ci + 8 * Gamma_star)

    A_gross = np.minimum(Ac, Aj)
    A_net = A_gross - Rd
    return A_net


# ---------------------------------------------------------------------------
# C4 photosynthesis model (von Caemmerer 2000, simplified)
# ---------------------------------------------------------------------------

def c4_photosynthesis(Ci, Vcmax=49.0, Vpmax=120.0, Kp=80.0,
                      Vpr=80.0, Rd=1.0, gbs=0.003, Jmax=400.0,
                      x=0.4, alpha=0.0, Kc=650.0, Ko=450.0, O=210.0,
                      gamma_star=1.93e-4):
    """
    Net CO2 assimilation A for a C4 plant given intercellular CO2 (Ci),
    using the von Caemmerer (2000) simplified model.

    Net A is the minimum of:
        Ac : enzyme-limited (PEPC and Rubisco co-limited)
        Aj : light/electron-transport-limited

    Defaults match Ainsworth, Sanz-Saez & Leisner (2025) Figure 2 for
    SoyFACE 2024 maize at 28 C: Vpmax = 120, plateau Asat ~= 48
    (Vcmax = 49 in this no-gm model gives plateau ~= 48).

    Kc=650 umol/mol, Ko=450 mmol/mol (von Caemmerer 2000)

    gamma_star is DIMENSIONLESS (von Caemmerer et al. 1994), characterizing
    the bundle-sheath CO2 compensation point. Os is taken as ambient O2 in
    umol/mol (matching CO2 units), valid when alpha=0.
    """
    Ci = np.asarray(Ci, dtype=float)

    # Bundle-sheath O2 in umol/mol (matching CO2 units); with alpha=0, Os ~= O.
    Os_umol = O * 1000.0  # mmol/mol -> umol/mol
    Rm = 0.5 * Rd  # mesophyll respiration (von Caemmerer 2000)

    # ---- Enzyme-limited rate (Ac), von Caemmerer 2000 eq. 4.21 ----
    Vp = np.minimum(Ci * Vpmax / (Ci + Kp), Vpr)

    a = 1.0 - (alpha * Kc) / (0.047 * Ko)
    b = -((Vp - Rm + gbs * Ci) +
          (Vcmax - Rd) +
          gbs * (Kc * (1.0 + O / Ko)) +
          (alpha / 0.047) * (gamma_star * Vcmax + Rd * Kc / Ko))
    c = ((Vcmax - Rd) * (Vp - Rm + gbs * Ci) -
         (Vcmax * gbs * gamma_star * Os_umol +
          Rd * gbs * Kc * (1.0 + O / Ko)))

    discriminant = b**2 - 4.0 * a * c
    if np.any(discriminant < 0):
        bad_idx = np.where(np.atleast_1d(discriminant) < 0)[0]
        raise ValueError(
            f"Negative discriminant in C4 enzyme-limited quadratic at "
            f"indices {bad_idx}. Parameter set may be non-physical."
        )
    Ac = (-b - np.sqrt(discriminant)) / (2.0 * a)

    # ---- Light-limited rate (Aj) ----
    J = Jmax  # saturating light assumed (standard A-Ci curve protocol)
    Vp_j = x * J / 2.0
    Vc_j = (1.0 - x) * J / 3.0

    a_j = 1.0 - (7.0 * gamma_star / 3.0) * (alpha / 0.047)
    b_j = -((Vp_j - Rm + gbs * Ci) +
            (Vc_j - Rd) +
            gbs * (7.0 / 3.0) * gamma_star * Os_umol +
            (alpha / 0.047) * (gamma_star * Vc_j + Rd * Kc / Ko))
    c_j = ((Vc_j - Rd) * (Vp_j - Rm + gbs * Ci) -
           (Vc_j * gbs * gamma_star * Os_umol +
            Rd * gbs * (7.0 / 3.0) * gamma_star * Os_umol))

    discriminant_j = np.maximum(b_j**2 - 4.0 * a_j * c_j, 0.0)
    Aj = (-b_j - np.sqrt(discriminant_j)) / (2.0 * a_j)

    A_net = np.minimum(Ac, Aj)
    return A_net


# ---------------------------------------------------------------------------
# Operating point: solve for Ci given Ca and stomatal conductance
# ---------------------------------------------------------------------------

def solve_operating_point(Ca, gs, photo_func, **kwargs):
    """
    Solve for operating Ci, A given atmospheric CO2 (Ca) and stomatal
    conductance gs.

    Supply: A = gs * (Ca - Ci) / 1.6
    Demand: A = photo_func(Ci, **kwargs)

    Raises ValueError on bracket failure rather than returning silently.
    """
    def residual(Ci_val):
        A_demand = photo_func(np.array([Ci_val]), **kwargs)[0]
        A_supply = gs * (Ca - Ci_val) / 1.6
        return A_demand - A_supply

    lo, hi = 1.0, Ca - 0.1
    f_lo, f_hi = residual(lo), residual(hi)
    if f_lo * f_hi > 0:
        scan = np.linspace(lo, hi, 50)
        vals = [residual(v) for v in scan]
        sign_changes = [i for i in range(len(vals) - 1)
                        if vals[i] * vals[i + 1] < 0]
        if not sign_changes:
            raise ValueError(
                f"solve_operating_point: no sign change in residual on "
                f"[{lo}, {hi}] for Ca={Ca}, gs={gs}. "
                f"residual({lo})={f_lo:.3f}, residual({hi})={f_hi:.3f}."
            )
        i = sign_changes[0]
        lo, hi = scan[i], scan[i + 1]

    Ci_op = brentq(residual, lo, hi)
    A_op = photo_func(np.array([Ci_op]), **kwargs)[0]
    return Ci_op, A_op


# ===========================================================================
# Stress factors
# ===========================================================================

def beta_swc(swc):
    """Verhoef-Egea linear soil-water stress factor, clipped to [0, 1]."""
    return float(np.clip((swc - SWC_WP) / (SWC_FC - SWC_WP), 0.0, 1.0))


def soy_biochem_stressed(params, beta):
    """C3 biochemistry under drought -- Vcmax,Jmax linear beta-ramp.

    scale = f_min + (1 - f_min) x beta
    beta=1 -> scale=1.0 (no stress)
    beta=0 -> scale=f_min=0.60 (severe stress retains 60% capacity)
    """
    scale = F_MIN_C3 + (1.0 - F_MIN_C3) * beta
    return {**params,
            'Vcmax': params['Vcmax'] * scale,
            'Jmax':  params['Jmax']  * scale}


def maize_biochem_stressed(params, beta):
    """C4 biochemistry under drought -- Vpmax beta-scaled.

    Direct linear scaling: beta=0 -> CCM fully shut down.
    Vcmax (Rubisco) and Jmax not scaled -- bundle-sheath protects Rubisco
    from cytosolic dehydration in C4 anatomy.
    """
    return {**params, 'Vpmax': params['Vpmax'] * beta}


# ===========================================================================
# Stomatal conductance -- anchored Medlyn USO model
# ===========================================================================

def _medlyn_gs_at_ca(Ca, vpd, g1, photo_func, **params):
    """Pure Medlyn gs at the analytical Ci/Ca = g1/(g1+sqrt(VPD)) point.

    Used to compute the dimensionless Ca-dependence kernel for soybean.
    Returns gs in mol H2O m^-2 s^-1.
    """
    r  = g1 / (g1 + np.sqrt(vpd))
    Ci = Ca * r
    A  = max(photo_func(np.array([Ci]), **params)[0], 0.0)
    return G0 + 1.6 * (1.0 + g1 / np.sqrt(vpd)) * A / Ca


def medlyn_ca_kernel(Ca, vpd, g1, photo_func, **params):
    """Dimensionless Medlyn USO Ca-kernel: gs(Ca) / gs(Ca_ref).

    Captures the Ca-dependent gs reduction predicted by USO theory holding
    the Medlyn slope g1 fixed.
    """
    gs_ca  = _medlyn_gs_at_ca(Ca,     vpd, g1, photo_func, **params)
    gs_ref = _medlyn_gs_at_ca(CA_REF, VPD_REF, g1, photo_func, **params)
    if gs_ref < 1e-9:
        return 1.0
    return gs_ca / gs_ref


def gs_soy(Ca, vpd, swc, params=None):
    """C3 soybean stomatal conductance -- anchored Medlyn.

    gs = beta(SWC) x gs_fc_obs_soy x kernel_medlyn_C3(Ca, VPD)

    The Medlyn kernel is computed using the leaf's current biochemical state
    (stressed or unstressed depending on what params is passed). This is
    self-consistent with the demand curve: a stressed leaf's stomata track
    A/Ca at the actual reduced Vcmax/Jmax, so the kernel reflects both the
    Ca-dependence and the drought-modified A-Ci curve geometry. beta scales
    the absolute gs magnitude independently via the hydraulic signal.
    """
    if params is None:
        params = SOY_PARAMS
    beta   = beta_swc(swc)
    kernel = medlyn_ca_kernel(Ca, vpd, G1_SOY, c3_photosynthesis, **params)
    return beta * GS_FC_OBS_SOY * kernel


def gs_maize(Ca, vpd, swc, params=None):
    """C4 maize stomatal conductance -- anchored Medlyn.

    gs = beta(SWC) x gs_fc_obs_maize x kernel_medlyn_C4(Ca, VPD, g1_maize=1.62)

    The Medlyn kernel is computed using the leaf's current biochemical state
    (stressed or unstressed depending on what params is passed). Under
    drought, reduced Vpmax shifts the C4 saturation knee leftward, making
    A more Ca-sensitive; the kernel captures this change, which is
    self-consistent with the demand curve used in the operating-point solve.
    beta scales the absolute gs magnitude independently via the hydraulic signal.
    """
    if params is None:
        params = MAIZE_PARAMS
    beta   = beta_swc(swc)
    kernel = medlyn_ca_kernel(Ca, vpd, G1_MAIZE, c4_photosynthesis, **params)
    return beta * GS_FC_OBS_MAIZE * kernel


# ===========================================================================
# CFE computation
# ===========================================================================

def _safe_op(Ca, gs, photo_func, **params):
    """Robust operating-point solve -- returns None on bracket failure."""
    try:
        return solve_operating_point(Ca, gs, photo_func, **params)
    except Exception:
        return None, None


def cfe_soy(swc, vpd, params):
    """Per-ppm CFE for soybean across CA_LOW -> CA_HIGH."""
    beta     = beta_swc(swc)
    stressed = soy_biochem_stressed(params, beta)
    gs_lo    = gs_soy(CA_LOW,  vpd, swc, params=stressed)
    gs_hi    = gs_soy(CA_HIGH, vpd, swc, params=stressed)
    _, A_lo = _safe_op(CA_LOW,  gs_lo, c3_photosynthesis, **stressed)
    _, A_hi = _safe_op(CA_HIGH, gs_hi, c3_photosynthesis, **stressed)
    if A_lo is None or A_lo <= 0:
        return None
    return 100.0 * (A_hi - A_lo) / A_lo / (CA_HIGH - CA_LOW)


def cfe_maize(swc, vpd, params):
    """Per-ppm CFE for maize across CA_LOW -> CA_HIGH."""
    beta     = beta_swc(swc)
    stressed = maize_biochem_stressed(params, beta)
    gs_lo    = gs_maize(CA_LOW,  vpd, swc, params=stressed)
    gs_hi    = gs_maize(CA_HIGH, vpd, swc, params=stressed)
    _, A_lo = _safe_op(CA_LOW,  gs_lo, c4_photosynthesis, **stressed)
    _, A_hi = _safe_op(CA_HIGH, gs_hi, c4_photosynthesis, **stressed)
    if A_lo is None or A_lo <= 0:
        return None
    return 100.0 * (A_hi - A_lo) / A_lo / (CA_HIGH - CA_LOW)


def compute_dcfe_vs_swc(swc_values, vpd, soy_mc, maize_mc):
    """DCFE = CFE_soy - CFE_maize across SWC, with paired-MC bands."""
    n_mc, n_swc = len(soy_mc), len(swc_values)
    cfe_s = np.full((n_mc, n_swc), np.nan)
    cfe_m = np.full((n_mc, n_swc), np.nan)
    for j, swc in enumerate(swc_values):
        for i, p in enumerate(soy_mc):
            v = cfe_soy(swc, vpd, p)
            if v is not None:
                cfe_s[i, j] = v
        for i, p in enumerate(maize_mc):
            v = cfe_maize(swc, vpd, p)
            if v is not None:
                cfe_m[i, j] = v
    diff = cfe_s - cfe_m

    def band(arr):
        return (np.nanmedian(arr,    axis=0),
                np.nanpercentile(arr, 10, axis=0),
                np.nanpercentile(arr, 90, axis=0))

    return band(cfe_s), band(cfe_m), band(diff)


# ===========================================================================
# Monte Carlo parameter sampling
# ===========================================================================
#
# Cultivar variation in C3 soybean is represented by sampling biochemical
# parameters from biologically realistic ranges. The key design choice:
# parameterize Jmax through J/V (the carboxylation-to-electron-transport
# balance) rather than as an independent uniform variate.
#
# Vcmax (capacity) and J/V (balance) are conceptually independent biological
# factors -- different cultivars and environments shift each independently.
# Sampling Vcmax and Jmax as independent uniform variates allows unphysical
# combinations (e.g. J/V = 1.2, which does not occur in healthy soybean):
# such samples produce A-Ci curves that plateau at very low Ci, generating
# a long lower tail in the CFE distribution that is a sampling artifact, not
# a physiological feature.
#
# Sampling scheme:
#     Vcmax ~ U[105, 145]
#     J/V   ~ U[1.45, 1.75]    (centred on 1.60)
#     Jmax  = (J/V) x Vcmax
#
# The J/V range covers documented soybean variation:
#   Bernacchi et al. (2005) Planta 220:434   -- J/V ~= 1.66 (SoyFACE in vivo)
#   Rogers et al. (2004)    PCE 27:449       -- seasonal J/V = 1.6-1.8 (SoyFACE)
#   Bishop et al. (2015)    PCE 38:1765      -- 18 cultivars, J/V ~= 1.5-1.7
#
# For maize, Vcmax / Vpmax / Jmax are sampled independently because (i) the
# C4 J/V relationship is less well-constrained than for C3 species, and
# (ii) the C4 plateau over the Ci range of interest is governed by Vpmax,
# not by the J/V balance.

JV_RANGE_SOY = (1.45, 1.75)


def sample_soy_params(n):
    """Sample C3 soybean params via independent Vcmax and J/V."""
    out = []
    for _ in range(n):
        v  = RNG.uniform(105, 145)
        jv = RNG.uniform(*JV_RANGE_SOY)
        out.append(dict(Vcmax=v, Jmax=jv * v, Rd=RNG.uniform(1.0, 2.0)))
    return out


def sample_maize_params(n):
    """Sample C4 maize params (independent -- see module note)."""
    return [dict(Vcmax=RNG.uniform(40, 55),
                 Vpmax=RNG.uniform(95, 135),
                 Kp=RNG.uniform(60, 100),
                 Jmax=RNG.uniform(350, 450)) for _ in range(n)]


def get_curve_band(photo_func, params_list, ci_max=700, n=300):
    """Return (Ci, median, p10, p90) across a list of parameter dicts."""
    Ci    = np.linspace(1, ci_max, n)
    A_all = np.array([photo_func(Ci, **p) for p in params_list])
    return (Ci,
            np.median(A_all,        axis=0),
            np.percentile(A_all, 10, axis=0),
            np.percentile(A_all, 90, axis=0))


# ===========================================================================
# Diagnostic table -- print operating points and CFE for verification
# ===========================================================================

def print_diagnostics(vpd=1.5):
    """Print a verification table -- should match observations."""
    print()
    print("=" * 78)
    print("DIAGNOSTIC TABLE -- operating points at VPD = {:.1f} kPa, Ca = {} ppm"
          .format(vpd, CA_REF))
    print("=" * 78)
    print(f"  {'Species':<8} {'SWC':>5}  {'beta':>5}  {'gs':>6}  {'Ci':>6}  "
          f"{'Ci/Ca':>6}  {'A':>6}  CFE(336->420)")
    print("  " + "-" * 76)

    for species, gs_func, photo_func, params_default, biochem_func, cfe_func in [
        ('Soy',   gs_soy,   c3_photosynthesis, SOY_PARAMS,   soy_biochem_stressed,   cfe_soy),
        ('Maize', gs_maize, c4_photosynthesis, MAIZE_PARAMS, maize_biochem_stressed, cfe_maize),
    ]:
        for swc in [0.31, 0.20, 0.15, 0.12]:
            beta = beta_swc(swc)
            stressed = biochem_func(params_default, beta)
            gs = gs_func(CA_REF, vpd, swc, params=stressed)
            try:
                Ci, A = solve_operating_point(CA_REF, gs, photo_func, **stressed)
            except Exception:
                print(f"  {species:<8} {swc:>5.2f}  {beta:>5.2f}  ERROR")
                continue
            cfe = cfe_func(swc, vpd, params_default)
            cfe_str = f"{cfe:+.4f}" if cfe is not None else "    N/A"
            print(f"  {species:<8} {swc:>5.2f}  {beta:>5.2f}  "
                  f"{gs:>6.3f}  {Ci:>6.1f}  {Ci/CA_REF:>6.3f}  "
                  f"{A:>6.2f}  {cfe_str} %/ppm")
        print()


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

def _run_tests():
    print("Running unit tests...")

    # Test 1: C3 at Gamma_star -> A ~ -Rd
    A_c3 = c3_photosynthesis(np.array([42.75]))[0]
    assert abs(A_c3 - (-1.5)) < 0.5, f"C3 at Gamma*: got {A_c3:.3f}"

    # Test 2: C3 at Ci=400 in saturating range
    A_c3_high = c3_photosynthesis(np.array([400]))[0]
    assert 20 < A_c3_high < 40, f"C3 A(Ci=400)={A_c3_high:.2f}"

    # Test 3: C4 nearly saturated at Ci=150 (canonical knee ~120)
    A_c4_150 = c4_photosynthesis(np.array([150]))[0]
    A_c4_high = c4_photosynthesis(np.array([400]))[0]
    ratio = A_c4_150 / A_c4_high
    assert ratio > 0.9, f"C4 A(150)/A(400)={ratio:.3f}"

    # Test 4: C4 plateau in field range (matches Ainsworth Fig 2: Asat~=48)
    assert 40 < A_c4_high < 55, f"C4 plateau={A_c4_high:.2f}"

    # Test 5: C3 monotonic
    Ci_test = np.array([200, 250, 300, 350])
    A_test = c3_photosynthesis(Ci_test)
    assert np.all(np.diff(A_test) > 0), "C3 monotonicity"

    # Test 6: Wet-maize operating point physical
    Ci_op, A_op = solve_operating_point(Ca=400, gs=0.25,
                                        photo_func=c4_photosynthesis)
    assert 100 < Ci_op < 250, f"Wet maize Ci={Ci_op:.1f}"
    assert 30 < A_op < 50, f"Wet maize A={A_op:.1f}"

    # Test 7: Ci/Ca for wet maize in canonical range
    ratio_ci_ca = Ci_op / 400
    assert 0.25 < ratio_ci_ca < 0.65, f"Ci/Ca={ratio_ci_ca:.2f}"

    print(f"  C3 A(42.75)        = {A_c3:.2f}   (target ~= -Rd = -1.5)")
    print(f"  C3 A(400)          = {A_c3_high:.2f}")
    print(f"  C4 A(150)          = {A_c4_150:.2f}")
    print(f"  C4 A(400)          = {A_c4_high:.2f}")
    print(f"  C4 A(150)/A(400)   = {ratio:.3f}   (target > 0.9)")
    print(f"  Wet-maize Ci/Ca    = {ratio_ci_ca:.2f}  (target ~= 0.3-0.5)")
    print("All tests passed.")


if __name__ == "__main__":
    _run_tests()
