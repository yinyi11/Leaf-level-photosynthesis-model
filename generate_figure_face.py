"""FACE-level diagnostic figure: gs, A, and delta-A vs SWC at 376->550 ppm."""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from photosynthesis_models import (
    # photosynthesis / operating point
    c3_photosynthesis, c4_photosynthesis, solve_operating_point,
    # constants
    SWC_FC,
    SOY_PARAMS, MAIZE_PARAMS,
    SOY_COLOR, MAIZE_COLOR, DIFF_COLOR,
    N_MC,
    # model functions
    beta_swc, soy_biochem_stressed, maize_biochem_stressed,
    gs_soy, gs_maize,
    sample_soy_params, sample_maize_params,
)

# ===========================================================================
# Configuration
# ===========================================================================
CA_AMB = 376    # SoyFACE / maize-FACE ambient
CA_ELE = 550    # SoyFACE / maize-FACE elevated
VPD    = 1.5    # kPa, midday FACE conditions

# Soil moisture sweep (panels a, b, c)
SWC_VALUES = np.linspace(0.13, 0.32, 60)

# SWC anchors for vertical reference lines
SWC_DROUGHT = 0.18         # the "moderate drought" reference

# ===========================================================================
# Panel (a) plotting constants
# ===========================================================================
PANEL_A_SWC       = [0.31, 0.18]   # well-watered, moderate drought
PANEL_A_SWC_LS    = ['-',  '--']
PANEL_A_VPD       = [1.0,  2.5]    # low, high VPD
PANEL_A_VPD_ALPHA = [1.0,  0.45]


# ===========================================================================
# Compute responses across SWC
# ===========================================================================

def compute_responses(swc_values, ca, params_default_soy, params_default_maize):
    """Return arrays of gs, A for both species at the given Ca across SWC."""
    n = len(swc_values)
    gs_s = np.zeros(n); A_s = np.zeros(n)
    gs_m = np.zeros(n); A_m = np.zeros(n)
    for i, swc in enumerate(swc_values):
        beta = beta_swc(swc)
        soy_str = soy_biochem_stressed(params_default_soy, beta)
        mz_str  = maize_biochem_stressed(params_default_maize, beta)
        gs_s[i] = gs_soy(ca, VPD, swc, params=soy_str)
        gs_m[i] = gs_maize(ca, VPD, swc, params=mz_str)
        try:
            _, A_s[i] = solve_operating_point(
                ca, gs_s[i], c3_photosynthesis, **soy_str)
        except Exception:
            A_s[i] = np.nan
        try:
            _, A_m[i] = solve_operating_point(
                ca, gs_m[i], c4_photosynthesis, **mz_str)
        except Exception:
            A_m[i] = np.nan
    return gs_s, A_s, gs_m, A_m


def compute_mc_band(swc_values, ca_amb, ca_ele, soy_mc, maize_mc, q_lo=10, q_hi=90):
    """Monte-Carlo bands for delta-A% across SWC under FACE treatment."""
    n_swc = len(swc_values)
    n_mc = len(soy_mc)
    dA_soy = np.full((n_mc, n_swc), np.nan)
    dA_mz  = np.full((n_mc, n_swc), np.nan)

    for j, p_s in enumerate(soy_mc):
        for i, swc in enumerate(swc_values):
            beta = beta_swc(swc)
            ps = soy_biochem_stressed(p_s, beta)
            try:
                gs_a = gs_soy(ca_amb, VPD, swc, params=ps)
                gs_e = gs_soy(ca_ele, VPD, swc, params=ps)
                _, A_a = solve_operating_point(ca_amb, gs_a, c3_photosynthesis, **ps)
                _, A_e = solve_operating_point(ca_ele, gs_e, c3_photosynthesis, **ps)
                if A_a > 0:
                    dA_soy[j, i] = 100.0 * (A_e - A_a) / A_a
            except Exception:
                pass

    for j, p_m in enumerate(maize_mc):
        for i, swc in enumerate(swc_values):
            beta = beta_swc(swc)
            pm = maize_biochem_stressed(p_m, beta)
            try:
                gs_a = gs_maize(ca_amb, VPD, swc, params=pm)
                gs_e = gs_maize(ca_ele, VPD, swc, params=pm)
                _, A_a = solve_operating_point(ca_amb, gs_a, c4_photosynthesis, **pm)
                _, A_e = solve_operating_point(ca_ele, gs_e, c4_photosynthesis, **pm)
                if A_a > 0:
                    dA_mz[j, i] = 100.0 * (A_e - A_a) / A_a
            except Exception:
                pass

    soy_med = np.nanmedian(dA_soy, axis=0)
    soy_lo  = np.nanpercentile(dA_soy, q_lo, axis=0)
    soy_hi  = np.nanpercentile(dA_soy, q_hi, axis=0)
    mz_med  = np.nanmedian(dA_mz, axis=0)
    mz_lo   = np.nanpercentile(dA_mz, q_lo, axis=0)
    mz_hi   = np.nanpercentile(dA_mz, q_hi, axis=0)

    return (soy_med, soy_lo, soy_hi), (mz_med, mz_lo, mz_hi)


# ===========================================================================
# Panel (a): gs vs Ca at a few VPD x SWC combinations
# ===========================================================================

def plot_gs_vs_ca(ax):
    ca_sweep = np.linspace(280, 600, 100)

    for swc, ls in zip(PANEL_A_SWC, PANEL_A_SWC_LS):
        for vpd, alpha in zip(PANEL_A_VPD, PANEL_A_VPD_ALPHA):
            beta      = beta_swc(swc)
            soy_str   = soy_biochem_stressed(SOY_PARAMS,    beta)
            maize_str = maize_biochem_stressed(MAIZE_PARAMS, beta)

            gs_s = np.array([gs_soy(Ca,   vpd, swc, params=soy_str)
                             for Ca in ca_sweep])
            gs_m = np.array([gs_maize(Ca, vpd, swc, params=maize_str)
                             for Ca in ca_sweep])

            ax.plot(ca_sweep, gs_s, color=SOY_COLOR,   linewidth=1.8,
                    linestyle=ls, alpha=alpha)
            ax.plot(ca_sweep, gs_m, color=MAIZE_COLOR, linewidth=1.8,
                    linestyle=ls, alpha=alpha)

    # Mark the FACE Ca treatment levels
    for Ca, label in [(CA_AMB, f'Ambient\n({CA_AMB})'),
                      (CA_ELE, f'Elevated\n({CA_ELE})')]:
        ax.axvline(Ca, color='gray', linewidth=0.8, linestyle=':', alpha=0.6)
        ax.text(Ca + 4, 0.005, label, fontsize=6.5, color='gray',
                va='bottom', ha='left')

    ax.set_xlabel(r'Atmospheric CO$_2$, $C_a$ ($\mu$mol mol$^{-1}$)', fontsize=10)
    ax.set_ylabel(r'Stomatal conductance, $g_s$'
                  '\n' r'(mol H$_2$O m$^{-2}$ s$^{-1}$)', fontsize=10)
    ax.set_xlim(280, 600)
    ax.set_ylim(bottom=0)

    handles = [
        Line2D([], [], color=SOY_COLOR,   linewidth=2.0, label='C3 Soybean'),
        Line2D([], [], color=MAIZE_COLOR, linewidth=2.0, label='C4 Maize'),
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label='Soil moisture'),
    ] + [
        Line2D([], [], color='gray', linewidth=1.4, linestyle=ls,
               label=f'  SWC = {swc:.2f}')
        for swc, ls in zip(PANEL_A_SWC, PANEL_A_SWC_LS)
    ] + [
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label='VPD'),
    ] + [
        Line2D([], [], color='gray', linewidth=1.4, alpha=a,
               label=f'  {vpd:.1f} kPa')
        for vpd, a in zip(PANEL_A_VPD, PANEL_A_VPD_ALPHA)
    ]
    leg = ax.legend(handles=handles, loc='upper right', frameon=False,
                    fontsize=7.5, handlelength=1.8, handletextpad=0.4,
                    labelspacing=0.15, bbox_to_anchor=(0.98, 0.98))
    for text in leg.get_texts():
        if text.get_text() in ('Soil moisture', 'VPD'):
            text.set_fontweight('bold')


# ===========================================================================
# Panel (b): A vs SWC at Ca=376 vs Ca=550
# ===========================================================================

def plot_A(ax, swc_values, A_s_amb, A_s_ele, A_m_amb, A_m_ele):
    ax.plot(swc_values, A_s_amb, color=SOY_COLOR,   linewidth=2.0, linestyle='-')
    ax.plot(swc_values, A_s_ele, color=SOY_COLOR,   linewidth=2.0, linestyle='--')
    ax.plot(swc_values, A_m_amb, color=MAIZE_COLOR, linewidth=2.0, linestyle='-')
    ax.plot(swc_values, A_m_ele, color=MAIZE_COLOR, linewidth=2.0, linestyle='--')

    # Fill between the two Ca curves to make the CFE gap visible
    ax.fill_between(swc_values, A_s_amb, A_s_ele,
                    color=SOY_COLOR,   alpha=0.12, linewidth=0)
    ax.fill_between(swc_values, A_m_amb, A_m_ele,
                    color=MAIZE_COLOR, alpha=0.12, linewidth=0)

    ax.set_xlabel(r'Root-zone soil water content (m$^{3}$ m$^{-3}$)', fontsize=10)
    ax.set_ylabel(r'Net photosynthesis, $A$ ($\mu$mol m$^{-2}$ s$^{-1}$)',
                  fontsize=10)
    ax.set_xlim(0.13, 0.33)
    ax.set_ylim(0, 50)

    handles = [
        Line2D([], [], color=SOY_COLOR,   linewidth=2.0, label='C3 Soybean'),
        Line2D([], [], color=MAIZE_COLOR, linewidth=2.0, label='C4 Maize'),
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label=r'Atmospheric CO$_2$'),
        Line2D([], [], color='gray', linewidth=1.8, linestyle='-',
               label=f'  Ambient ({CA_AMB} ppm)'),
        Line2D([], [], color='gray', linewidth=1.8, linestyle='--',
               label=f'  Elevated ({CA_ELE} ppm)'),
    ]
    leg = ax.legend(handles=handles, loc='upper left', frameon=False,
                    fontsize=7.5, handlelength=1.8, handletextpad=0.4,
                    labelspacing=0.15)
    for text in leg.get_texts():
        if text.get_text() == r'Atmospheric CO$_2$':
            text.set_fontweight('bold')


# ===========================================================================
# Panel (c): delta-A % (FACE CFE response) vs SWC, with MC bands
# ===========================================================================

def plot_dA(ax, swc_values, soy_band, mz_band):
    soy_med, soy_lo, soy_hi = soy_band
    mz_med,  mz_lo,  mz_hi  = mz_band

    ax.fill_between(swc_values, soy_lo, soy_hi, color=SOY_COLOR,
                    alpha=0.18, linewidth=0)
    ax.fill_between(swc_values, mz_lo,  mz_hi,  color=MAIZE_COLOR,
                    alpha=0.18, linewidth=0)

    ax.plot(swc_values, soy_med, color=SOY_COLOR,   linewidth=2.2,
            label=r'C3 Soybean')
    ax.plot(swc_values, mz_med,  color=MAIZE_COLOR, linewidth=2.2,
            label=r'C4 Maize')

    ax.axhline(0, color='gray', linewidth=0.7, linestyle='--', alpha=0.5)

    # Observed FACE values at well-watered conditions (SWC ~= FC = 0.31)
    # Soybean: +15-25% (Ainsworth & Long 2005 meta-analysis; Long et al. 2006)
    #          plotted as midpoint +/-range error bar
    # Maize:   ~0%, not significant (Leakey et al. 2006b Table II)
    obs_swc = SWC_FC   # well-watered anchor
    ax.errorbar(obs_swc, 20, yerr=[[5], [5]], fmt='o',
                color=SOY_COLOR, markersize=6, capsize=4, linewidth=1.5,
                zorder=5, label='Obs. soy (A&L 2005)')
    ax.errorbar(obs_swc, 0, yerr=[[0], [3]], fmt='s',
                color=MAIZE_COLOR, markersize=6, capsize=4, linewidth=1.5,
                zorder=5, label='Obs. maize (Leakey 2006b)')

    ax.set_xlabel(r'Root-zone soil water content (m$^{3}$ m$^{-3}$)', fontsize=10)
    ax.set_ylabel(r'$\Delta A$ at FACE elevated CO$_2$ (%)', fontsize=10)
    ax.set_xlim(0.13, 0.33)
    ax.set_ylim(-2, 50)

    ax.legend(loc='upper right', fontsize=8, frameon=False,
              handlelength=2.0, bbox_to_anchor=(0.98, 0.98))
    ax.text(0.03, 0.97,
            f'Ca: {CA_AMB} -> {CA_ELE} ppm\nVPD = {VPD:.1f} kPa',
            transform=ax.transAxes, fontsize=8, ha='left', va='top',
            color='gray')


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("Computing FACE-level responses across SWC ...")
    _, A_s_amb, _, A_m_amb = compute_responses(
        SWC_VALUES, CA_AMB, SOY_PARAMS, MAIZE_PARAMS)
    _, A_s_ele, _, A_m_ele = compute_responses(
        SWC_VALUES, CA_ELE, SOY_PARAMS, MAIZE_PARAMS)

    print("Sampling Monte-Carlo cultivars and computing delta-A% bands ...")
    soy_mc   = sample_soy_params(N_MC)
    maize_mc = sample_maize_params(N_MC)
    soy_band, mz_band = compute_mc_band(
        SWC_VALUES, CA_AMB, CA_ELE, soy_mc, maize_mc)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))

    plot_gs_vs_ca(axes[0])
    plot_A (axes[1], SWC_VALUES, A_s_amb,  A_s_ele,  A_m_amb,  A_m_ele)
    plot_dA(axes[2], SWC_VALUES, soy_band, mz_band)

    for i, letter in enumerate(['(a)', '(b)', '(c)']):
        axes[i].text(-0.13, 1.02, letter, transform=axes[i].transAxes,
                     fontsize=12, fontweight='bold', va='bottom', ha='left')

    plt.tight_layout()
    out_path = os.path.join(_HERE, 'figure_face.png')
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {out_path}")

    # ---- Diagnostic table at key SWC values ----
    print()
    print("=" * 80)
    print(f"FACE diagnostic table (Ca: {CA_AMB} -> {CA_ELE} ppm, VPD={VPD:.1f} kPa)")
    print("=" * 80)
    print(f"  {'Species':<8} {'SWC':>5}  {'gs_a':>6} {'gs_e':>6} {'Dgs%':>6} | "
          f"{'A_a':>6} {'A_e':>6} {'DA%':>6}")
    print("  " + "-" * 76)
    for swc in [0.31, 0.27, 0.23, 0.20, 0.18, 0.15]:
        beta = beta_swc(swc)
        for label, photo, gs_func, params, biochem in [
            ('Soy',   c3_photosynthesis, gs_soy,   SOY_PARAMS,   soy_biochem_stressed),
            ('Maize', c4_photosynthesis, gs_maize, MAIZE_PARAMS, maize_biochem_stressed),
        ]:
            ps = biochem(params, beta)
            gs_a = gs_func(CA_AMB, VPD, swc, params=ps)
            gs_e = gs_func(CA_ELE, VPD, swc, params=ps)
            _, A_a = solve_operating_point(CA_AMB, gs_a, photo, **ps)
            _, A_e = solve_operating_point(CA_ELE, gs_e, photo, **ps)
            print(f"  {label:<8} {swc:>5.2f}  "
                  f"{gs_a:>6.3f} {gs_e:>6.3f} "
                  f"{100*(gs_e-gs_a)/gs_a:>+5.1f}% | "
                  f"{A_a:>6.2f} {A_e:>6.2f} "
                  f"{100*(A_e-A_a)/A_a:>+5.1f}%")
        print()


if __name__ == '__main__':
    main()
