"""Figure 4 -- A-Ci curves, DCFE vs SWC, and gs vs Ca for soybean and maize."""

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
    CA_LEVELS, CA_LABELS, CA_COLORS, CA_LOW, CA_HIGH, CA_REF,
    GS_FC_OBS_SOY, GS_FC_OBS_MAIZE,
    G1_SOY, G1_MAIZE, G0,
    SWC_WP, SWC_FC,
    F_MIN_C3,
    SOY_PARAMS, MAIZE_PARAMS,
    SOY_COLOR, MAIZE_COLOR, DIFF_COLOR,
    N_MC, RNG,
    # model functions
    beta_swc, soy_biochem_stressed, maize_biochem_stressed,
    _medlyn_gs_at_ca, medlyn_ca_kernel, gs_soy, gs_maize,
    _safe_op, cfe_soy, cfe_maize, compute_dcfe_vs_swc,
    sample_soy_params, sample_maize_params, get_curve_band,
    print_diagnostics,
)

# ===========================================================================
# Plot-only constants
# ===========================================================================
SWC_LEVELS = [0.31, 0.21]          # well-watered, moderate drought
SWC_COLORS = ['#4393c3', '#f4a582']
SWC_LS     = ['-', '--']
SWC_LABELS = ['Well-watered (0.31)',
              'Moderate drought (0.21)']

VPD_REF    = 1.5   # kPa -- for panel (a) operating points and (b) DCFE-vs-SWC (FACE midday)


PANEL_C_VPD    = [1.0, 2.5]
PANEL_C_VPD_LS = ['-', '--']
PANEL_C_SWC    = [0.31, 0.21]
PANEL_C_ALPHA  = [1.0, 0.45]


# ===========================================================================
# Panel (a): A-Ci curves with operating points
# ===========================================================================

def plot_aci(ax, soy_mc, maize_mc):
    Ci = np.linspace(1, 700, 300)

    # Well-watered demand curves with MC bands
    A_soy_all   = np.array([c3_photosynthesis(Ci, **p) for p in soy_mc])
    A_maize_all = np.array([c4_photosynthesis(Ci, **p) for p in maize_mc])
    A_soy_med   = np.median(A_soy_all,   axis=0)
    A_maize_med = np.median(A_maize_all, axis=0)

    ax.fill_between(Ci, np.percentile(A_soy_all,   10, axis=0),
                        np.percentile(A_soy_all,   90, axis=0),
                    color=SOY_COLOR,   alpha=0.13, linewidth=0, zorder=1)
    ax.fill_between(Ci, np.percentile(A_maize_all, 10, axis=0),
                        np.percentile(A_maize_all, 90, axis=0),
                    color=MAIZE_COLOR, alpha=0.13, linewidth=0, zorder=1)
    ax.plot(Ci, A_soy_med,   color=SOY_COLOR,   linewidth=2.4, zorder=3)
    ax.plot(Ci, A_maize_med, color=MAIZE_COLOR, linewidth=2.4, zorder=3)

    # Per-SWC stressed curves (median params), used to read off A at operating Ci
    swc_curves = {}
    for swc in SWC_LEVELS:
        beta = beta_swc(swc)
        soy_str   = soy_biochem_stressed(SOY_PARAMS,    beta)
        maize_str = maize_biochem_stressed(MAIZE_PARAMS, beta)
        A_s = c3_photosynthesis(Ci, **soy_str)   if beta < 1.0 else A_soy_med
        A_m = c4_photosynthesis(Ci, **maize_str) if beta < 1.0 else A_maize_med
        swc_curves[swc] = (interp1d(Ci, A_s, kind='linear', fill_value='extrapolate'),
                           interp1d(Ci, A_m, kind='linear', fill_value='extrapolate'),
                           A_s, A_m)

    # Plot stressed demand curves (only for SWC < FC)
    for i_swc, swc in enumerate(SWC_LEVELS):
        if beta_swc(swc) >= 1.0:
            continue
        _, _, A_s, A_m = swc_curves[swc]
        ax.plot(Ci, A_s, color=SOY_COLOR,   linewidth=1.2,
                linestyle=SWC_LS[i_swc], alpha=0.55, zorder=2)
        ax.plot(Ci, A_m, color=MAIZE_COLOR, linewidth=1.2,
                linestyle=SWC_LS[i_swc], alpha=0.55, zorder=2)

    # Operating points: supply lines + markers
    for i_ca, Ca in enumerate(CA_LEVELS):
        era_color = CA_COLORS[i_ca]
        ax.scatter([Ca], [0], s=18, color=era_color, zorder=4,
                   edgecolor='black', linewidth=0.5)

        for i_swc, swc in enumerate(SWC_LEVELS):
            beta         = beta_swc(swc)
            soy_str      = soy_biochem_stressed(SOY_PARAMS,    beta)
            maize_str    = maize_biochem_stressed(MAIZE_PARAMS, beta)
            gs_s         = gs_soy(Ca,   VPD_REF, swc, params=soy_str)
            gs_m         = gs_maize(Ca, VPD_REF, swc, params=maize_str)

            ci_s, _ = _safe_op(Ca, gs_s, c3_photosynthesis, **soy_str)
            ci_m, _ = _safe_op(Ca, gs_m, c4_photosynthesis, **maize_str)
            if ci_s is None or ci_m is None:
                continue

            soy_interp, maize_interp, _, _ = swc_curves[swc]
            a_s = float(soy_interp(ci_s))
            a_m = float(maize_interp(ci_m))

            ax.plot([Ca, ci_s], [0, a_s], color=SOY_COLOR, linewidth=0.7,
                    linestyle=SWC_LS[i_swc], alpha=0.7, zorder=2)
            ax.plot([Ca, ci_m], [0, a_m], color=MAIZE_COLOR, linewidth=0.7,
                    linestyle=SWC_LS[i_swc], alpha=0.7, zorder=2)

            marker = ['o', 's', '^'][i_swc]
            ax.scatter([ci_s], [a_s], s=28, marker=marker,
                       color=SOY_COLOR, edgecolor=era_color,
                       linewidth=1.0, zorder=6)
            ax.scatter([ci_m], [a_m], s=28, marker=marker,
                       color=MAIZE_COLOR, edgecolor=era_color,
                       linewidth=1.0, zorder=6)

    ax.set_xlabel(r'Intercellular CO$_2$, $C_i$ ($\mu$mol mol$^{-1}$)', fontsize=10)
    ax.set_ylabel(r'Net photosynthesis, $A$ ($\mu$mol m$^{-2}$ s$^{-1}$)', fontsize=10)
    ax.set_xlim(0, 700)
    ax.set_ylim(0, 50)

    handles = [
        Line2D([], [], color=SOY_COLOR,   linewidth=2.2, label='C3 Soybean'),
        Line2D([], [], color=MAIZE_COLOR, linewidth=2.2, label='C4 Maize'),
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label='Soil moisture'),
    ] + [
        Line2D([], [], color='gray', linewidth=1.0,
               linestyle=SWC_LS[i], marker=['o','s','^'][i], markersize=4,
               label=f'  {lbl}')
        for i, lbl in enumerate(SWC_LABELS)
    ] + [
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label=r'Atmospheric CO$_2$'),
    ] + [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
               markersize=5, markeredgecolor='black', markeredgewidth=0.4,
               label=f'  {lbl} ({ca} ppm)')
        for c, lbl, ca in zip(CA_COLORS, CA_LABELS, CA_LEVELS)
    ]
    leg = ax.legend(handles=handles, loc='lower right', frameon=False,
                    fontsize=7.5, handlelength=1.8, handletextpad=0.4,
                    labelspacing=0.15, bbox_to_anchor=(1.0, 0.01))
    for text in leg.get_texts():
        if text.get_text() in ('Soil moisture', r'Atmospheric CO$_2$'):
            text.set_fontweight('bold')


# ===========================================================================
# Panel (b): DCFE vs SWC
# ===========================================================================

def plot_dcfe(ax, soy_mc, maize_mc):
    swc_values = np.linspace(0.175, 0.35, 60)
    (s_med, s_lo, s_hi), \
    (m_med, m_lo, m_hi), \
    (d_med, d_lo, d_hi) = compute_dcfe_vs_swc(swc_values, VPD_REF,
                                              soy_mc, maize_mc)

    ax.fill_between(swc_values, s_lo, s_hi, color=SOY_COLOR,
                    alpha=0.18, linewidth=0, zorder=1)
    ax.fill_between(swc_values, m_lo, m_hi, color=MAIZE_COLOR,
                    alpha=0.18, linewidth=0, zorder=1)
    ax.fill_between(swc_values, d_lo, d_hi, color=DIFF_COLOR,
                    alpha=0.15, linewidth=0, zorder=1)

    ax.plot(swc_values, s_med, color=SOY_COLOR,   linewidth=2.2,
            zorder=3, label=r'CFE$_{\mathrm{soy}}$')
    ax.plot(swc_values, m_med, color=MAIZE_COLOR, linewidth=2.2,
            zorder=3, label=r'CFE$_{\mathrm{corn}}$')
    ax.plot(swc_values, d_med, color=DIFF_COLOR,  linewidth=2.4,
            zorder=4, label=r'$\Delta$CFE (soy $-$ corn)')

    # Mark field capacity inside the visible range (WP=0.10 is off-plot)
    ax.axvline(SWC_FC, color='gray', linewidth=1.0, linestyle=':', alpha=0.7)
    ax.text(SWC_FC - 0.003, 0.02, f'FC\n({SWC_FC})',
            fontsize=6.5, color='gray', va='bottom', ha='right',
            transform=ax.get_xaxis_transform())

    ax.axhline(0, color='gray', linewidth=0.7, linestyle='--', alpha=0.5)
    ax.set_xlabel(r'Root-zone soil water content (m$^{3}$ m$^{-3}$)', fontsize=10)
    ax.set_ylabel(r'CO$_2$ sensitivity (% ppm$^{-1}$)', fontsize=10)
    ax.set_xlim(0.175, 0.35)
    ax.set_ylim(0, 0.45)

    ax.legend(loc='upper right', fontsize=8, frameon=False,
              handlelength=2.0, bbox_to_anchor=(0.98, 0.98))
    ax.text(0.03, 0.97, f'VPD = {VPD_REF:.1f} kPa',
            transform=ax.transAxes, fontsize=8, ha='left', va='top',
            color='gray')


# ===========================================================================
# Panel (c): gs vs Ca at two VPD x two SWC
# ===========================================================================

def plot_gs_vs_ca(ax):
    ca_sweep = np.linspace(280, 460, 80)

    for swc, alpha in zip(PANEL_C_SWC, PANEL_C_ALPHA):
        for vpd, ls in zip(PANEL_C_VPD, PANEL_C_VPD_LS):
            beta = beta_swc(swc)
            soy_str   = soy_biochem_stressed(SOY_PARAMS,    beta)
            maize_str = maize_biochem_stressed(MAIZE_PARAMS, beta)

            gs_s_arr = np.array([gs_soy(Ca, vpd, swc, params=soy_str)
                                 for Ca in ca_sweep])
            gs_m_arr = np.array([gs_maize(Ca, vpd, swc, params=maize_str)
                                 for Ca in ca_sweep])

            ax.plot(ca_sweep, gs_s_arr, color=SOY_COLOR,   linewidth=1.8,
                    linestyle=ls, alpha=alpha)
            ax.plot(ca_sweep, gs_m_arr, color=MAIZE_COLOR, linewidth=1.8,
                    linestyle=ls, alpha=alpha)

    for Ca, clr in zip(CA_LEVELS, CA_COLORS):
        ax.axvline(Ca, color=clr, linewidth=0.6, linestyle=':',
                   alpha=0.55, zorder=1)

    ax.set_xlabel(r'Atmospheric CO$_2$, $C_a$ ($\mu$mol mol$^{-1}$)', fontsize=10)
    ax.set_ylabel('Stomatal conductance, $g_s$\n'
                  r'(mol H$_2$O m$^{-2}$ s$^{-1}$)', fontsize=10)
    ax.set_xlim(280, 460)
    ax.set_ylim(bottom=0)

    handles = [
        Line2D([], [], color=SOY_COLOR,   linewidth=2.0, label='C3 Soybean'),
        Line2D([], [], color=MAIZE_COLOR, linewidth=2.0, label='C4 Maize'),
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label='VPD'),
        Line2D([], [], color='gray', linewidth=1.5, linestyle='-',
               label=f'  {PANEL_C_VPD[0]:.1f} kPa'),
        Line2D([], [], color='gray', linewidth=1.5, linestyle='--',
               label=f'  {PANEL_C_VPD[1]:.1f} kPa'),
        Line2D([], [], color='none', label=''),
        Line2D([], [], color='none', label='Soil moisture'),
        Line2D([], [], color='gray', linewidth=1.5, alpha=PANEL_C_ALPHA[0],
               label=f'  SWC = {PANEL_C_SWC[0]:.2f} (wet)'),
        Line2D([], [], color='gray', linewidth=1.5, alpha=PANEL_C_ALPHA[1],
               label=f'  SWC = {PANEL_C_SWC[1]:.2f} (dry)'),
    ]
    leg = ax.legend(handles=handles, loc='upper right', frameon=False,
                    fontsize=7.5, handlelength=1.8, labelspacing=0.15,
                    bbox_to_anchor=(0.98, 0.98))
    for text in leg.get_texts():
        if text.get_text() in ('VPD', 'Soil moisture'):
            text.set_fontweight('bold')


# ===========================================================================
# Main
# ===========================================================================

def main():
    soy_mc   = sample_soy_params(N_MC)
    maize_mc = sample_maize_params(N_MC)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    plot_aci(axes[0],       soy_mc, maize_mc)
    plot_dcfe(axes[1],      soy_mc, maize_mc)
    plot_gs_vs_ca(axes[2])

    for i, letter in enumerate(['(a)', '(b)', '(c)']):
        axes[i].text(-0.10, 1.02, letter, transform=axes[i].transAxes,
                     fontsize=10, fontweight='bold', va='bottom')

    fig.tight_layout(w_pad=1.5)
    out_path = os.path.join(_HERE, 'figure4_combined.png')
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {out_path}")

    print_diagnostics(vpd=VPD_REF)


if __name__ == "__main__":
    main()
