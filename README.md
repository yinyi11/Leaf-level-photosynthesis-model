# Differential CO₂ Fertilization Effect: A Leaf-Level Photosynthesis Model

A theoretical model diagnosing why soybean (C3) and maize (C4) respond differently to rising atmospheric CO₂, and why that difference collapses under drought. The model couples biochemical photosynthesis (Farquhar–von Caemmerer–Berry for C3; von Caemmerer 2000 for C4) to an anchored Medlyn stomatal sub-model with soil-water stress.

## Repository structure

```
photosynthesis_models.py       # All model code and parameters (import this)
generate_figure4_combined.py   # Figure 4: A-Ci curves, ΔCFE vs SWC, gs vs Ca
generate_figure_face.py        # FACE figure: gs vs Ca, A vs SWC, ΔA% vs SWC
```

## Quick start

```bash
pip install numpy scipy matplotlib
python generate_figure4_combined.py   # → figure4_combined.png
python generate_figure_face.py        # → figure_face.png
python photosynthesis_models.py       # run unit tests
```

---

## 1. What the model does

Three independent physiological controls modulate leaf-level net photosynthesis (*A*) in response to atmospheric CO₂ (*C*a), VPD, and root-zone soil water content (SWC):

| Mechanism | Species | Drought response |
|---|---|---|
| **(I) Stomatal aperture** (*g*s) | Both | Anchored Medlyn USO × β(SWC) |
| **(II) C3 biochemical capacity** | Soybean | β-ramped *V*cmax & *J*max (Flexas framework) |
| **(III) C4 CCM capacity** (*V*pmax) | Maize | β-scaled PEPC activity |

The differential CO₂ fertilization effect (CFE) between species, and its collapse under drought, emerge from the interaction of these three mechanisms with species-specific A–*C*i curve geometry. No species-asymmetric stomatal rule is required.

---

## 2. Stomatal conductance

### 2.1 Functional form — anchored Medlyn USO

For both species:

$$g_s(C_a, \text{VPD}, \text{SWC}) = \beta(\text{SWC}) \times g_{s,\text{fc,obs}} \times \kappa(C_a, \text{VPD}, g_1)$$

where the dimensionless Ca-kernel is:

$$\kappa(C_a) = \frac{g_{s,\text{Medlyn}}(C_a, \text{VPD}, g_1)}{g_{s,\text{Medlyn}}(C_{a,\text{ref}}, \text{VPD}_\text{ref}, g_1)}, \quad C_{a,\text{ref}} = 370\ \mu\text{mol mol}^{-1}$$

$$g_{s,\text{Medlyn}} = 1.6\left(1 + \frac{g_1}{\sqrt{\text{VPD}}}\right)\frac{A}{C_a}$$

The kernel captures the Ca-dependence of optimal stomatal behaviour. The Medlyn kernel is evaluated at the leaf's current biochemical state (stressed or unstressed), making it self-consistent with the demand curve used in the operating-point solve.

### 2.2 Why anchor magnitude empirically

The Medlyn USO formulation is theoretically robust on the *shape* of the *g*s(*C*a, VPD) response, but absolute *g*s inherits errors in *A* from the no-mesophyll-conductance photosynthesis model (~15–30% underestimate). Anchoring *g*s,fc to direct field observations sidesteps this, while the kernel preserves the theoretical Ca-dependence.

### 2.3 Empirical anchors (well-watered, *C*a,ref = 370 ppm, VPD = 1.5 kPa)

| Species | *g*s,fc,obs | Source |
|---|---|---|
| C3 Soybean | 0.55 mol H₂O m⁻² s⁻¹ | Ainsworth & Long (2005) SoyFACE meta-analysis |
| C4 Maize | 0.30 mol H₂O m⁻² s⁻¹ | Leakey et al. (2006b) Table II, DOY 190/201/215 mean |

The soy anchor yields *C*i/*C*a ≈ 0.77 in this no-*g*m model vs. observed ~0.73 — a small offset attributable to omitted mesophyll conductance.

### 2.4 Medlyn slope *g*₁

*g*₁ is taken from Lin et al. (2015) PFT-level fits to a global database (314 species, 56 field studies):

| Species | *g*₁ (kPa⁰·⁵) | PFT |
|---|---|---|
| Soybean | 5.79 ± 0.64 | C3 cropland |
| Maize | 1.62 ± 0.13 | C4 grassland (no separate C4 cropland PFT in Lin's database) |

These are literature-default values chosen prior to any comparison with FACE observations — not values tuned to match outcomes.

### 2.5 Why the same Medlyn rule produces different Ca-responses

**C3 soybean** (*g*₁ = 5.79, *C*i/*C*a ≈ 0.77): operates on the *rising* portion of its A–*C*i curve. As *C*a rises, both *C*i and *A* rise nearly proportionally, so *g*s,Medlyn ∝ *A*/*C*a is only weakly *C*a-dependent. Kernel change: −4.5% over 336→420 ppm; −20% over 376→550 ppm.

**C4 maize** (*g*₁ = 1.62, *C*i/*C*a ≈ 0.37): operates on the *saturated plateau* of its A–*C*i curve. *A* is nearly *C*a-independent, so *g*s,Medlyn ∝ 1/*C*a. Kernel change: −20% over 336→420 ppm; −32% over 376→550 ppm.

### 2.6 Validation against FACE observations

Comparison at ambient (376 ppm) → elevated (550 ppm), well-watered, VPD = 1.5 kPa:

| Variable | Soybean model | Soybean obs | Maize model | Maize obs |
|---|---|---|---|---|
| Δ*g*s | −20% | −16% (Bernacchi 2007), −22% (A&L 2005 meta) | −32% | −34% (Leakey 2006b Tab. II) |
| Δ*C*i/*C*a | preserved | preserved | preserved | preserved (*p* = 0.52) |
| Δ*A* | +22% | +15–25% (SoyFACE meta) | +1% | ~0% (n.s.) |

All predicted responses fall within observed ranges for both species.

---

## 3. Soil-water stress — Verhoef-Egea linear β

$$\beta(\text{SWC}) = \text{clip}\!\left(\frac{\text{SWC} - \text{SWC}_\text{wp}}{\text{SWC}_\text{fc} - \text{SWC}_\text{wp}},\ 0,\ 1\right)$$

Standard linear β from Verhoef & Egea (2014) and Egea et al. (2011); the same form used in CLM, JULES, and ORCHIDEE.

**Anchors** (Corn Belt silt loam / Mollisol, Saxton & Rawls 2006; Allen et al. 1998 FAO-56):

| Parameter | Value | Definition |
|---|---|---|
| SWC_wp | 0.10 m³ m⁻³ | Permanent wilting point |
| SWC_fc | 0.31 m³ m⁻³ | Field capacity |

β multiplies **both** *g*s (stomatal closure) and biochemical capacity (§4). Using the same β for both pathways captures the leading-order leaf-hydraulic / metabolic coupling without introducing a second free parameter.

---

## 4. Biochemical drought response

### 4.1 C3 soybean — *V*cmax and *J*max

$$V_{c\max,\text{eff}} = V_{c\max} \times [f_\min + (1 - f_\min)\,\beta], \quad f_\min = 0.60$$

$$J_{\max,\text{eff}} = J_{\max} \times [f_\min + (1 - f_\min)\,\beta]$$

Following Flexas & Medrano (2002): stomatal limitation dominates mild–moderate stress; metabolic limitations (*V*cmax, *J*max decline) emerge only at severe stress. The floor *f*min = 0.60 is consistent with observations that even severe drought reduces *V*cmax/*J*max to ~40–60% of well-watered values rather than to zero (Zhou et al. 2014). The central conclusion is robust to *f*min ∈ [0.40, 0.70].

### 4.2 C4 maize — *V*pmax (CCM capacity)

$$V_{p\max,\text{eff}} = V_{p\max} \times \beta$$

PEPC activity declines under water stress in C4 grasses (Ripley et al. 2007; Ghannoum 2009; Bellasio 2024). Direct β scaling (no *f*min) so β = 0 fully shuts down the CCM. *V*cmax (Rubisco) and *J*max are **not** scaled — Rubisco resides in bundle-sheath cells, anatomically protected from cytosolic dehydration.

---

## 5. Monte Carlo parameter uncertainty

Cultivar variation is represented by sampling from biologically realistic ranges:

**Soybean:**

| Parameter | Distribution | Units |
|---|---|---|
| *V*cmax | U[105, 145] | µmol m⁻² s⁻¹ |
| *J*/*V* | U[1.45, 1.75] | — |
| *R*d | U[1.0, 2.0] | µmol m⁻² s⁻¹ |
| *J*max | = (*J*/*V*) × *V*cmax | µmol m⁻² s⁻¹ |

*J*max is parameterised through the *J*/*V* ratio (not sampled independently) to avoid unphysical combinations (e.g. *J*/*V* = 1.2) that would distort the CFE distribution. The *J*/*V* range covers documented soybean variation: Bernacchi et al. (2005) ≈ 1.66; Rogers et al. (2004) 1.6–1.8; Bishop et al. (2015) 1.5–1.7.

**Maize:**

| Parameter | Distribution | Units |
|---|---|---|
| *V*cmax | U[40, 55] | µmol m⁻² s⁻¹ |
| *V*pmax | U[95, 135] | µmol m⁻² s⁻¹ |
| *K*p | U[60, 100] | µmol mol⁻¹ |
| *J*max | U[350, 450] | µmol m⁻² s⁻¹ |

---

## 6. Central mechanism: why ΔCFE collapses under drought

**Well-watered (β = 1):** Maize *C*i ≈ 137 µmol mol⁻¹, above the C4 saturation knee (~75–90 in this model). *A* is on the plateau → raising *C*a leaves *A* unchanged → CFE_maize ≈ 0. Soybean operates on the rising portion at *C*i ≈ 286 → CFE_soy ≈ +0.27 %/ppm. ΔCFE ≈ +0.25 %/ppm.

**Under drought (β < 1):** Two mechanisms push maize into the CO₂-sensitive regime:
1. *g*s ↓ → *C*i ↓ via supply-demand balance
2. *V*pmax ↓ → A–*C*i curve shifts down and its knee shifts left

Together these move maize *C*i from above the knee to at or below it. Once on the curved portion, raising *C*a produces a real *A* response → CFE_maize > 0. Soybean *A* decreases in magnitude but remains on the rising portion; CFE_soy stays roughly constant. **ΔCFE collapses because CFE_maize rises while CFE_soy is stable.**

---

## 7. References

| Citation | Journal | Details |
|---|---|---|
| Ainsworth & Long (2005) | *New Phytol.* | 165: 351–372 |
| Ainsworth & Rogers (2007) | *Plant Cell Environ.* | 30: 258–270 |
| Ainsworth, Sanz-Saez & Leisner (2025) | *Phil. Trans. R. Soc. B* | 380: 20240230 |
| Allen et al. (1998) | FAO Irrigation Drainage Paper 56 | |
| Bellasio (2024) | *New Phytol.* | CCM drought response |
| Bernacchi et al. (2001) | *Plant Cell Environ.* | 24: 253–259 |
| Bernacchi et al. (2005) | *Planta* | 220: 434–446 |
| Bernacchi et al. (2007) | *Plant Physiol.* | 143: 134–144 |
| Bishop et al. (2015) | *Plant Cell Environ.* | 38: 1765–1774 |
| De Kauwe et al. (2015) | *Geosci. Model Dev.* | 8: 431–452 |
| Egea et al. (2011) | *Agric. For. Meteorol.* | 151: 1369–1384 |
| Farquhar, von Caemmerer & Berry (1980) | *Planta* | 149: 78–90 |
| Flexas & Medrano (2002) | *Ann. Bot.* | 89: 183–189 |
| Ghannoum (2009) | *Ann. Bot.* | 103: 635–644 |
| Leakey et al. (2006a) | *Plant Cell Environ.* | 29: 1794–1800 |
| Leakey et al. (2006b) | *Plant Physiol.* | 140: 779–790 |
| Leakey et al. (2009) | *J. Exp. Bot.* | 60: 2859–2876 |
| Lin et al. (2015) | *Nat. Clim. Change* | 5: 459–464 |
| Long et al. (2006) | *Science* | 312: 1918–1921 |
| Medlyn et al. (2011) | *Glob. Change Biol.* | 17: 2134–2144 |
| Ripley et al. (2007) | *J. Exp. Bot.* | 58: 1351–1363 |
| Rogers et al. (2004) | *Plant Cell Environ.* | 27: 449–458 |
| Saxton & Rawls (2006) | *Soil Sci. Soc. Am. J.* | 70: 1569–1578 |
| Verhoef & Egea (2014) | *Agric. For. Meteorol.* | 191: 22–32 |
| von Caemmerer (2000) | *Biochemical Models of Leaf Photosynthesis* | CSIRO Publishing |
| Zhou et al. (2014) | *Tree Physiol.* | 34: 1035–1046 |
