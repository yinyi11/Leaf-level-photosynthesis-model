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

Ainsworth, E.A. & Long, S.P. (2005). What have we learned from 15 years of free-air CO₂ enrichment (FACE)? A meta-analytic review of the responses of photosynthesis, canopy properties and plant production to rising CO₂. *New Phytologist*, 165, 351–372.

Ainsworth, E.A., Sanz-Saez, A. & Leisner, C.P. (2025). Crop photosynthesis in a high-CO₂ world: insights from two decades of free-air CO₂ enrichment. *Philosophical Transactions of the Royal Society B*, 380, 20240230.

Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998). *Crop Evapotranspiration: Guidelines for Computing Crop Water Requirements*. FAO Irrigation and Drainage Paper 56. FAO, Rome.

Bellasio, C. (2024). Drought responses of C4 photosynthesis and the CO₂-concentrating mechanism. *New Phytologist*. https://doi.org/10.1111/nph.19536

Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis, A.R. & Long, S.P. (2001). Improved temperature response functions for models of Rubisco-limited photosynthesis. *Plant, Cell & Environment*, 24, 253–259.

Bernacchi, C.J., Leakey, A.D.B., Heady, L.E., Morgan, P.B., Dohleman, F.G., McGrath, J.M., Gillespie, K.M., Wittig, V.E., Rogers, A., Long, S.P. & Ort, D.R. (2005). Hourly and seasonal variation in carbon assimilation, stomatal conductance and transpiration of soybean grown under free-air CO₂ enrichment. *Planta*, 220, 434–446.

Bernacchi, C.J., Morgan, P.B., Ort, D.R. & Long, S.P. (2007). The growth of soybean under free air CO₂ enrichment (FACE) stimulates photosynthesis while decreasing in vivo Rubisco capacity. *Plant Physiology*, 143, 134–144.

Bishop, K.A., Betzelberger, A.M., Long, S.P. & Ainsworth, E.A. (2015). Is there potential to adapt soybean (*Glycine max* Merr.) to future [CO₂]? An analysis of the yield response of 18 genotypes in free-air CO₂ enrichment. *Plant, Cell & Environment*, 38, 1765–1774.

De Kauwe, M.G., Medlyn, B.E., Zaehle, S., Walker, A.P., Dietze, M.C., Wang, Y.-P., Luo, Y., Jain, A.K., El-Masri, B., Hickler, T., Wårlind, D., Weng, E., Parton, W.J., Thornton, P.E., Wang, S., Prentice, I.C., Asao, S., Smith, B., McCarthy, H.R., Iversen, C.M., Hanson, P.J., Warren, J.M., Oren, R. & Norby, R.J. (2015). Where does the carbon go? A model–data intercomparison of vegetation carbon allocation and turnover processes at two temperate forest free-air CO₂ enrichment sites. *Geoscientific Model Development*, 8, 431–452.

Egea, G., Verhoef, A. & Vidale, P.L. (2011). Towards an improved and more flexible representation of water stress in coupled photosynthesis–stomatal conductance models. *Agricultural and Forest Meteorology*, 151, 1369–1384.

Farquhar, G.D., von Caemmerer, S. & Berry, J.A. (1980). A biochemical model of photosynthetic CO₂ assimilation in leaves of C3 species. *Planta*, 149, 78–90.

Flexas, J. & Medrano, H. (2002). Drought-inhibition of photosynthesis in C3 plants: stomatal and non-stomatal limitations revisited. *Annals of Botany*, 89, 183–189.

Ghannoum, O. (2009). C4 photosynthesis and water stress. *Annals of Botany*, 103, 635–644.

Leakey, A.D.B., Bernacchi, C.J., Dohleman, F.G., Ort, D.R. & Long, S.P. (2006). Will photosynthesis of maize (*Zea mays*) in the US Corn Belt increase in future [CO₂] rich atmospheres? An analysis of diurnal courses of CO₂ uptake under free-air concentration enrichment (FACE). *Plant Physiology*, 140, 779–790.

Lin, Y.-S., Medlyn, B.E., Duursma, R.A., Prentice, I.C., Wang, H., Baig, S., Eamus, D., de Dios, V.R., Mitchell, P., Ellsworth, D.S., de Beeck, M.O., Wallin, G., Uddling, J., Tarvainen, L., Linderson, M.-L., Cernusak, L.A., Nippert, J.B., Ocheltree, T.W., Tissue, D.T., Martin-StPaul, N.K., Rogers, A., Warren, J.M., De Angelis, P., Hikosaka, K., Han, Q., Onoda, Y., Gimeno, T.E., Barton, C.V.M., Bennie, J., Bonal, D., Bosc, A., Löw, M., Macinins-Ng, C., Rey, A., Rowland, L., Setterfield, S.A., Tausz-Posch, S., Zaragoza-Castells, J., Broadmeadow, M.S.J., Drake, J.E., Freeman, M., Ghannoum, O., Hutley, L.B., Kelly, J.W., Kikuzawa, K., Kolari, P., Koyama, K., Limousin, J.-M., Meir, P., Lola da Costa, A.C., Mikkelsen, T.N., Salinas, N., Sun, W. & Wingate, L. (2015). Optimal stomatal behaviour around the world. *Nature Climate Change*, 5, 459–464.

Long, S.P., Ainsworth, E.A., Leakey, A.D.B., Nösberger, J. & Ort, D.R. (2006). Food for thought: lower-than-expected crop yield stimulation with rising CO₂ concentrations. *Science*, 312, 1918–1921.

Massad, R.-S., Tuzet, A. & Bethenod, O. (2007). The effect of temperature on C4-type leaf photosynthesis parameters. *Plant, Cell & Environment*, 30, 1191–1204.

Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Prentice, I.C., Barton, C.V.M., Crous, K.Y., De Angelis, P., Freeman, M. & Wingate, L. (2011). Reconciling the optimal and empirical approaches to modelling stomatal conductance. *Global Change Biology*, 17, 2134–2144.

Ripley, B.S., Gilbert, M.E., Ibrahim, D.G. & Osborne, C.P. (2007). Drought constraints on C4 photosynthesis: stomatal and metabolic limitations in C3 and C4 subspecies of *Alloteropsis semialata*. *Journal of Experimental Botany*, 58, 1351–1363.

Rogers, A., Allen, D.J., Davey, P.A., Morgan, P.B., Bernacchi, C.J., Mahoney, J., Gielen, B., Fiscus, E.L., Stitt, M., Hendrey, G.R., Athanasiou, A., Frenck, G., Zhu, X.-G., DeLucia, E.H., Ort, D.R. & Long, S.P. (2004). Leaf photosynthesis and carbohydrate dynamics of soybeans grown throughout their life-cycle under free-air carbon dioxide enrichment. *Plant, Cell & Environment*, 27, 449–458.

Saxton, K.E. & Rawls, W.J. (2006). Soil water characteristic estimates by texture and organic matter for hydrologic solutions. *Soil Science Society of America Journal*, 70, 1569–1578.

Verhoef, A. & Egea, G. (2014). Modeling plant transpiration under limited soil water: comparison of different plant and soil hydraulic parameterizations and preliminary implications for their use in land surface models. *Agricultural and Forest Meteorology*, 191, 22–32.

von Caemmerer, S., Evans, J.R., Hudson, G.S. & Andrews, T.J. (1994). The kinetics of ribulose-1,5-bisphosphate carboxylase/oxygenase in vivo inferred from measurements of photosynthesis in leaves of transgenic tobacco. *Planta*, 195, 88–97.

von Caemmerer, S. (2000). *Biochemical Models of Leaf Photosynthesis*. CSIRO Publishing, Collingwood, Australia.

von Caemmerer, S. (2013). Steady-state models of photosynthesis. *Plant, Cell & Environment*, 36, 1617–1630.

Zhou, S., Duursma, R.A., Medlyn, B.E., Kelly, J.W.G. & Prentice, I.C. (2014). How should we model plant responses to drought? An analysis of stomatal and non-stomatal responses to water stress. *Agricultural and Forest Meteorology*, 182–183, 204–214.
