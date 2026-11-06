# Scientific model, validation and references


---

# Evolving giant-star surface: statistical model and limits

This statistical model replaces the evolved-star solar-atlas extrapolation. It does not modify the accepted MIST/MESA track, C1 evolutionary playback coordinate, pulsation, projected radius, camera, gas trajectories, ionization, opacity, intrinsic stellar luminosity or v4 nebular imaging. The selected model passed native diagnostics and final publication review; see VALIDATION.md.

## Why a different representation is necessary

The former adapter sampled the same solar radiation maps at coordinates proportional to an instantaneous pressure-height magnification. Its cubic repair moderated that magnification, but still enlarged and contracted existing features. A continuous magnification curve cannot provide independent cell turnover, deformation or fragmentation. The new model has no giant-stage image atlas or magnification variable.

## Literature constraints and their domain

[Tremblay et al. (2013), A&A 557, A7](https://arxiv.org/abs/1307.2810) provide approximate characteristic size and Eulerian intensity-decay fits:

\[
d_{\rm fit}[\mathrm{km}]=13.5g^{-1}(T_{\rm eff}-300\log g)^{1.75}10^{0.05[\mathrm{Fe/H}]}
\]
\[
\tau_{\rm fit}[\mathrm{s}]=2.08g^{-1}(T_{\rm eff}-300\log g)^{1.75}10^{0.05[\mathrm{Fe/H}]}
\]

Here g is in cgs and the solar-metallicity track uses [Fe/H]=0. We mark the relevant non-degenerate, cool-star domain as 3600–6750 K and log g 1–4.5, subject to the paper's convective-domain condition. These fits have appreciable scatter and do not determine a full surface flow. The white-dwarf portion of that paper's grid does not validate their use on hot, low-gravity post-AGB stars.

At the coolest, lowest-gravity RGB/AGB epochs, these conditions fail. The code then labels a **pressure-height/sound-speed continuation**, anchored at 3600 K / log g 1: d∝T/g and τ∝√T/g. This is a physically motivated extrapolation, not a validated extension of the empirical fit. Every epoch records the method and the in-grid flag.

[Ludwig & Kučinskas (2012), A&A 547, A118](https://doi.org/10.1051/0004-6361/201220264) find strong giant granulation, approximately 18% intensity contrast and a typical 5 Gm diameter in a 3660 K / log g 1 model. Their size statistic and individual stellar model are not identical to the CIFIST spectral characteristic. We use these results to constrain plausible contrast and scale ranges, not to claim an exact match to that atmosphere.

[Freytag, Liljegren & Höfner (2017), A&A 600, A137](https://arxiv.org/abs/1702.05433) distinguish long-lived giant structures, shorter-lived fine surface structure, pulsation and shocks. Their extended-atmosphere intensity need not directly trace convection at a sharply defined photosphere. Our weak large-cell plus faster fine-cell populations are statistical representations of that hierarchy; shocks, molecular transfer, dynamical atmospheric extension and self-consistent convection are not calculated.

## Persistent spherical cells

`src/giant_convection.py` creates three populations with reproducible seeds 6006, 106009 and 206012. Birth times are continuous, spatial directions are isotropic, and individual sizes/lifetimes have a finite distribution. Random numbers are drawn once when creating the population; evaluation never regenerates per-frame noise.

- Main granules use the literature-informed newborn diameter d/R
- A weaker population has larger structures and approximately 3.3 times longer representative turnover
- A smaller population has approximately half the main scale and 0.6 times its turnover, providing evolving internal detail

At full-star display scale, tiny physical cells cannot be resolved. Newborn main-cell diameter has an explicit 0.040 radian resolution floor and 0.28 radian ceiling. Fine cells have a 0.024 radian floor; larger structures are bounded at 0.65 radian. Size scatter is applied afterward. These are spatial visualization choices and are recorded separately from physical reference lengths. At 4 seconds the main-cell floor is about 10.5 times the instantaneous fitted angular scale; it becomes unnecessary near the RGB maximum. Earlier cells can retain smaller angular dimensions than the current instantaneous target because their identity is not rescaled with the star.

**A living cell never changes size because R(t) changes.** Its stellar-coordinate direction follows a smooth individual tangent trajectory. The reference drift is 0.23 d_fit/τ_fit, of order 1–2 km/s for the selected parameters; its displayed motion shares the cell's spatial/time enhancements. This is prescribed statistical transport, not a solved velocity field.

Sites have independently phased elliptical shape deformation and a smooth sine-squared life envelope. A newly strengthening site subdivides existing territory; as a site fades, neighboring cells fill its area. This produces continuous statistical splitting/replacement and boundary motion without forcing an entire surface to expand or reverse. Site IDs and trajectories remain available for measurement.

At a surface direction, the two most competitive local site scores define a continuous cellular interior/boundary field. Warm interiors have a radial brightness gradient; close site-score ties form cooler lanes. The three populations create nested but non-repeating spatial detail. No cube faces, fixed chart windows, low-order harmonic pattern, image zoom, turbulence-noise layer or static texture tile is used.

## Representative atmospheric clock

Literal tracking over the millions of evolutionary years compressed into seconds would erase every individual convective feature. Each newborn site instead receives a representative display turnover

\[
\tau_d=0.70+0.20\log_{10}(1+\tau_{\rm fit}/1000\,\mathrm{s})\quad\mathrm{playback\ seconds}
\]

and a lifetime drawn around 2.6 τ_d. The populations' clock ratios differ as described above. This deliberately compressed clock is independent of camera, pulsation and the secular stellar-age mapping. It is not a measured relaxation time for the evolving MIST star. Published Eulerian autocorrelation times, material-site lifetimes and the measured movie's intensity autocorrelation are different quantities and are reported separately.

Birth rates vary with the expected area per cell and lifetime. Existing cells retain their identity; changed stellar parameters influence newborn statistics rather than remapping all old coordinates. The model can therefore respond gradually without the former coherent map breathing.

## Radiation and luminosity

The positive field f is normalized over the complete sphere, not the current camera hemisphere:

\[
\langle f\rangle_{4\pi}=1,\qquad T_{\rm local}^4=T_{\rm eff}^4 f.
\]

Thus the statistical thermal fluctuations preserve the prescribed global bolometric luminosity. The adopted radial pulsation still changes R and T at fixed L. Surface emission continues through the existing spectral/CIE infrastructure and normalized limb response; camera exposure, color treatment, detector response and optical PSFs are unchanged. The model does not apply a new false-color temperature map. This surface-radiation closure is not fed back into the accepted wind or ionizing-photon calculation: their source remains the secular stellar track, and the hot/unresolved source uses the inherited compact-star branch.

Each population's cellular heat field is standardized to full-sphere zero mean
and unit variance, giving z_main, z_large and z_fine. The adopted contrast is
explicitly prescribed:

\[
f_0=\exp\{A z_{\rm main}+0.045z_{\rm large}+0.105z_{\rm fine}\},
\quad A=0.18+0.020S[(2.5-\log g)/1.5],\quad f=f_0/\langle f_0\rangle.
\]

S is the clipped cubic smoothstep. These coefficients specify statistical
bolometric contrast; they are not a calculated convective heat flux or a fit
to the selected MIST star. The stellar-coordinate surface field contains no
independent false-color or sharpening operation. Its visible contrast follows
the existing temperature-dependent spectrum and detector/optical pipeline.

The spherical field is sampled on a fixed periodic longitude/colatitude grid. Longitude interpolation wraps continuously. Each pole is one physical point, with one value shared by all longitudes. Exact integration weights for the piecewise-linear colatitude basis make full-sphere normalization independent of viewing direction. This sampling grid is not an image being stretched with R.

The main-sequence opening retains No. 004. A smooth linear-radiance source handoff occurs as R increases from 1.6 to 3 solar radii. This is a disclosed transition between surface models, not an encoded image dissolve. Beyond that handoff no solar texture is used for the giant. The hot/unresolved remnant continues to use the accepted No. 003 path.



---

# No.006: literature review and implementation decision

Review date: 4 November 2026. Primary papers were downloaded and their relevant
methods/results read; URLs, PDF hashes and page counts are in
`research/papers/manifest.json`. Those local research copies are not a proposed
redistribution licence. This report precedes the v3 numerical experiment.

## Scenario and its important limit

Keep the **same nonrotating 1.000 M☉, Z=0.0142857 MIST/MESA sequence** used in v2,
including its native mass loss and surface abundances. Do not splice a more
massive, brighter central star into it. Its final mass is 0.539832 M☉. From our
post-AGB reference (Teff=7152 K), the existing interpolated track first reaches
10/20/30/50/100 kK at approximately 95/506/1535/5002/14623 yr. These are measured
from **our** track, not numbers copied from another grid.

The opening retains v2's explicit No004 reference match: at4.57Gyr, native
MIST luminosity/temperature are about1.106L☉/5848K. Small offsets produce
1L☉/5772K at the opening, with radius calculated from Stefan–Boltzmann. Those
offsets smoothly vanish over the next1Gyr; mass/composition and all late
giant/post-AGB values remain native. This is not a new solar-calibrated MESA
calculation. `data/evolution/clock-v3.csv` records the selected film's secular
time mapping and unpulsed stellar properties; its atmospheric clock is separate.

[Miller Bertolami 2016](https://doi.org/10.1051/0004-6361/201526577), Table 3
and §4, gives strongly mass-dependent H-burning crossing times. Its 1 M☉
Z=0.02 example ends at 0.5281 M☉, with transition/crossing times of
9.14/24.9 kyr; the Z=0.01 example ends at 0.5319 M☉ and is appreciably slower
(36.1/63.0 kyr). They are useful bounds on sensitivity, **not interchangeable
with our different MIST structure**. Updated tracks evolve faster than older
grids, but this does not establish a unique fate for the actual Sun.

[Gesicki, Zijlstra & Miller Bertolami 2018](https://doi.org/10.1038/s41550-018-0453-9),
Figs 1–3 and the final two paragraphs before Methods, specifically place the
Sun near the low-mass PN limit. Their 1 M☉ case starts H ionization near 5 kyr
and [O III] near 12.5 kyr; expansion makes it faint. A helium-burning departure
can fail to produce a visible PN. Their 0.10–0.15 M☉ idealized shells and
maximum-nebula luminosities must **not** replace our conserved wind. We retain
a faint, H-burning, model-dependent scenario, and disclose increased detector
sensitivity and a fixed line-color composite rather than increase intrinsic
luminosity.

## Literature-to-feature decisions

| Feature / mechanism | Verified supporting location | Applicability and implementation |
|---|---|---|
| Pulsation-assisted dusty AGB loss | [Höfner & Olofsson 2018](https://doi.org/10.1007/s00159-017-0106-5), §1.3, Eqs 1–4, Fig.2 and §2 wind measurements | Dust opacity can exceed gravity after shock levitation. Observed winds commonly 5–15 km/s, broader 3–30 km/s; condensation roughly 1000–1500 K. Retain 12 km/s, the native MIST mass-loss integral, radiative-equilibrium dust and explicit opacity assumption. We do not solve grain nucleation or infer mass loss from luminosity alone |
| Irregular atmospheric mass loading | [Höfner & Freytag 2019](https://arxiv.org/abs/1902.04074), Table 1, Figs 2–4; [Freytag & Höfner 2023](https://arxiv.org/abs/2301.11836), wind models and conclusions | Non-spherical shocks/dust clouds are supported in M-type AGB atmospheres. The 2019 1 M☉/7049 L☉ example has radius478 R☉ and Teff2417 K; it **does not calculate wind acceleration**. The 2023 models do include dust driving. Our lighter late star is not that exact model: launch inhomogeneities are prescribed, statistically correlated source conditions, not imported simulation outputs or a prediction of surviving knots |
| AGB pulsation / atmospheric extension | [Freytag, Liljegren & Höfner 2017](https://arxiv.org/abs/1702.05433), model grid, period-radius discussion and Figs 1–6; 2019 Table1 period507d at372 R☉/1 M☉ | A slowly sampled representative radial pulsation can be scaled by the dynamical-time relation P∝R^(3/2)M^(-1/2). Any small amplitude and separate readable atmosphere clock remain prescribed. No claim of computed eigenmodes or new bolometric energy; secular R/L/M stay on MIST |
| Inner wind rim, pressure-driven outer structure, halo | [Villaver, Manchado & García-Segura 2002](https://arxiv.org/abs/astro-ph/0208323), §2, §3.1, Figs 5,11,17–18 | Their low-mass case shows early ionization-driven expansion and multiple shells; halos can extend to2.3pc and be >2000 times fainter than the bright shell. Their old 1 M☉ progenitor leaves a **0.569 M☉** core, unlike ours. Adopt the mechanisms and compare velocity orders, not their exact ages/radii/luminosities. Replace the single Gaussian shell with conservative radial material cells responding to thermal and bubble pressure |
| Corrugation, compression and nonuniform fronts | [Toalá & Arthur 2014](https://doi.org/10.1093/mnras/stu1360), Eqs1–4, §4.1, Fig.5 and AppendixA for their1.0–0.569 model | Their 2D radiation hydro demonstrates thin-shell, Rayleigh–Taylor and shadowing instabilities; a1.5–0.597 example is not the Sun. We calculate radial compression and differential sector response to prescribed AGB loading, plus direct-ray shielding. **No lateral hydro, nonlinear RT growth, mixing-layer filaments or conduction is resolved**. Do not paint those structures on to imitate Fig.5 |
| Dense knots / cometary tails | [O'Dell, Henney & Ferland 2005](https://doi.org/10.1086/430803), abstract, §1 and knot radiation/flow models; [García-Segura et al.2006](https://arxiv.org/abs/astro-ph/0606205), §3/Figs1–3 | Helix knots have molecular cores near10^6cm⁻³ and photoevaporative ionized faces. The 2006 experiment abruptly stops a selected fast wind at1kyr and cannot fix knot numbers/sizes at its resolution. **Do not adopt that shutdown, add Helix knots, or claim calculated photoevaporation**. Direct-ray low-ionization shadows may appear only where calculated columns shield photons; unresolved diffuse radiation is a caveat |
| Bipolar lobes, jets, strong global shaping | [Jones & Boffin2017](https://doi.org/10.1038/s41550-017-0117), morphology/binarity review | No companion is in this scenario. Exclude jets, bipolar lobes, a binary torus and strong imposed global magnetic shaping. Modest single-star asymmetry is an initial-condition assumption, not measured solar asymmetry |
| Hα, Hβ and collisionally excited oxygen/nitrogen | [Luridiana, Morisset & Shaw2015](https://doi.org/10.1051/0004-6361/201323152), §2.1–2.3, Eqs1–3, §6 limitations | Use actual PyNeb recombination tables and n-level statistical equilibrium, with recorded atomic files, temperature/density lookup and an interpolation check. j(Hα)=ne n(H+)ε_Hα; j(5007)=ne n(O²⁺)ε_5007, etc. PyNeb supplies neither ionic populations nor a self-consistent temperature/ionizing spectrum |
| Ionizing spectra, heating and multi-element transfer | [Ferland et al.2017](https://arxiv.org/abs/1705.10877), code/atomic database and small-model accuracy discussion | Cloudy shows what a complete thermal/ionization calculation entails. **We are not running Cloudy**. Use frequency groups, H/He absorption, threshold-resolved O/N equilibrium approximations and explicit heating ceilings. No claim of Cloudy accuracy or NLTE stellar spectrum |
| Photoionization cross sections | [Verner et al.1996](https://doi.org/10.1086/177435), Eq1 and Table1; [author's numerical table](https://www.pa.uky.edu/~verner/photo.html) | Use the published analytic ground-state fits rather than arbitrary constant O/N cross sections. Recombination, charge exchange, unresolved spectral structure and thermal balance remain separately documented approximations |

## Efficient numerical plan

1. Keep all v2 evidence and source adapters. Archive the selected v2 manifest.
2. Prepare fixed Lagrangian mass cohorts from the actual slow-wind history.
   Evolve independent radial columns with pressure gradients and a shared hot
   bubble; use a standard compressive artificial-viscosity shock closure.
   Photoionized pressure must enter the force calculation, not only the color.
   Radiative heating/cooling are approximate; no closed total-fluid-energy claim.
3. Compare a small radial refinement and track mass, bubble-energy input/work,
   velocity, minimum cell width and ionization behavior. If this inexpensive
   model fails, preserve its evidence and do not promote it for rendering.
4. Calculate line coefficients with PyNeb. Use the same gas density in absorption,
   ionization and emission; evaluate H/He competition and O/N stages locally.
   Keep intrinsic line luminosities separate from fixed display gains.
5. Sample the continuous 3D density/line field without particle beads or contour
   overlays. Resulting arcs are projected compressed gas, not drawn curves.
   No promise of genuine resolved 3D filament instability at this cost.
6. First make three representative frames, then a low-resolution full32s proof
   of concept. Reuse actual No.004/003 source functions (see linked audit).
   Do not start1080×1920/30fps production in this stage.

Cost target: seconds–minutes for atomic tables and startup checks; benchmark the
radial-sector evolution before extending it (aim <30min per accepted small run);
roughly10–30min for a288–360px-wide preview based on previous measured renders.
These are estimates, not promises. No high-resolution3D RHD/MHD is proposed.

## Interpretation boundary

Calculated: published-track interpolation, conserved cohort masses, approximate
radial force/pressure evolution, direct photon transport, atomic line emissivity,
opacity integration and spectral/optical rendering.

Prescribed: unresolved angular launch statistics, terminal wind speeds, dust
opacity, simple thermal/ionic closures, grey dust reprocessing, pulsation
amplitude, and omitted lateral fluid physics. The source maps remain a rescaled
solar-morphology surrogate for giants, not an AGB radiation-hydrodynamic output.

Presentation: separate atmosphere clock, nonuniform secular time compression,
camera travel, detector sensitivity, line composite and IR display mapping.
None increases the modeled gas mass, ionizing photons or line luminosity.

## Reviewed interface extension and its limit

The calculated contact acceleration gives intervals of positive instability
growth. A local incompressible estimate, `γ_l² = A a sqrt[l(l+1)] / R`, diagnoses
susceptibility. [Sharp1984](https://doi.org/10.1016/0167-2789(84)90510-4), §4 and
Table1, separates linear growth from nonlinear3D curtains/bubbles and discusses
viscosity, compressibility and geometry. His departures from linear behavior
near0.1–0.4wavelength do **not** predict our final amplitudes.

The prototype tests a bounded RT-informed contact closure. Seeds are the fine
AGB loading left unresolved by the coarse pressure grid; physical-length growth
uses the measured acceleration and a density-based Atwood estimate. A smooth
bound near0.1wavelength and exponential radial extension are **prescribed
nonlinear closures**, not solved fluid dynamics. There is no fresh random field,
drawn filament, extra mass or independent emitter. Material volumes, density,
absorption, photon inventory and line radiation are recomputed on the displaced
cells. Bubble volume is preserved; the additional kinetic-energy allocation is
checked against1% of wind work. Lateral backreaction, mode coupling, ablation and
molecular knot formation remain unresolved. This is an approximate
instability-informed reconstruction, not a validated3D instability simulation.

The current preview ends at16kyr for the mature PN, before the old45kyr fading
halo. Its trajectory is calculated to20kyr. The same stellar track and intrinsic
luminosities remain in use. The final source is a hot post-AGB core, not yet a
cold DA white dwarf. A fixed100× narrowband-to-continuum detector sensitivity
makes the low-mass nebula readable; intrinsic luminosities are reported raw.

## Corrections discovered by validation

The first new radial trial used instantaneous H equilibrium. At500yr it ionized
0.308M☉ although even the **entire** emitted photon count could ionize at most
0.122M☉ without dust or recombination. Its snapshot bookkeeping closed but its
front history was physically impossible. Both96/144cell trial outputs remain
recorded and rejected. The selected calculation instead conserves an evolving
ion inventory and includes recombinations, dust, escape and new ionizations.

The independent radiation rays cap the approximate He allowance by the actual
number of photons above24.59eV. The coarse force calculation retains a fixed12%
grey allowance; this early heating difference (bounded by12% of QH) is an
additional approximation, not fully coupled multi-species RHD. Ionized/neutral
temperatures are fixed10kK/100K. The hot-core spectrum is Planck/grey, not NLTE.
PyNeb improves the **line coefficients**, not those underlying closures.

The sharpened preview also exposed an undersampled-OTF axial cross around the
unresolved core. A nonnegative Airy PSF integrated8×8 within each central pixel,
with normalized Moffat wings and the original No003 convolution methods,
corrects detector sampling. A hue-preserving display curve avoids whitening all
line colors. These changes do not remove stellar extinction or add luminosity.

The encoded base preview exposed a second imaging problem: an optical
photospheric exposure could not show the dust's infrared reprocessing during
envelope obscuration. The corrected preview uses an independently exposed,
warm false-color2–40μm channel during9–18s (dimensionless exposure8000 at full
weight). Outside that smoothly blended window the relative IR gain remains0.12.
It is a multispectral composite, not true visible color or new emitted energy.
An intrinsic point-source exposure reference preserves the compact core's
readability after emergence; its calculated foreground extinction is still
applied. Line-power camera framing replaces scattering-halo framing. Same-epoch
mass, photons and intrinsic line/IR power are checked against the base render.


---

# Final targeted validation — GIANT CELLS selection

**Decision: accept the existing GIANT CELLS movie without another render.**

The selected artifact is 1080 × 1920, 30 fps, 32 s, H.264/Rec.709 limited-range,
with 960 unique frames. Its SHA-256 remains
`e72f76440db97cb60e9c47a62c426d9107bb2629890677313c0e0902df998fcb`.
This review changes documentation and packaging, not stellar evolution, gas,
ionization, radiative inputs, camera, exposure or the selected movie.

## Surface implementation and the two questioned intervals

The evolved surface is generated from persistent sites on the sphere. Site
births, angular sizes, lifetimes, tangent drift and elliptic deformation are
deterministic for the frozen seeds. Living sites retain their birth sizes;
current stellar radius is absent from their coordinate transform. Birth
statistics do respond to the stellar parameters. Main, weaker larger and
shorter-lived fine populations compete locally to form bright interiors and
dark boundaries. New sites subdivide territories and fading sites are replaced
continuously. This represents statistical turnover, not calculated fluid flow.

The giant field has no cube-face atlas, four fixed chart regions or repeated
solar-image magnification. Its periodic longitude and common pole values pass
the existing continuity checks. The actual No.004 opening hands off before
the questioned intervals: surface-model weight is exactly one at every newly
checked epoch from 5.7 through 8.0 s. There is no source-model switch there.

**5.7–6.3 s:** the accepted stellar trajectory passes the RGB maximum and warms
rapidly afterward. Evaluated temperatures change from about 3,107 K at 5.9 s
to 4,285 K at 6.3 s, including the existing representative radial pulsation.
The normalized bolometric field changes relatively modestly, while the optical
band is much more temperature sensitive. Cell birth/replacement and the
smoothly changing prescribed contrast also contribute. All surviving site
IDs retain exactly the same angular birth sizes across this interval.

**7.7–8.0 s:** cooling and increasing apparent stellar diameter make the
pattern warmer, more contrasty and larger in pixels. Newborn-site statistics
also change with the evolving gravity. This is the declared projection and
statistical response; living coordinates are not globally zoomed. Surviving
sites again retain exactly the same angular birth sizes. The unwrapped
pulsation phase passes through zero at eight seconds without resetting.

| Playback | Effective T including pulse | Bolometric field RMS | Optical CIE-Y field RMS |
| --- | ---: | ---: | ---: |
| 5.7 s | 3,646 K | 22.34% | 38.09% |
| 5.9 s | 3,107 K | 22.26% | 43.59% |
| 6.3 s | 4,285 K | 21.13% | 31.20% |
| 7.7 s | 4,006 K | 22.13% | 34.79% |
| 8.0 s | 3,422 K | 22.44% | 40.60% |

These are full-sphere source-field statistics before limb response, opacity,
PSF, exposure and tone mapping. They are not contrasts measured from a real
giant atmosphere. The camera-registered encoded patches have consecutive-frame
correlations of at least 0.962 and 0.971 in the two respective intervals.
Direct inspection of all 31 consecutive encoded frames shows continuous
structure, with no texture replacement or four-region reversal.

## Luminosity, spectrum and clocks

- All 1,920 saved physical samples satisfy
  `4πR²σT_eff⁴ = L_adopted` to a maximum relative difference of 1.12 × 10⁻¹⁵
- Eight newly evaluated full-resolution statistical maps satisfy
  `⟨f⟩ = 1`, `T_local⁴ = T_eff⁴ f` and surface-integrated luminosity to
  within 1.34 × 10⁻¹⁵ of the adopted track luminosity
- The original grey limb operator, independently integrated over emission
  angle, reproduces its supplied monochromatic flux to floating-point precision
- Broadband interpolation/ratio emission was compared with direct spectral
  integration at five local-temperature quantiles and five emission angles
  per epoch. Maximum CIE-Y difference is 0.0655%, within the declared 1%
  numerical-consistency check. This does not validate a grey atmosphere
- All local temperatures at these epochs lie inside the existing spectral LUT
- Physical time is strictly monotonic; global pulsation phases match the
  original formula exactly. Camera, source and statistical atmospheric clocks
  remain separate

The initial, explicitly documented No.004 reference calibration is retained.
This is the accepted calibrated MIST rendering track, not a newly evolved
stellar model. The convection closure does not feed back into the accepted
wind or ionizing-photon calculation.

## Encoded six-second transition and complete film

Fresh measurements from the final encoded frames give:

| Frame / PTS | Encoded limb diameter |
| --- | ---: |
| 179 / 5.966667 s | 814.25865 px |
| 180 / 6.000000 s | 813.50762 px |

The change is **−0.75103 px**. The former ~143.84 px rollback has not returned.
Calculated projected radii match the preceding selected production exactly.

The complete encoded movie was decoded again. Frame count, exact 30 fps PTS,
frame order, absence of duplicated frames, source-chunk identity and the
1,920-sample physical checks pass. The 540 untouched decoded frames, including
the entire 16–32 s compact-core/nebula sequence, remain identical. Existing
mass/photon closure errors stay below 1.9 × 10⁻¹⁵ and bubble-energy closure
below 8.5 × 10⁻¹⁵. No limits were relaxed.

Visual review used actual decoded 0.5 s storyboards spanning the complete film,
31 consecutive 30 fps frames in the questioned intervals, native frames
179–180 and the new cover. Earlier native consecutive-frame join strips cover
all replacement boundaries. This is an encoded-frame review, not a claim of
real-time playback inspection. No new camera/exposure jump, texture reset,
four-region pattern or geometry rollback was identified. The accepted faint
emergence interval and comparatively smooth nebula are preserved.

The maximum non-text saturated-channel fraction is 2.18 × 10⁻⁷: no broad
surface/nebula clipping was found. HDR sensitivity, false-color lines, warm
infrared presentation, PSF and exposure remain display choices, not additional
physical luminosity.

## Physical limits of this acceptance

Flux normalization, deterministic motion and image continuity are necessary
consistency checks. They do not establish realistic convection by themselves.
Cell morphology, drift, splitting, lifetime, contrast and the representative
clock are prescribed. Fit extrapolation at low gravity, angular-resolution
floors and the scale hierarchy are documented. No convective velocity field,
heat transport, shocks, nonadiabatic envelope or global 3D RHD is solved.

The spectrum is Planck/grey with fixed limb coefficient `u = 0.6`; molecular
opacity, atmospheric extension and non-LTE giant/hot-core transfer are not
calculated. Nebular gas temperatures, wind speeds, angular structure and ionic
closures retain their existing approximations. A bright nebula is not assured
for the actual Sun. Acceptance means the chosen physics-informed visualization
is internally consistent and the specific production artifacts are resolved.

## Reproducible evidence

- `validation/publication-targeted-science.json`
- `validation/publication-radius-179-180.json`
- `validation/publication-encoded-photometry.json`
- `validation/giant-cells-final-checks.json`
- `output/publication-review/critical-*-actual-encoded.jpg`
- `output/publication-review/film-*.jpg`
- `production/RUNTIME.json` and the frozen production/surface revisions

The cover is actual encoded frame 780 at 26.0 s, exported at native resolution
without added graphics. The curated release omits rejected experiments,
unselected caches, papers, machine-specific paths and private records. Large
No.004 solar-source datasets remain external dependencies with exact hashes
and documented setup; only the two selected nebular caches are included in
the complete bundle. Nothing is published externally by this workflow.

## Included encoded-frame evidence

- [5.7–6.3 s consecutive frames](METHODS.md)
- [7.7–8.0 s consecutive frames](METHODS.md)
- [0–8 s](METHODS.md) · [8–16 s](METHODS.md) · [16–24 s](METHODS.md) · [24–32 s](METHODS.md)


---

# Scientific references and their actual roles

This file covers the original prototype and the current research preview.
The v3 additions below and `LITERATURE-IMPLEMENTATION-V3.md` describe the
selected calculation; earlier approximate coefficients are historical only.

- **Choi et al. (2016), MESA Isochrones and Stellar Tracks I** — [paper](https://arxiv.org/abs/1604.08592), [official grids](https://mist.science/model_grids.html). Actual numerical stellar properties in the preview are from the nonrotating solar-metallicity MIST v1.2 1 M☉ track, calculated with MESA revision 7503. Download URL and SHA-256 are recorded separately. This is reused published evolution, not a fresh MESA run.
- **Dotter (2016), MESA Isochrones and Stellar Tracks 0** — [paper](https://arxiv.org/abs/1601.05144). Equivalent-evolutionary-point sampling context. EEP interpolation is not a new high-cadence calculation of atmospheric pulses.
- **Miller Bertolami (2016), A&A 588, A25** — [paper](https://arxiv.org/abs/1512.04129), [CDS model description](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/588/A25?format=html). Comparison of low-mass post-AGB heating times and mass-loss/mixing uncertainties. The published 1 M☉ Z=0.02 example has a different final mass/timescale from the adopted MIST sequence. The CDS numeric download timed out; no values were fabricated and no track splice was made.
- **Gesicki, Zijlstra & Miller Bertolami (2018)** — [paper](https://arxiv.org/html/1805.02643). PN visibility near the Sun's initial mass, declining envelope density, AGB-exit phase dependence and the difference between an ionized and a bright nebula. This work motivates explicitly treating the Sun's outcome as uncertain, rather than guaranteeing a Helix-like object.
- **Schröder & Smith (2008)** — [paper](https://arxiv.org/abs/0801.4031). Independent solar giant-phase study demonstrating sensitivity to mass-loss assumptions; not a source of replacement numerical points for the MIST movie.
- **Höfner & Olofsson (2018)** — [review](https://doi.org/10.1007/s00159-017-0106-5). Physical role of pulsation, convection, atmospheric levitation and radiation pressure on dust in slow AGB winds. Our terminal wind model does not claim to solve those mechanisms.
- **Freytag, Liljegren & Höfner (2017)** — [3D RHD paper](https://arxiv.org/abs/1702.05433). Large convective cells and radial pulsation in AGB atmospheres. Qualitative constraints only; no CO5BOLD snapshots have been imported and no AGB RHD calculation was run.
- **Höfner & Freytag (2019)** — [dust-forming 3D AGB model](https://arxiv.org/abs/1902.04074). Clumpy dust and non-spherical shocks motivate inhomogeneous launch conditions. The preview's specific clumps are explicitly assumed, not borrowed simulation results.
- **Kwok, Purton & Fitzgerald (1978)** — [original interacting-winds paper](https://adsabs.harvard.edu/pdf/1978ApJ...219L.125K). Fast wind interacting with prior slow wind, rather than an impulsive stellar explosion.
- **Perinotto et al. (2004)** — [radiation-hydrodynamics study](https://www.aanda.org/articles/aa/pdf/2004/06/aah4695.pdf). Physical comparison for wind compression, ionization, rim/shell structure and the limitations of a shared-bubble sector approximation. The original numerical code/output is not used.
- **Storey & Hummer (1995)** — [recombination calculations](https://doi.org/10.1093/mnras/272.1.41). Case-B recombination physics. The original prototype used fixed10,000K Hα/Hβ coefficients. The selected v3 uses actual PyNeb table lookups over electron temperature/density, while retaining a prescribed10,000K gas temperature.
- **Luridiana, Morisset & Shaw (2015), PyNeb** — [paper](https://arxiv.org/abs/1410.6662), [official documentation](https://morisset.github.io/PyNeb_devel/). PyNeb1.1.32 was installed locally and run to generate the selected v3 emissivity tables, with recorded atomic files and interpolation checks. It calculates atomic emissivities, not the full photoionization or energy balance.
- **Jacob, Schönberner & Steffen (2013)** — [visibility and expansion study](https://arxiv.org/abs/1307.6189). Expansion/visibility times and distinction between material motion and ionization-front appearance.

- **Ivezić & Elitzur (1997), dusty radiative-transfer scaling** — [paper](https://arxiv.org/abs/astro-ph/9612164). Context for the importance of opacity, spectral reprocessing and radiative equilibrium in an AGB dust shell. The preview uses a much simpler grey diffusion/escape approximation, **not DUSTY or the paper’s exact numerical solutions**.

## Reused code, numerical color data and fonts

SiderAtlas No.004 and curated No.005: original MIT code, retained notice. The copied video encoder is unchanged. The spectral integration uses the unmodified CIE 1931 2° color-matching-function table; CIE attribution/CC BY-SA 4.0 rights are retained. Latin Modern fonts remain under their GUST Font License. No NASA/Hubble/Helix/Ring footage or imagery was used. No unconfirmed Bifrost data redistribution is involved in the new package.

- **Mathur etal.(2011), granulation in red giants** — [paper](https://doi.org/10.1088/0004-637X/741/2/119). Pressure-height and turnover-time scalings motivate the evolved sampling scales. The No004 source morphology remains a surrogate; these scalings do not validate a new giant convection calculation.
# Research-led v3 additions

Primary literature and the exact implementation decisions are mapped in
[LITERATURE-IMPLEMENTATION-V3.md](METHODS.md).

- Miller Bertolami2016, A&A588 A25 — [published paper](https://doi.org/10.1051/0004-6361/201526577), Table3; sensitivity of the post-AGB scenario, not a replacement mass history
- Gesicki, Zijlstra & Miller Bertolami2018, Nature Astronomy2,580 — [paper](https://doi.org/10.1038/s41550-018-0453-9); faint near-threshold Sun-like PN, not its maximum-nebula luminosity imposed on our gas
- Höfner & Olofsson2018, A&A Review26,1 — [review](https://doi.org/10.1007/s00159-017-0106-5); dust-driving criteria, wind/condensation ranges
- Villaver, Manchado & García-Segura2002, ApJ581,1204 — [paper](https://arxiv.org/abs/astro-ph/0208323), Figs5/11; ionization pressure, rims and faint halos; different older remnant mass
- O'Dell, Henney & Ferland2005, ApJ630,996 — [paper](https://doi.org/10.1086/430803); dense molecular cores needed for Helix knots, not copied into this model
- Jones & Boffin2017, Nature Astronomy1,0117 — [review](https://doi.org/10.1038/s41550-017-0117); many complex morphologies require binary shaping; no binary is introduced
- Luridiana, Morisset & Shaw2015, A&A573 A42 — [paper](https://doi.org/10.1051/0004-6361/201323152); actual PyNeb n-level/Case-B emissivity machinery, not ionization/thermal equilibrium
- Ferland et al.2017, RMxAA53,385 — [paper](https://arxiv.org/abs/1705.10877); reference for what full Cloudy calculations would include; Cloudy is **not** run here
- Toalá & Arthur2014, MNRAS443,3486 — [paper](https://doi.org/10.1093/mnras/stu1360), Eqs1–4/§4.1/AppendixA; pressure, instability and shadowing mechanisms, not borrowed simulation geometry
- Höfner & Freytag2019, A&A623 A158 — [paper](https://arxiv.org/abs/1902.04074), Table1/Figs2–4; clumpy shocked AGB atmospheres; that model does not drive a wind
- Freytag, Liljegren & Höfner2017, A&A600 A137 — [paper](https://arxiv.org/abs/1702.05433); AGB convection/pulsation scale context
- Freytag & Höfner2023, A&A669 A155 — [paper](https://arxiv.org/abs/2301.11836); dust-driven clumpy winds; qualitative source conditions, not imported3D output
- García-Segura et al.2006, ApJ646 L61 — [paper](https://arxiv.org/abs/astro-ph/0606205); fast-wind shutdown/fragmentation experiment, whose abrupt shutdown and Helix identification are not adopted
- Sharp1984, Physica D12,3–18 — [review](https://doi.org/10.1016/0167-2789(84)90510-4); linear versus nonlinear RT scope, not a prediction of saturation amplitude
- Verner et al.1996, ApJ465,487 — [paper](https://doi.org/10.1086/177435) and [author's data](https://www.pa.uky.edu/~verner/photo.html); actual photoionization cross-section fits
- Historical RR and low-temperature DR — [author's tables with references/formulas](https://www.pa.uky.edu/~verner/rec.html); Aldrovandi & Péquignot1973, Shull & van Steenberg1982, Arnaud & Rothenflug1985, Nussbaumer & Storey1983/1986/1987

The generated PyNeb table records exact data filenames/version in
`data/atomic/pyneb-v3.json`: H recombination Storey & Hummer1995; OIII collision
strengths Storey, Sochi & Badnell2014; NII Tayal2011; A-values Froese Fischer
et al.2004 and Storey & Zeippen2000. Both density dependence and collisional
de-excitation come from those tables. Temperature remains prescribed.

Downloaded papers/packages are local research inputs, excluded from publication.
Their presence establishes provenance, not a redistribution licence. The selected publication version is linked in the episode README.

## Selected statistical giant-star surface

- Tremblay et al. 2013, A&A 557, A7 — [primary paper](https://arxiv.org/abs/1307.2810); characteristic granulation size and intensity-decay fits, with explicit grid-domain flags and separately identified low-gravity continuation
- Ludwig & Kučinskas 2012, A&A 547, A118 — [primary paper](https://doi.org/10.1051/0004-6361/201220264); representative red-giant contrast and spatial-scale context, not an imported atmosphere or raw simulation sequence
- Freytag, Liljegren & Höfner 2017, A&A 600, A137 — [primary paper](https://arxiv.org/abs/1702.05433); giant/fine-structure hierarchy and atmospheric-timescale context; its convection, shocks and molecular radiative transfer are not solved here

The implemented statistical model, fit applicability, spatial floors,
representative clock, luminosity normalization and native encoded checks are
in [GIANT-CELLS-MODEL.md](METHODS.md). No additional RHD dataset is
redistributed or claimed as a validated source field.

