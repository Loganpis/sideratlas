# Technical methods and attribution

Selected final-model documentation, consolidated without changing the methods.

---

# Physical foundation and scope — v2.2

This is a scientific visualization of a representative ZZ Ceti star. It is not
an observed star reconstruction, a telescope image, a new radiation-hydrodynamic
solution, or a nonlinear calculation of convection interacting with pulsations.

The preserved WDEC structure has M=0.5931574054 solar masses,
R=9018.13376175 km, Teff=11290.0200413 K and log(g)=7.985750095 (cgs).
The 435-point parametrized C/O-core structure and six GYRE 9.1.1 adiabatic
interior eigenvalues are unchanged. The modes include gravitational-potential
perturbations. Selected (ell,k,m) are (1,4,0), (1,6,1), (1,8,-1), (2,6,1),
(2,8,2), (2,10,-2). Their periods at infinity are approximately
315.813,397.376,488.202,247.476,286.946,348.720 seconds.

The model uses real standing spherical harmonics, not radial eigenfunction
nodes painted onto the surface. Surface flux amplitudes, m assignments,
phases, equal mode variance and the selected excitation scenario are prescribed.
Adiabatic GYRE does not predict the emergent temperature amplitudes, driving,
saturation or a validated convective transfer function.

F/F0 = 1 + sum_j A_j H_j(n) cos[2pi(t-delay)/P_j + phase_j]

Teff = T0 (F/F0)^(1/4)

The original temperature field has linear sphere/time RMS 40.234 K. Its physical
amplitudes and atmosphere cache are unchanged; only the display mapping is
stronger. A fixed +/-160 K color range is checked at every ray. No extrapolated
atmosphere temperatures, deformation, rotation, spots or flares are introduced.

## Radiation and gravity

Saved Koester pure-H 1D LTE spectra are integrated with CIE color matching
functions at 0.25 nm from 360 to 830 nm. Claret DA-3D broadband limb fits supply
an approximate angular law. This hybrid is not matched angle- and wavelength-
resolved 3D time-dependent transfer. The direct Schwarzschild ray map retains
weak bending, travel delays and spectral frequency shifts (g^5 for wavelength
specific intensity). Surface z=0.0000971155; GYRE local frequencies are redshifted
once. These small effects are not visually magnified.

## New cellular layer

See [CONVECTION.md](METHODS.md) for the independently calibrated statistical
surrogate. It modulates grey intensity, not calculated mode temperatures. It
is a separate visualization of convection statistics, with an explicitly
magnified spatial scale and a different display clock. Multiplication by the
mode radiance produces a coherent surface image but is not a calculated
physical interaction between convection and oscillations.

## Camera and time

The preserved B camera eases slowly from inclination 55 to 70 degrees and
azimuth -25 to +25 degrees; roll=-20 degrees. It does not represent stellar
rotation. The star remains centered on black, with only the three SiderAtlas
labels. The 32-second film spans 1721.6 seconds of the mode clock (53.8x faster).
The cellular clock advances 0.05 local seconds per display second (20x slower),
with local ray delay included. The two clocks deliberately reveal otherwise
incompatible timescales; this is not one simultaneous physical time sequence.

The copied fixed-observer light curve (inclination 55 deg, azimuth 0) contains
only the original physical mode radiation. It excludes camera motion, enlarged
cell flicker, false colors, tone mapping and halo gains. It must not be inferred
from the enhanced movie's total pixel brightness.

## Retained numerical evidence and limits

The original GYRE integration refinement changes the six local periods by at
most 6.84e-8 s between the fine and finer integrations on the same 435-point
background. The supplied WDEC regression test passed. Thinning the equilibrium
background to 218 points changes a selected period by as much as 1.0173 s.
That is a background sensitivity test, not a full stellar-structure convergence
proof or physical frequency uncertainty. The parametrized composition profile,
mode excitation and static atmosphere remain significant assumptions. The root
`validation/mode-convergence.json` retains the per-mode evidence; no newly
converged stellar equilibrium or nonadiabatic driving calculation is claimed.


---

# Statistically constrained convection visualization

## What this is

A deterministic-seed, evolving cellular surrogate informed by DA radiation-
hydrodynamic literature. It is **not** a validated hydrodynamic flow, a recovered
CO5BOLD snapshot, a prediction of optical granulation, or output of the failed
local convection solver. No failed-flow outputs are imported or packaged.

Tremblay et al. (2013, A&A 557 A7), Table 1, provides neighboring DA models.
The interpolation retained in `validation/literature-constraints.json` suggests
a characteristic granulation length 0.8523 km and vertical bolometric intensity
RMS 17.2679%. Equation 16 gives an approximate intensity e-folding time 0.1105 s,
with intrinsic scatter and domain limitations. This is an empirical estimate,
not a directly measured lifetime for this exact star. The log(g)=8 rows at
11005 and 11529 K supply an interpolated bright-area fraction of 54.89%.
Neighbor interpolation is not a new matched RHD model.

These are bolometric intensity statistics. They do not uniquely imply local
Teff fluctuations, optical contrasts, emergent angle-dependent spectra, gas
velocities or convective fluxes. We do not make those identifications. Conditional
Mach-based velocity estimates in the retained audit are context only; the
surrogate does not solve or reproduce a fluid velocity field.

## Cellular construction

About 8192 finite-lived centers cover the sphere. A seeded stationary Poisson
birth process, lognormal lifetimes, smooth sin-squared growth/decay weights and
prescribed tangent drifts create continuously evolving spherical weighted
Voronoi cells. Two smallest weighted squared chord distances define an interior-
versus-boundary gap. A fixed analytic profile makes bright interiors and dark
lanes; no noise image, random per-frame texture or magnetic spot layer is added.
The solar granulation driver of Gudiksen & Nordlund (2005) is a geometric
precedent only. Its solar numerical scales are not used as DA constraints.

A profile power and offset are fitted to bright area and contrast. A dimensionless
time unit is fitted to the empirical e-fold time. Held-out directions and epochs
check contrast and decay. The fitted intensity field G is positive and normalized
to unit whole-sphere area mean at each atlas time. This removes global source
mean drift; it is not observer auto-exposure or a prediction of true granular
photometric variability. The finite birth catalog does not wrap periodically.

The held-out calibrated RMS is about 17.24%, bright fraction 54.90% and decay
0.1108 s. These are statistical calibration checks, **not validation of RHD**.
Spatial power peaks vary across finite patches; a median angular characteristic
wavelength 0.05210 rad defines the spatial mapping. At this star's radius it
corresponds to about 470 km: an explicit **551.3x magnification** of the 0.852 km
reference scale. Physical granules would be only 0.05 native pixels. The full-
star visual cell count is consequently much smaller than the physical count.

## Transfer, time and limitations

The existing calculated mode radiance is multiplied by G as a grey bolometric-
contrast proxy. There is no local cellular temperature calculation or new
atmosphere spectrum. The same mean limb law is used for all cells. Nonadiabatic
phase response, mode-dependent granulation statistics, fluid velocity, Doppler
signals and nonlinear mode/convection coupling are absent. Birth, drift,
lifetime distribution and morphology remain prescribed. Matching a few moments
does not establish the true spatial spectrum, probability distribution or
correlation function of a DA atmosphere.

The local convection display clock advances 0.05 s per displayed second; mode
physics advances 53.8 s. Light travel delays use each layer's own clock. Temporal
interpolation is at 0.005 local seconds on a 1025x2048 spherical atlas. Spatial
bilinear interpolation and the optical PSF smooth the narrowest boundaries;
`validation/cellular-numerics.json` quantifies this approximation. The intensity
contrast is additionally enhanced in display log-luminance by exponent 4, after
normalized optical filtering. It is not 4x the physical temperature fluctuation.

## Numerical display losses

Independent equal-area samples of the 1025x2048 atlas have about 16.63%
relative intensity RMS, versus 17.25% in the direct analytic cell field at the
same test time. Bilinear spatial interpolation smooths the narrowest lanes.
Halving the local time step changes intensity by 0.00089 RMS. In a native-pixel
image patch, after a one-pixel Gaussian proxy for core smoothing, the absolute
intensity RMS difference from direct evaluation is 0.00903. These are measured
numerical approximation errors, not proofs of hydrodynamic convergence.

The prescribed mean center lifetime is 3 dimensionless units, about 0.741 local
seconds after calibration. This is distinct from the approximately 0.11-second
intensity correlation decay: growing boundaries and changing cell dominance
can decorrelate intensity before a center disappears. Neither center lifetimes
nor tangent drift speeds are independently measured DA fluid quantities.


---

# Fixed display and luminous imaging

Every changing large-scale structure follows the six preserved mode phases.
There is no artificial surface deformation or displacement of eigenfunctions.
Fine cellular changes follow the reproducible statistical model described in
CONVECTION.md. They are not extra GYRE modes or calculated fluid trajectories.

## Pulsation diagnostic

A fixed +/-160 K scale encodes the calculated mode temperature perturbation:

u = asinh(clamp(DeltaT,-160,160)/4) / asinh(160/4)

The palette is fixed blue to blue-white (all RGB anchors are in config.json).
The palette blend strength is 96%. It is false color, not the true visible color of a resolved white dwarf. The
physical temperature amplitudes, mode frequencies and phases are unmodified.
After filtering, a gain 2.7*u is added to log-luminance to reveal pulsation
interference. There is no frame-dependent histogram or automatic exposure.

## Optical calculation versus display enhancements

The physical spectral RGB field is filtered with a normalized polychromatic
circular-aperture Airy OTF (first zero at 550 nm =3.8 pixels at 1080 width),
Gaussian jitter sigma=0.6 pixels, and a broad Moffat beta=2.5 component with
alpha=85 pixels. Broad-scatter fraction is 0.12. Parameters scale with image
resolution. The prescribed virtual instrument is not any named telescope,
absolute angular resolution or physical distance.

Before applying display enhancements, padded linear CIE-Y flux conservation
is measured for the core/halo mixture. That check does not assert conservation
of the final diagnostic RGB image. Cropping, false colors, contrast, gain,
nonlinear tone mapping and display gamut compression change displayed flux.

A sixth source channel holds Y*G for the grey granular proxy. The core-filtered
ratio filtered(Y*G)/filtered(Y) modulates display log-luminance with exponent 4.
It is not an additional thermal spectrum or physical temperature amplitude.
The core display is 0.94*sigmoid(log[1.7*(1-0.12)*Ycore] + 2.7*u*coverage
+ 4*log(Gfiltered)*coverage). A fixed sixfold display gain makes the broad
halo more visible. Halo RGB comes from the calculated DA spectrum; the surface
false-color palette does not color the outer halo. The remaining luminance
headroom blends the mapped halo into the core without washing out the texture.
A fixed luminance-preserving gamut compression produces valid sRGB values.

The halo is an imaging/display effect, not a physical stellar atmosphere
extending hundreds of kilometers beyond the surface. The final image is a
scientific visualization, not an apparent telescope image or calibrated
photometric measurement. The caption discloses these choices outside the Reel.

## Sampling, motion and format

Native final frames are 1080x1920, with 2x2 source sampling, 8x8 limb sampling,
analytic pixel/circle coverage and two Gauss temporal samples. Cellular atlases
use 1025x2048 nodes and 0.005 local-second steps. The normalized mean field is
interpolated on the sphere with periodic longitude, not wrapped at the movie
boundary. B's slow camera is unchanged. Three labels use Latin Modern fonts.

The preview is eight seconds at 540x960/15fps, using the first eight seconds of
the same mode and camera schedule. The final film is one continuous 32-second,
30fps sequence with 960 distinct time samples. It is not a forced seamless
loop or four repetitions. H.264 High, 4:2:0, Rec.709 limited-range metadata,
explicit color conversion, high-quality compression and faststart are used.


---

# Primary references and input provenance

## Interior model and oscillations

- Bischoff-Kim & Montgomery (2018), *WDEC: A Code for Modeling White Dwarf Structure and Pulsations*, AJ 155, 187 — [DOI](https://doi.org/10.3847/1538-3881/aab70e)
- Townsend & Teitler (2013), *GYRE: an open-source stellar oscillation code based on a new Magnus Multiple Shooting scheme*, MNRAS 435, 3406 — [DOI](https://doi.org/10.1093/mnras/stt1533)
- Unmodified [GYRE 9.1.1](https://github.com/rhdtownsend/gyre/tree/v9.1.1), commit `d9ebbd1f3ca941042c5c9bd667102d9cf30c035b`; its supplied `models/wdec/m0.593/evolved` equilibrium and original WDEC controls are retained
- Wu & Goldreich (1999), *Gravity Modes in ZZ Ceti Stars. II. Eigenvalues and Eigenfunctions*, ApJ 519, 783 — [DOI](https://doi.org/10.1086/307412)
- Van Grootel et al. (2012), *The instability strip of ZZ Ceti white dwarfs. I. Introduction of time-dependent convection*, A&A 539, A87 — [DOI](https://doi.org/10.1051/0004-6361/201118371)
- [GYRE's documented convection closures](https://gyre.readthedocs.io/en/stable/ref-guide/osc-equations/conv-effects.html): used for a diagnostic sensitivity study, not as validated DAV convective driving

## Atmosphere and granulation

- Koester (2010), *White Dwarf Spectra and Atmosphere Models*, Mem. S.A.It. 81, 921 — [preprint](https://arxiv.org/abs/0812.0482)
- Published Koester pure-H flux spectra, distributed by the [SVO theoretical spectra service](https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=koester2); exact files, URLs, parameters and hashes are in the retained `data/*spectra-provenance.json` files
- Tremblay & Bergeron (2009), hydrogen line profiles, ApJ 696, 1755 — [DOI](https://doi.org/10.1088/0004-637X/696/2/1755)
- Claret et al. (2020), DA-3D limb darkening, A&A 641, A157 — [DOI](https://doi.org/10.1051/0004-6361/202038436), [CDS tables](https://cdsarc.cds.unistra.fr/ftp/J/A+A/641/A157/); saved four-term coefficients supply the broadband angular approximation
- Tremblay et al. (2013a), *Pure-hydrogen 3D model atmospheres of cool white dwarfs*, A&A 552, A13 — [preprint](https://arxiv.org/abs/1302.2013), [DOI](https://doi.org/10.1051/0004-6361/201220813)
- Tremblay et al. (2013b), *Granulation properties of giants, dwarfs, and white dwarfs from the CIFIST 3D model atmosphere grid*, A&A 557, A7 — [preprint](https://arxiv.org/abs/1307.2810); Table 1 supplies the neighboring RHD statistics used to constrain the statistical cellular approximation, not a matched new RHD solution
- Tremblay et al. (2013c), *Spectroscopic analysis of DA white dwarfs with 3D model atmospheres*, A&A 559, A104 — [preprint](https://arxiv.org/abs/1309.0886), [DOI](https://doi.org/10.1051/0004-6361/201322318); informs the distinction between full RHD, mean 3D structures and 1D spectra
- Hollands et al. (2026), *Any colour you like: the diverse chromaticities of white dwarfs* — [DOI](https://doi.org/10.1093/rasti/rzag069); investigated as a possible angle-resolved spectral improvement, but no matching private/request-only atmosphere grid is claimed or used

- Kupka, Zaussinger & Montgomery (2018), *Mixing and overshooting in surface convection zones of DA white dwarfs: first results from ANTARES*, MNRAS 474, 4660 — [DOI](https://doi.org/10.1093/mnras/stx3119); resolved grey-RHD context for the rejected local-flow audit

## Magnetic compatibility

- Tremblay et al. (2015), *On the Evolution of Magnetic White Dwarfs*, ApJ 812, 19 — [preprint](https://arxiv.org/abs/1509.05398), [DOI](https://doi.org/10.1088/0004-637X/812/1/19)
- Gentile Fusillo et al. (2018), *Can magnetic fields suppress convection in the atmosphere of cool white dwarfs? A case study on WD2105−820*, MNRAS 473, 3693 — [preprint](https://arxiv.org/abs/1710.02151), [DOI](https://doi.org/10.1093/mnras/stx2584)
- Rui, Fuller & Hermes (2025), *Supersensitive Seismic Magnetometry of White Dwarfs*, ApJ 981, 72 — [preprint](https://arxiv.org/abs/2410.20557)

These papers motivate excluding imposed magnetic spots from an unchanged
nonmagnetic mode calculation. They do not establish a universal prohibition
on magnetism and pulsation coexisting.

## Propagation, color and optics

- Pechenick, Ftaclas & Cohen (1983), *Hot spots on neutron stars: the near-field gravitational lens* — [DOI](https://doi.org/10.1086/161498); spherical direct-ray propagation context, used here in the very weak white-dwarf field
- [CIE 1931 2-degree color matching functions, official dataset](https://doi.org/10.25039/CIE.DS.xvudnb9b), CC BY-SA 4.0; retained unmodified with its metadata
- [W3C CSS Color 4 conversion definitions](https://www.w3.org/TR/css-color-4/#color-conversion-code), XYZ/linear-sRGB matrix and sRGB transfer conventions
- [Photutils circular-aperture Airy PSF](https://photutils.readthedocs.io/en/stable/api/photutils.psf.AiryDiskPSF.html)
- [Photutils normalized Moffat PSF](https://photutils.readthedocs.io/en/stable/api/photutils.psf.MoffatPSF.html)
- Trujillo et al. (2001), *The effects of seeing on Sersic profiles. II. The Moffat PSF* — [preprint](https://arxiv.org/abs/astro-ph/0109067)

The optics are an explicitly prescribed virtual instrument, not the PSF of
any named telescope. Latin Modern fonts retain the bundled GUST/LPPL notices.
Original source code is MIT; external data retain their individual credits
and terms. No NASA footage, generative texture, decorative turbulence or
invented magnetic map is used.

## Cellular geometry precedent

- Gudiksen & Nordlund (2005), *An Ab Initio Approach to the Solar Coronal Heating Problem*, ApJ 618, 1020 — https://doi.org/10.1086/426063; https://arxiv.org/abs/astro-ph/0407266. Weighted Voronoi granulation driver: geometric precedent only. No solar physical scales, magnetic model or coronal heating solution is used.


---

# Scientific data credits and terms

## CIE color-matching functions

Creator and rights holder: International Commission on Illumination (CIE)

Citation: CIE (2019), *Colour-matching functions of CIE 1931 standard
colorimetric observer*, CIE, Vienna, DOI
[10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b)

[Official dataset landing page](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer)

The original CSV and metadata, downloaded on 2026-10-17, are redistributed
**unmodified**:

- [Official CSV](https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv)
- [Official metadata](https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv_metadata.json)

CSV SHA-256: `fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1`

MD5: `17cca777db64b17170f06f67ce9d3ab7`

Both match the original metadata. The table has 471 rows over 360–830 nm,
with wavelength, x-bar, y-bar and z-bar columns. Original source: CIE
018:2019, Table 6.

The metadata identifies [Creative Commons Attribution-ShareAlike 4.0
International](https://creativecommons.org/licenses/by-sa/4.0/). Retain this
attribution, DOI, license link and modification status when sharing the
dataset. Adaptations of the dataset remain subject to those terms. These
data have their own license and are not relicensed under the repository's
software MIT license. CIE does not endorse this project.

The code interpolates the original matching functions and integrates
observer-frame spectral radiance to XYZ, scaled by 683 lm/W. The W3C
[color-conversion definitions](https://www.w3.org/TR/css-color-4/#color-conversion-code)
supply the standard XYZ/linear-sRGB matrix and sRGB transfer convention.
False-color diagnostics, the virtual PSF and the fixed display response are
separate declared visualization operations, described in MODEL.md.

## Atmosphere, interior inputs and fonts

Published Koester pure-H atmosphere spectra retain their authorship and SVO
service provenance. SVO distributes its Koester models with the author's
permission. The original Claret DA-3D table and ReadMe are checksum-identified external inputs in this consolidated collection. Exact
input URLs, parameters and SHA-256 values are in the `data/` provenance JSON
files. Their source attribution and terms remain distinct from the code license.

The supplied GYRE/WDEC equilibrium and controls are credited to their original
authors and GYRE's GPLv3 distribution. The GPL COPYING and WDEC license notices
are retained under `licenses/`, with the public source revision in SOURCES.md.
The compiled GYRE solver and MESA SDK are not redistributed. Calculated
eigenvalues and selected numerical results are supplied with their provenance.

The granulation assessment uses attributed numerical statistics from
Tremblay et al. (2013), not paper figures, a copied movie or a fabricated
simulation snapshot. Latin Modern fonts retain their bundled GUST/LPPL
notices. Citation and data availability do not imply third-party endorsement.

## v2.2 cellular surrogate

`v220/data/cellular-catalog.npz` is generated by original MIT-licensed project code with seed 220003. It is not a CO5BOLD/RHD dataset. `cellular-calibration.json` records statistical fits to neighboring DA bolometric summaries from Tremblay et al. (2013), with full references in SOURCES.md. Solar weighted-Voronoi work is cited as a geometric precedent only; no external texture, flow snapshot or failed local-RHD result is used.

