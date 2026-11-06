# Technical methods and attribution

Selected final-model documentation, consolidated without changing the methods.

---

# No.004 — golden photosphere and calculated ultraviolet plasma

## What this film is

One 32-second scientific visualization of a representative 1 M☉, 1 R☉, nominal 5772 K G2 V star. It combines **independent local published radiation-MHD models** in prescribed positions on a sphere. It is not a newly solved global magnetohydrodynamic star, a reconstruction of a particular observed Sun, or a literal visible-light telescope image. No new jet trajectories, flame shapes, magnetic lines, particles or spot masks are drawn.

The physical photospheric calculations and full-source geometry from the last photospheric preview are preserved. The underlying ultraviolet calculation from the actual spicule-producing model is reused unchanged. This revision refines fixed display grading and red-orange channel presentation, increases projection quadrature for all viewing angles, and replaces the camera scan with four uniform orbits. Historical rejected coronal composites are not movie inputs.

## Photospheric radiation and magnetic structure

Quiet/network radiation comes from 42 time-dependent Bifrost `en024048_hion` atmospheric states. The already completed inclined LTE H⁻ continuum transfer uses T, electron density and neutral hydrogen, 17 wavelengths over 380–800 nm, five μ angles and four azimuths. A **fixed** initial spectral calibration to TSIS-1 HSRS v2 supplies missing mean line blanketing. It is not a new complete solar atmosphere spectrum. Limb darkening and its local angular variation follow that transfer. The source has 47.6 km horizontal cells; output radiation maps sample at 95.2 km.

The complete photospheric magnetic complex comes from Cheung et al.'s **evolving numerical MURaM outgoing-gray-intensity maps**, paired with magnetic slices. Both magnetic concentrations stay in their source relative geometry. Their dark cores, resolved radial structures and bright concentrations are calculated quantities. Stronger contrast does not create new filaments. Native 192 km sampling cannot resolve the finest solar penumbral filaments. No finer grid or new radiative MHD was run.

The active map's spectral proxy is T_b = 5772 (I_gray/I_quiet)^(1/4), with the measured solar reference spectrum multiplied by a Planck spectral ratio and corrected for optical/bolometric fraction. It is a brightness-temperature approximation, **not the actual sunspot temperature or directional spectrum**. Its limb law is the quiet-star mean. Quantitative facular hot-wall radiation and Wilson-depression visibility have not been solved. The Bifrost network does retain independently calculated inclined transfer.

The quiet atlas uses overlapping fixed local-source transforms: quarter turns, periodic translations and fixed 0–300 s initial phase offsets on six spherical charts. There is no animated procedural noise or source sliding. The one coherent active map has a prescribed smooth outer join. Local models do not exchange magnetic flux, energy or mass across joins. Spot coverage and global magnetic topology are reconstruction assumptions.

## Actual spicule-producing plasma

The public **Bifrost `en096014_gol`** model of Martínez-Sykora et al. includes radiative MHD, conduction, Hall and ambipolar terms. Its ejections arise from the published plasma dynamics and magnetic interactions. The archive is **2.5D**, with x and z structure and all vector components; the third direction is invariant. It does not calculate a 3D population of individual cylindrical jets.

We retain 71 actual snapshots, numbered 300–370. Temperatures and electron densities are read in SI units from byte-verified FITS slabs spanning approximately 0.6–14 Mm above the source photosphere. The native horizontal spacing is 13.89 km and z is nonuniform. Actual header timestamps are used, including unusually close outputs; integer snapshot number is never treated as a uniform clock. Source thermodynamics retain their own effective temperature and LTE ionization treatment. They are not adjusted to 5772 K.

The source validation independently retrieves temperature, density, electron density, horizontal/vertical velocities and magnetic components at snapshot 325. Its diagnostic selects real cool upflowing material without using that selection to draw anything. The detailed heights, speeds, temperatures and field strengths are in `../validation/plasma-source.json`. This checks units and source content; it does not establish a fresh MHD convergence study or identify every bright strand as a jet.

## Ultraviolet radiation calculation

We calculate the **C IV 1548.187 + 1550.772 Å doublet**, separate from the optical continuum. Heated spicule interfaces, transition-region surfaces and other C IV-emitting structures can contribute. Cool spicule cores need not be visible. No Hα, Mg II, Ca II, He II 304 Å, broadband-visible spicule radiation or full corona is claimed.

The scalar emissivity is

`j_CIV = n_e² G(T)` in erg cm⁻³ s⁻¹ sr⁻¹

`q_1u = 8.629132 × 10⁻⁶ Υ_1u(T) exp[-ΔE_1u/(kT)] / [g_1 sqrt(T)]`

CHIANTI **11.0.2** atomic energies, collision strengths, radiative branches and equilibrium ion fractions are read directly. Burgess–Tully descaling follows the documented CHIANTI formulas, with observed energies used where available and theoretical energies as fallback. Ground-state collisional excitation and radiative cascades through the supplied model produce both resonance lines. At the contribution-function maximum, cascades change the doublet coefficient by about **0.43%**. No guessed Gaussian temperature-response curve is used.

Carbon abundance is 12 + log(n_C/n_H) = 8.43 and n_H/n_e = 1/1.2 in the emitting phase. The calculation assumes Maxwellian electrons, low-density excited-level populations, **coronal-approximation ionization equilibrium** and optically thin radiation. Using atomic files from v11 does **not** mean its advanced density/charge-transfer ionization model was implemented. Those effects and time-dependent carbon ionization can materially change UV brightness and formation temperatures. Resonance scattering and line opacity are omitted. This is an explicitly approximate emission synthesis, not validated absolute UV spectrophotometry.

Emissivity is calculated on the native T and n_e fields **before** 4×4 spatial averaging. Vertical cell-width weights conserve its volume integral. Averaging T first would create false intermediate-temperature radiation and is not done. The numerical data retain physical emissivity before any display enhancement.

## 2.5D geometry, projection and occultation

Two finite windows of the **same** local spicule sequence are placed at fixed prescribed stellar limb locations. This reuse does not imply independent 3D events or global magnetic connectivity. The invariant source direction receives a **6 Mm native finite depth**, with a 1 Mm cosine edge; that uncalculated depth is a display reconstruction assumption. The source x edges taper over 2 Mm, and the truncated 12–14 Mm top tapers smoothly. These define the finite domain window, not individual jet shapes.

The model is uniformly magnified **10×** in x, invariant depth and height. Its C IV coefficient is divided by ten in the enlarged coordinates, preserving a uniformly enlarged straight column's intensity. Gnomonic spherical wrapping has its actual volume Jacobian. Positive bilinear voxel projection integrates `j dV / pixel_area`, with 48 depth samples. The opaque photosphere blocks far-side material. No emissivity is placed inside the opaque star.

The projected surface is externally viewed and the star stays centered. True solar-scale heights and apparent full-frame velocities must be recovered by removing the disclosed magnification; enlarging the plasma does not make it a giant physical eruption. No large flare or coronal jet absent from the selected data is imposed. Larger UV structures can be model loops or heated boundaries rather than ejected mass.

## Clock and camera

All local source clocks advance at **20 physical seconds per playback second**, giving **640 physical seconds in 32 seconds**. Fixed different starting epochs are necessary because independent models are being assembled. Quiet tiles retain their fixed initial phase offsets. Optical background starts at 3852.4 s, the active photosphere at its original 26731.935 s epoch, and the GOL component at 974 s.

Physical timestamps bracket every sample. Linear interpolation of **radiation/emissivity**, rather than invented fluid velocities, produces display frames between recorded states. The C IV layers are thin and move rapidly: native-map alternate-output holdout errors are large. They are not concealed as proof of temporal convergence. All available selected outputs are retained, but subcadence motion is not independently validated. Source velocities must be read from MHD fields, not measured from a dissolving emissivity contour.

The atmospheric projector includes relative propagation delay `(R − z_observer)/c` in the assembled display geometry. This is a consistent display timing convention, not a globally exact light-travel model for magnified independent regions. Optical propagation uses the existing weak-field stellar-ray calculation. Gravitational redshift is negligible visually and is not exaggerated.

A separately prescribed 25.38-day rigid rotation advances only about 0.105° during the film. The camera has fixed inclination **67°** and unwrapped azimuth **−15° + 45° t**, with t in playback seconds. This gives exactly **four full 360° orbits in 32 s**, one every **8 s**. Its radius, projection, image center and apparent stellar radius remain constant. This accelerated viewpoint is independent of physical stellar rotation and plasma evolution. Frames are sampled at temporal midpoints `(i + 0.5)/30`; there is no duplicated endpoint, azimuth reset, source-time reset or repeated eight-second video segment. The evolving physical state does not return to its initial condition at 32 s. Each output is an instantaneous source/view sample; finite camera-shutter motion blur is not integrated. At the requested orbit speed, small structures can move several output pixels between frames.

Surface atlas transforms and the coherent magnetic region are anchored in stellar coordinates. The two UV windows also have stellar-coordinate centers established once, at initialization. All three layers use the same camera basis and physical rotation. Far-side photospheric structures become hidden; far-side atmospheric voxels are occulted by the opaque sphere when appropriate. No spot is copied to the rear and no UV sheet is turned toward the camera to disguise its finite-depth approximation. Consequently some views reveal that the source is 2.5D rather than a fully calculated 3D jet population.

## Deliberate grading, contrast and imaging

Every setting is fixed across frames; there is no frame-by-frame auto-exposure.

- Both photospheric sources: **10× spatial magnification**
- Quiet radiation contrast: `I_mean(μ) [I/I_mean(μ)]^3.25`, applied as a scalar to calculated spectral RGB
- Active radiation contrast: exponent **3.0** about a Gaussian **σ = 5 native pixels = 960 km** baseline from that same intensity frame; no generated small-scale structure
- Golden grade: let `r = Y / Y_initial_mean(μ)`, using the **fixed initial calculated limb law** including the unchanged redshift factor. Linearly interpolate RGB gain vectors **(1.38, 0.75, 0.075)** at r = **0.65**, **(1.28, 0.90, 0.13)** at r = **1.0**, and **(1.08, 1.00, 0.65)** at r = **1.9**, clamping beyond the endpoints. Multiply the calculated RGB by that vector and rescale each pixel to retain its exact **pre-grade CIE Y**. This continuous, fixed LUT gives dimmer lanes amber chromaticity and brighter calculated granules a warm-white highlight, without adding texture or changing the physical intensity array. μ here is the regular-pixel ray emission cosine, including the negligible weak-field radius correction; partial limb pixels use the same display reference. Physical spectra/RGB remain available upstream. The final chromaticities are intentionally changed; it would be false to say the graded film preserves the original RGB ratios or is true visible color
- C IV false-color linear RGB ratio **(1, 0.23, 0.035)**, normalized to unit Y; optical and UV physical intensities have **independent display exposures**, not a measured broadband channel ratio
- Fixed UV radial visibility gain `min[20, exp(max(0, projected_radius/R−1)/0.04)]`; this enhances calculated faint off-limb radiation and does not create plasma
- UV display intensity `0.8 asinh(I_enhanced/10⁴)/asinh(10)`, with I in erg cm⁻² s⁻¹ sr⁻¹; not a calibrated instrument response
- Normalized spectrum-weighted Airy core: first zero **1.2 px at 550 nm and 1080 output width**; Gaussian jitter **0.2 px**; normalized β = 2.5 Moffat wing with α = **105 px**, fraction **6%**
- The **additional** halo display gain of **5.5** is restricted by a smooth radial mask from 0.85–1.10 R so that it does not veil the disk. The normalized physical PSF still acts everywhere. This extra halo is a display effect, not emitted coronal plasma or a flux-conserving instrument response
- A common-channel smooth shoulder `out = 2 HDR / [1 + max(2 HDR)]` avoids clipping the gold disk to a uniform value. It changes displayed brightness nonlinearly while retaining each pixel's graded RGB ratio

The sphere occupies 81% of the portrait width, leaving room for the magnified limb structures. The black background contains no stars or decorative rays. The only labels are THE MAIN-SEQUENCE STAR, @sideratlas and NO.004, in Latin Modern fonts.

## Validation scope

See the machine-readable reports under `../validation/`. Tests cover source identity/timing/units, finite positive emission, conservative spatial averaging, depth quadrature, spectral-luminance preservation, photospheric output sampling, source temporal holdouts, PSF normalization and crop loss, and full video decode/color metadata. The inherited optical transfer and MURaM intensity/magnetic tests are retained separately.

These tests do **not** validate the independent models as one global star, establish native MHD grid convergence, recover unresolved penumbral filaments, calculate non-LTE chromospheric line transfer, establish equilibrium C IV photometry, or prove smooth subcadence UV advection. The scientific merit is real local MHD-driven dynamics with disclosed forward-model and display approximations.

## Measured checks for this release

The 71 selected GOL outputs cover **970–1630 s**, including irregular spacings of **0.05–10 s**. Near-duplicate states are physically close in the emissivity data; timestamps are not repaired to a guessed ten-second grid. The final 974–1614 s source-clock interval, including relative propagation offsets, lies inside the archive.

Float32-stored, spatially reduced UV emissivity retains the original volume integral to **3.4 × 10⁻⁹** relative accuracy. Whole-orbit quadrature was checked at eleven source/view times from 0 to 32 s. At 1080 output width, increasing invariant-depth quadrature from **48 to 96** samples changes the UV image by at most **2.282% L1**, with much smaller integrated-flux differences. This check is reused after verifying equality of every physical projection configuration field; presentation grading does not change the projected physical emission. A clean-package validation recomputes it if the historical report is absent. The delivered native 1080 render uses one regular photospheric sample plus analytic limb-pixel area and 8×8 limb sampling. Its one-versus-four-ray check is repeated with the new grading and reported in `../validation/scientific-checks.json`. These are post-processing tests, not native MHD convergence.

In contrast, omitting alternate GOL outputs changes native emissivity maps by median **116% L1**. Thin transition-region interfaces move across many source pixels between outputs. This is a substantive temporal-reconstruction limitation. All 71 selected native states are retained; interpolated frames must not be used for precision front-speed measurements. A separate projected holdout report records the effect at presentation resolution and does not hide it behind the bright photospheric disk.

The optical continuum and C IV layers share a chosen **presentation PSF** after assembly. The UV layer is not an image from a calibrated UV telescope with the optical wavelength's diffraction law. Color, co-registration and matched presentation resolution are deliberate composite-imaging choices. No calibrated optical-to-UV flux ratio is inferred from the final colors.

At the final 1080 projection resolution, four alternate-output C IV holdouts have **21.9–39.0% raw projected L1** differences. A declared 1.2 px Gaussian reduces these to **9.8–19.4%**, with integrated UV flux differences of **−3.8% to +7.5%**. This is a test of omitting native states, not a bound on ten-second interpolation error or validation of subcadence velocities. The final film retains every selected native output.

Refining native emissivity averaging from 4×4 to 2×2 at snapshot 325 changes the final raw projection by **1.90% L1** and the declared 1.2 px Gaussian-smoothed projection by **0.526% L1**, with integrated UV flux difference **0.00624%**. This supports the retained output averaging; it is not a new MHD mesh convergence result.


---

# Sources, rights and exact reuse

## Calculated plasma models

- **Gudiksen et al. (2011)**, *The stellar atmosphere simulation code Bifrost*, A&A 531, A154 — https://doi.org/10.1051/0004-6361/201116520
- **Hayek et al. (2010)**, radiative transfer with scattering — https://doi.org/10.1051/0004-6361/201014210
- **Carlsson et al. (2016)**, *A publicly available simulation of an enhanced network region of the Sun*, A&A 585, A4 — https://doi.org/10.1051/0004-6361/201527226; numeric source https://sdc.uio.no/vol/simulations/en024048_hion/atmos/
- **Cheung et al. (2019, online 2018)**, *A comprehensive three-dimensional radiative magnetohydrodynamic simulation of a solar flare*, Nature Astronomy 3, 160–166 — https://doi.org/10.1038/s41550-018-0629-3; numerical data https://purl.stanford.edu/dv883vb9686, **CC BY 3.0**. Only the photospheric numerical outgoing radiation and magnetic slices enter this movie; no source movie is used
- **Rempel (2017)**, *Extension of the MURaM Radiative MHD Code for Coronal Simulations* — https://doi.org/10.3847/1538-4357/834/1/10
- **Vögler et al. (2005)**, MURaM magnetoconvection — https://doi.org/10.1051/0004-6361:20041507
- **Martínez-Sykora, De Pontieu, Hansteen, Rouppe van der Voort, Carlsson & Pereira (2017)**, *On the generation of solar spicules and Alfvénic waves*, Science 356, 1269–1272 — https://doi.org/10.1126/science.aah5412
- **Martínez-Sykora, De Pontieu, Carlsson, Hansteen, Nóbrega-Siverio & Gudiksen (2017)**, *Two-dimensional Radiative Magnetohydrodynamic Simulations of Partial Ionization in the Chromosphere. II. Dynamics and Energetics of the Low Solar Atmosphere* — https://doi.org/10.3847/1538-4357/aa8866; full methods and public-data description https://arxiv.org/html/1708.06781; actual numeric source https://sdc.uio.no/vol/simulations/en096014_gol/atmos/
- **Martínez-Sykora et al. (2018)**, *Impact of Type II Spicules in the Corona: Simulations and Synthetic Observables* — https://arxiv.org/abs/1805.06475. Establishes the source model's multi-thermal behavior and discusses approximate equilibrium optically thin observables. It does not validate this project's finite-depth reconstruction

**Required MURaM attribution:** M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa, J. Martínez-Sykora, A. Sainz Dalda, M. L. DeRosa, A. Malanushenko, V. Hansteen, B. De Pontieu, M. Carlsson, B. Gudiksen and S. W. McIntosh. **Provided by the Stanford University Libraries.** Changes: radiative/color approximation, local-source sampling, spherical reconstruction, magnification, enhanced contrast, prescribed viewing geometry, optical imaging and video export. Original numerical data retain CC BY 3.0; this project's MIT code license does not replace it.

Bifrost data are publicly provided by ITA/University of Oslo and the IRIS modeling project for scientific reuse with citation of the code and experiment papers. No blanket MIT license is assigned to these data. The large source FITS files are retrieved by the scripts, not redistributed in the package. Exact variable, snapshot, byte range, header and payload hashes are retained in the source manifest.

## Radiation, color and imaging

- **CHIANTI atomic database**, freely available with requested scientific acknowledgement — https://www.chiantidatabase.org/ and https://www.chiantidatabase.org/chianti_download.html. Atomic files used are version **11.0.2**, from the public SolarSoft mirror https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/dbase/. The exact C IV source-file checksums are retained in the saved coefficient provenance; this collection supplies the calculated coefficient table rather than duplicating the raw atomic tables. CHIANTI is a collaborative project involving George Mason University, University of Michigan, University of Cambridge and NASA Goddard
- **Dere et al. (1997)**, initial CHIANTI database, A&AS 125, 149; **Dere et al. (2023)**, version 10.1 revised ionization/recombination rates — https://arxiv.org/abs/2305.15221
- **Dufresne et al. (2024)**, CHIANTI XVIII, version 11 advanced ionization models — https://doi.org/10.3847/1538-4357/ad6765. **Limit:** this film uses the density-independent `chianti.ioneq` coronal approximation, not the advanced charge-transfer/density model; the paper explains material differences this can cause
- **Asplund et al. (2009)**, solar abundance reference — https://doi.org/10.1146/annurev.astro.46.060407.145222. Adopted carbon abundance is 8.43
- **Coddington et al. (2023)**, TSIS-1 HSRS v2 — https://doi.org/10.1029/2022EA002637; source spectral file and hash in the optical manifest
- **CIE 1931 standard observer**, CC BY-SA 4.0 data — https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer; provider metadata retained
- H⁻ coefficient tables follow Geltman (1962), Stilley & Callaway (1970) and the **MIT-licensed Lightweaver** implementation — https://github.com/Goobley/Lightweaver/blob/master/Source/Background.cpp. The original MIT notice is in `licenses/Lightweaver-MIT.txt`
- **Latin Modern** fonts retain their GUST Font License in `fonts/`

The source papers establish methods and model behavior. They do not independently validate this visualization's atmosphere assembly, display grading, interpolated frames, absolute C IV photometry or new numerical convergence claims. No scientists were contacted and no observed solar movie, NASA animation or generative-image asset is used as a texture.


---

# Native active-region data and what it adds

[Cheung et al. (2019; published online 2018)](https://doi.org/10.1038/s41550-018-0629-3) released the numeric MURaM outputs through the [Stanford Digital Repository](https://purl.stanford.edu/dv883vb9686), **CC BY 3.0**. Required institutional attribution: Provided by the Stanford University Libraries. Credit: M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa, J. Martínez-Sykora, A. Sainz Dalda, M. L. DeRosa, A. Malanushenko, V. Hansteen, B. De Pontieu, M. Carlsson, B. Gudiksen, and S. W. McIntosh. This visualization modifies their data through normalization, temporal/spatial interpolation, spherical projection, documented contrast enhancement and optical imaging.

Only outgoing **photospheric** gray intensity and photospheric magnetic slices are used. The source study includes a corona and a flare, but those channels are not drawn. This is a visualization of the selected photospheric sequence, not a simulation of prominences or flare radiation.

`I_out_every200steps.tar` contains 300 frames with actual timestamps from 26729.535 to 30668.799 s. Its SHA-256 is `90c7afc88f93d4092f63152c5aaa3019e51e9b275db14edc2152e5f5bcb0482e`. The layout is the documented four-float header plus a 512 × 256 little-endian float32 image. The paper's 98.304 × 49.152 Mm domain gives 192 km sampling. Initial exploratory coordinates based on 384 km were discarded before the delivered preview and final render; only the confirmed 192 km scale is used.

Matching `tau_slice_0.100` fields have ten planes. The supplied authors' reader identifies planes 5,6,7 as horizontal x, **vertical z**, horizontal y, with multiplication by √(4π) giving gauss. Exact member names, timestamps, byte ranges and hashes are checked during retrieval. This prevents accidentally interpreting a horizontal field component as the vertical one.

For the chosen first-frame patch:

| Radiation class, defined only for validation | Median total field | Native mean radiation / quiet |
|---|---:|---:|
| Dark cores, I < 0.5 I_quiet | 2.97 kG | 0.306 |
| Magnetized intermediate-intensity structure | 1.29 kG | 0.724 |
| Bright points, I > 1.1 I_quiet and vertical field >500 G | 1.45 kG | 1.185 |

The intermediate-intensity class is not a perfect penumbral segmentation. Fine filaments are limited by the source grid. These classes validate the magnetic association of features; **they do not define rendering masks**. The rendered radiation is sampled directly from the full intensity map.

After 320 s, the weak-field patch's fixed-position intensity correlation drops to about 0.078, while the whole patch remains more correlated (about 0.70) because the large dark cores persist. This is the desired multiscale behavior: changing convection around more persistent magnetic structure. It is not an imposed lifetime or a synthetic noise model.

Holding out alternating actual snapshots and interpolating their neighbors gives mean RMS error about 1.29% of the fixed quiet reference, with the largest tested frame-average RMS about 1.83%. Individual local errors can be greater. Neutral contrast-one sampling recovers native outgoing gray intensity after the declared spectral factor with relative luminance error below 2 × 10⁻⁸ in the tested ray.

The gray radiation and native magnetic geometry are calculated. The full-star placement, patch blending, 10× magnification, display contrast, thermal color proxy thermal optical-band correction and active-region limb approximation are prescribed. The native vertical bright points should not be presented as a quantitative prediction of directional facular contrast near the limb. No unmodeled extended magnetic filaments are added to the photosphere.

![Native radiation and matching magnetic field](METHODS.md)

Figure: unenhanced native MURaM gray radiation and the matching vertical field, generated from the cited CC BY 3.0 data. The colored field map is documentation only and never appears in the Reel.


---

# Validation and remaining uncertainty

The numerical reports accompany the final pipeline. A successful numerical check does not make the assembled sphere a globally coupled MHD solution.

| Check | Result | Meaning |
|---|---:|---|
| Native GOL states retained | 71, 970–1630 s | Actual irregular timestamps bracket the film |
| C IV emissivity volume integral after stored averaging | maximum relative error 3.4 × 10⁻⁹ | Spatial averaging conserves calculated emission |
| UV depth samples, 48 vs 96 across eleven views/times | maximum 2.282% projected L1 | Oblique-view quadrature checked; original identical-geometry report retained with hash |
| UV emissivity averaging, 4×4 vs 2×2 | 1.90% raw L1, 0.526% after the declared Gaussian | Retained averaging is adequate at presentation resolution in the tested state |
| Photospheric sampling, 1 vs 4 regular rays per output pixel | 1.533% post-imaging RMS relative to disk mean | Native 1080 delivery sampling has been checked |
| Golden grade | maximum absolute linear-Y error 1.8 × 10⁻¹⁵ | Deliberate chromaticity change preserves pre-tone-map luminance |
| Cached source/radiation reproduction | completed successfully | Retrieved/cached inputs and calculated UV arrays match the manifest |
| Native UV alternate-output holdout | median 116% L1 | **Subcadence emissivity motion is not independently validated** |
| Four UV projected holdouts | 21.9–39.0% raw L1; 9.8–19.4% with the declared Gaussian | Fine UV changes cannot be described as precision continuous fluid trajectories |

The temporal holdouts intentionally omit real recorded outputs. They demonstrate the importance of retaining every selected state; they do not provide a proven bound on the remaining ten-second interpolation errors. The final animation linearly interpolates calculated emissivity and radiation, not a newly integrated velocity field.

The snapshot-325 source check finds real cool upflows above 2 Mm at roughly 30–45 km/s in the selected diagnostic class. The bright C IV phase is generally around 0.07–0.20 MK, with typical electron densities near 10⁹–10¹⁰ cm⁻³. Its velocity distribution includes downward and slowly moving material. The diagnostic finds **no cells meeting its additional >50 km/s, T<0.2 MK, density>10⁻¹³ kg m⁻³ cut at that state**. This is not a survey of the complete simulation, and no larger jet is fabricated to compensate. A bright UV strand must not automatically be called an ejected mass parcel.

Inherited optical tests include 0.035% RMS change from 17 to 33 transfer wavelengths, 1.58% RMS vertical formal-integration refinement, 2.28% RMS at an independent oblique angle, and 1.14% RMS for a 40-second optical midpoint interpolation. MURaM photospheric temporal holdouts average 1.29% RMS of the fixed quiet reference. These concern optical radiation; they must not be substituted for the larger UV errors.

The source's isolated ancillary Bifrost TEFF header anomaly is retained and flagged in the inherited source report. It never sets the emission temperature or exposure; local thermodynamic arrays supply the radiation. No native MHD mesh convergence or global energy transport test is claimed from these archive subsets.

The final export report verifies 960 frames, exactly 32 seconds, 1080 × 1920, H.264 High, Rec.709 limited YUV420p, complete decoding, nonidentical adjacent frames and the normalized PSF checks. Its recorded SHA256 identifies the delivered movie. Encoding follows one high-depth-to-H.264 conversion; source physics and display grading remain separate in the reproducible pipeline.

## Four-orbit revision

The 960 frame centers sample a continuous unwrapped azimuth at exactly 45°/s, or 1.5° per frame. The boundaries at 8, 16 and 24 seconds are ordinary consecutive frame pairs, not cuts or animation resets. Each eight-second interval contains exactly 240 frames. The analytic start/end camera transformation agrees to 1.5 × 10⁻¹⁵; actual quiet-surface and active-map samples return identically when physical time is held fixed. The source clocks and independent stellar rotation are not frozen in the film.

The new fixed chromatic LUT preserves pre-tone CIE Y to 1.8 × 10⁻¹⁵ in the tested full-resolution frame. Its amber lanes and warm-white highlights follow existing calculated intensity. One-versus-four-ray photospheric sampling changes the post-imaging result by 1.533% RMS relative to disk mean. This is presentation convergence, not source MHD convergence.

`validation/reused-inputs.json` records SHA256 hashes of the existing processed radiation inputs used here. `validation/scientific-checks.json` records the configuration-matched inherited UV quadrature and fresh camera/anchoring/grade tests. The final encoded-camera report checks the 240-frame orbit boundaries and four nonidentical evolving rotations. None of these checks removes the disclosed 2.5D geometry or native UV cadence limitation.

The original report behind the reused whole-orbit UV quadrature is retained as `validation/inherited-UV/whole-orbit-projection-original.json`; its SHA256 matches the fresh report's provenance record. Only its physical UV projection results are reused. Its unrelated historical palette test is not the golden film's grading or imaging model. A clean reproduction of the new validator recomputes the projection test directly.


---

# Third-party material

- `continuum-coefficients.json`: numerical tables adapted from the MIT-licensed Lightweaver `Source/Background.cpp`, copyright Christopher M. J. Osborne (2019). The MIT notice is retained in `licenses/Lightweaver-MIT.txt`. Scientific origins: Geltman (1962), Stilley & Callaway (1970), Mihalas (1978)
- `CIE_xyz_1931_2deg.csv`: unmodified CIE 1931 2° color-matching functions, DOI [10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b), International Commission on Illumination, 2019. CC BY-SA 4.0 as specified in the accompanying official metadata. Numerical integration and display conversion are this project's calculations
- Latin Modern font files: GUST Font License, retained under `fonts/`
- Bifrost atmosphere subsets: retrieved from the University of Oslo's public archive; citation/public-use statement described in `docs/SOURCES.md`. Raw source arrays are not included in the distributable package
- TSIS-1 HSRS v2: LASP reference spectrum, Coddington et al. (2023); retrieved by `fetch_auxiliary.py`; the exact reference file is included for reproducibility. Institutional provenance, official URL and SHA-256 are retained in `source-manifest.json`

The MIT license for this project's original code does not relicense any third-party scientific data.

- Active-region MURaM data: Cheung et al. (2019; online 2018), Stanford Digital Repository, [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/), [source](https://purl.stanford.edu/dv883vb9686). Attribution, author credit, and modifications are detailed in `docs/ACTIVE-REGION-ATTRIBUTION.md`. Numeric gray radiation and magnetic fields are downloaded from this source; the documentation figure is derived from those data. The source authors’ movie footage is not used

CHIANTI atomic files and the public GOL model retain their original provider rights; see docs/SOURCES.md.

