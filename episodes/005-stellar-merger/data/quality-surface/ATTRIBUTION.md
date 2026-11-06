# Derived morphology assets

These arrays are a statistical rendering surrogate, **not giant-star RHD or
calculated donor magnetic fields**. Source temperature, gravity and spectra
are discarded. Mean stellar properties are from the local MESA donor.

## Granulation

`morphology.npz` and `mesoscales.npz` derive from SiderAtlas No.004's calculated
Bifrost continuum intensity sequence. The upstream atmosphere is the University
of Oslo Bifrost public archive, Gudiksen et al. (2011) and Carlsson et al. (2016),
DOI https://doi.org/10.1051/0004-6361/201527226. The archive's public-use/citation
terms are retained in the accompanying No.004 research provenance. This is not
an invented MIT license for the scientific data.

Changes: angle selection, log-intensity mean removal, per-map normalization,
int16 quantization, coarse filtering, temporal interpolation, rotated overlapping
spherical charts, giant-scale rescaling, magnification and contrast enhancement.
Source file hashes and exact model/display factors are in `provenance.json`.

## Illustrative activity

`illustrative-active-region.npz` derives from Cheung et al. (2019; online 2018),
Stanford Digital Repository, https://purl.stanford.edu/dv883vb9686,
**CC BY 3.0**, https://creativecommons.org/licenses/by/3.0/.

**Provided by the Stanford University Libraries**

Credit: M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa,
J. Martínez-Sykora, A. Sainz Dalda, M. L. DeRosa, A. Malanushenko, V. Hansteen,
B. De Pontieu, M. Carlsson, B. Gudiksen, and S. W. McIntosh.

Changes: every tenth outgoing-intensity frame, normalization to a quiet
reference, ratio clipping, rescaling to an illustrative giant spot temperature,
angular placement/blending, temporal interpolation and fading, spectral
recalculation and display contrast. The solar magnetic arrays are not used.
No author endorsement is implied.

## Other source material

Generic No.004 rendering code carries its retained MIT notice in
`third_party/sideratlas-no004-rendering/`. The CIE 1931 observer response is
International Commission on Illumination (2019), DOI 10.25039/CIE.DS.xvudnb9b,
CC BY-SA 4.0; the numerical CSV is unmodified. Latin Modern fonts carry the
included GUST Font License. These notices are independent of any license for
original SiderAtlas code.
