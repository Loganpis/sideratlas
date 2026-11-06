# Attribution and rights

The MIT [LICENSE](LICENSE) covers original software, documentation and author-generated numerical caches, not all scientific data, fonts or media. Rendered media retain SiderAtlas's rights and applicable third-party attribution; no blanket MIT media license is asserted. No author or institution endorsement is implied.

## SiderAtlas No.004 code

Generic spherical sampling, grading, optical imaging and export algorithms are adapted from [SiderAtlas No.004](https://github.com/Loganpis/sideratlas/tree/main/episodes/004-main-sequence-star). Its original MIT copyright and license are retained in [licenses/No004-MIT.txt](licenses/No004-MIT.txt). Changes are import/path portability, selected-scene assembly and stellar-merger use; solar thermodynamics are not transplanted as giant-star physics.

## Latin Modern fonts

GUST e-foundry / Latin Modern Roman. The bundled font files are unmodified and retain the complete [GUST Font License](fonts/GUST-FONT-LICENSE.txt). This notice is independent of the code license.

## CIE observer table

International Commission on Illumination, **CIE 2019**, [Colour-matching functions of CIE 1931 standard colorimetric observer](https://doi.org/10.25039/CIE.DS.xvudnb9b). The numerical table is unmodified and remains **[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)**, not MIT. Its original [metadata and rights statement](data/spectral/CIE_xyz_1931_2deg.csv_metadata.json) are supplied.

## MURaM-derived illustrative activity

[Cheung et al. dataset](https://purl.stanford.edu/dv883vb9686), **[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)**. **Provided by the Stanford University Libraries.**

Credit: M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa, J. Martínez-Sykora, A. Sainz Dalda, M. L. DeRosa, A. Malanushenko, V. Hansteen, B. De Pontieu, M. Carlsson, B. Gudiksen and S. W. McIntosh.

The released derived intensity-ratio template selects every tenth frame, normalizes to a quiet reference, clips ratios, rescales to an illustrative giant surface, interpolates in time, projects spherically and applies the documented spectral/display treatment. It is not a calculation of the donor's magnetic field. No raw magnetic cubes are redistributed.

## Bifrost-derived morphology

Source: the [Oslo enhanced-network Bifrost archive](https://sdc.uio.no/vol/simulations/en024048_hion/atmos/), Gudiksen et al. (2011) and Carlsson et al. (2016). The prior provenance records public-use/citation terms, not a confirmed explicit redistribution license. Therefore the new release does **not** redistribute morphology.npz or mesoscales.npz and does not assign them an invented open license.

Their hashes, original source-map provenance and derivation recipe are provided for users with the applicable inputs and rights. Required transformations include intensity normalization, log-contrast extraction, quantization, coarse filtering, temporal interpolation, rotated overlapping spherical charts, giant-scale rescaling and disclosed enhancement. See [reproduction limits](METHODS.md) and [surface attribution](data/quality-surface/ATTRIBUTION.md).

## MESA and author-generated numerical outputs

MESA release 26.4.1 supplied previously calculated one-dimensional stellar structure. [MESA's source](https://github.com/MESAHub/mesa) is LGPL licensed; this repository does not redistribute MESA source, EOS tables, SDK or compiled libraries. The supplied profile/global parameters are model outputs. Their provenance does not license MESA itself under MIT or claim a new stellar-evolution calculation was reproduced during publication.

The permitted render-input asset also contains author-generated trajectory, material-support, force-interpolation, launch-reference and initialization arrays. These author-generated numerical outputs are supplied under the MIT license with attribution to SiderAtlas; third-party exclusions and notices above remain applicable. The two external Bifrost-derived maps are explicitly outside that asset.

NumPy, SciPy, Numba, Pillow and FFmpeg are separately installed dependencies; their libraries/binaries are not bundled. The selected source does not include StarSmasher, Athena++, SEAGen or MacLeod solver code. Their earlier investigation is not grounds to imply their use in the final trajectory integration.
