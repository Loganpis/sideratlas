# Reproduction and release checks

The exact accepted movies are authoritative artifacts. Consolidation preserves their bytes; it does not claim six newly rerendered films or new validation of their astrophysical approximations. Per-episode provenance pins the original selected source revision and final SHA-256.

| Episode | Practical reproduction scope |
|---|---|
| 001 | Source/configuration and numerical checks included; a full new ray trace and C++ build are required for the film; full regeneration not repeated here |
| 002 | Saved photon scene renders without a field solve; 63 software/numerical checks pass; optional GPL native source is a separate release asset |
| 003 | Saved modes and cellular model; original calculated radiation caches supplied separately; raw provider tables not needed for final rendering; bounded frame parity and statistical checks tested; no fresh GYRE solve |
| 004 | Source/configuration preserved; large external MHD inputs must be reconstructed; full transfer and native film not repeated here |
| 005 | Cached gas workflow tested originally; two Bifrost-derived surface inputs remain unredistributed; public clone alone cannot render the complete scene |
| 006 | Frozen radial/photon histories are separate assets; external MIST/rate/No.004 inputs required; rebuilt Metal projector passed a bounded CPU comparison; full movie not rerendered here |

Consolidation measurements: [validation.json](validation.json).

## What was checked

All six movies were matched against accepted publication manifests, fully decoded, and checked for 1080 × 1920, 30 fps, 960 frames and 32 seconds. Thumbnails are extracted from those exact exports. Source selection uses release decisions and numerical manifests, not file dates. The latest GIANT CELLS export includes the accepted No.006 radius repair; the nebular calculation was not restarted.

The shared encoder’s bytes and generated command match the original No.004 implementation. Its copies in 004/005/006 are relative symlinks. Nine No.004 modules are unchanged and now provide No.006’s solar adapter; only its sibling-directory locator changed. No.003 applicable source remains a fingerprint-verified vendor subset in No.006. Unmodified fonts are content-deduplicated with original notices. Scientific routines are not consolidated across incompatible models or licenses.

The No.003 144 × 256 source-frame comparison was pixel-identical at global frame 181. No.006’s rebuilt projector passed the original arithmetic tolerances in a bounded 27-second-state comparison; its original native-resolution qualification remains separately archived.

No.002’s scene metadata had an unused machine-local path prefix removed. Its numerical arrays are exactly equal to the accepted original scene; the old and portable archive checksums are both recorded. Archived validation reports describe original runs; current consolidation checks are separately identified. There is no claim that a passing software test validates agreement with observations.

No movies, numerical-history ZIPs, native executables, development environments, intermediate frames or raw external solar cubes are tracked. Provider terms remain distinct from code licenses. External-source availability is not guaranteed. Source downloads are checksum checked and cannot silently substitute different data.

## Clean checkout

Clone with Git symlink support (tested target: macOS/Python 3.12). Use a separate environment for each episode; their NumPy and optional native dependencies differ. FFmpeg/FFprobe with libx264 and zscale are external dependencies.

```sh
git clone https://github.com/Loganpis/sideratlas.git
cd sideratlas
python3 tools/check_release.py
```

A fresh local Git clone passed the complete static audit and all 63 pulsar tests. Its external white-dwarf input restoration also passed individual-file hashes. These checks used the existing Python 3.12 environment; six fresh dependency environments and six full native rerenders were not performed.

Each episode README provides its own rendering procedure. The input restorer never launches a physical solver. It refuses changed existing inputs, validates archive and individual-file hashes, and extracts only listed members. No.006 also refuses automatic replacement of a missing selected ionization history.

Original donor/SPH experiments and rejected visual candidates remain in local development storage, outside this repository. The individual repositories are superseded after verification of selected final inputs in this unified repository. Only the accepted numerical inputs were migrated; no repository backups or historical archives were created for deletion. Normal commits do not purge older Git blobs; original source tags remain historical freezes, and current reproduction instructions use main.

## Permanent-cleanup qualification

The exact selected pulsar field, original calculated white-dwarf integration caches and permitted stellar-merger input asset are hosted by the corresponding unified releases. No active input URL depends on an individual SiderAtlas repository. Original provider attributions and unresolved external surface-data limits remain explicit. Scientific callable bodies are preserved; the white-dwarf fixed-atmosphere data interface is checked against original pixels without access to raw model tables.
