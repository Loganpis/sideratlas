# Shared components

`video_export.py` is the byte-identical RGB48 → Rec.709/H.264 encoder used by episodes 004, 005 and 006. Relative source symlinks preserve imports and source-byte hashes. Its MIT notice is in LICENSE. Encoder command generation is checked against the selected original source.

`fonts/` contains unmodified, content-deduplicated font files and original GUST/LPPL/OFL notices. Episode-relative symlinks preserve their expected paths. These asset licenses remain distinct from software licensing.

Episode 006 imports the canonical episode 004 solar renderer through a path-only adapter change. Its nine reused scientific modules match the original selected bytes. A small frozen No.003 source subset remains under episode 006/vendor because its AST loader verifies exact historical callable bodies. No GPL pulsar code is shared into the MIT episodes.
