"""Selected SiderAtlas No.007 implementation. kpc, Myr and solar masses."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
CONFIG=json.loads((ROOT/'config.json').read_text())
