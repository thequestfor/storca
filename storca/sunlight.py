"""Acquisition and provenance for the default terrestrial reference spectrum."""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from urllib.request import urlopen


# pvlib distributes this machine-readable copy from NREL's ASTM G173-03 Excel
# table.  The retained output records both the immediate and original source.
PVLIB_AM15G_URL = "https://raw.githubusercontent.com/pvlib/pvlib-python/main/pvlib/data/ASTMG173.csv"
NREL_AM15G_URL = "https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls"


def fetch_default_am15g(destination: Path) -> Path:
    """Fetch AM1.5 global-tilt irradiance and normalize it for STORCA."""
    with urlopen(PVLIB_AM15G_URL, timeout=30) as response:
        lines = response.read().decode("utf-8").splitlines()
    # The upstream copy includes a one-line title before its CSV header.
    rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:]))))
    if not rows:
        raise RuntimeError("Downloaded AM1.5 spectrum is empty")
    fields = rows[0].keys()
    wavelength = next((field for field in fields if "Wvlgth" in field or "wavelength" in field.lower()), None)
    global_tilt = next((field for field in fields if "Global" in field or "global" in field.lower()), None)
    if not wavelength or not global_tilt:
        raise RuntimeError("Downloaded AM1.5 spectrum lacks wavelength/global-tilt columns")
    output_rows = [{"wavelength_nm": row[wavelength], "irradiance_W_m2_nm": row[global_tilt]} for row in rows]
    destination = Path(destination); destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "irradiance_W_m2_nm"])
        writer.writeheader(); writer.writerows(output_rows)
    destination.with_suffix(destination.suffix + ".metadata.json").write_text(json.dumps({
        "source": NREL_AM15G_URL, "machine_readable_copy": PVLIB_AM15G_URL,
        "spectrum": "ASTM G173-03 AM1.5 global hemispherical / 37-degree tilt",
    }, indent=2) + "\n")
    return destination
