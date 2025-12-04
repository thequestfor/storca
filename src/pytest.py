import requests
import json
import re

APPEARANCE_KEYWORDS = [
    "color", "colour", "appearance", "odor", "odour",
    "liquid", "solid", "gas", "fuming", "oily", "crystalline",
    "powder", "solution", "amorphous", "viscous"
]


def get_sids(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sids/JSON"
    r = requests.get(url)
    data = r.json()
    return data.get("InformationList", {}).get("Information", [{}])[0].get("SID", [])


def get_sid_record(sid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/JSON"
    r = requests.get(url)
    return r.json()


def scan_json_for_strings(obj):
    """Yield every string found in the JSON recursively."""
    if isinstance(obj, dict):
        for v in obj.values():
            yield from scan_json_for_strings(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from scan_json_for_strings(item)
    elif isinstance(obj, str):
        yield obj


def scan_json_for_appearance(obj):
    """Return any strings that contain physical-appearance keywords."""
    found = []
    for s in scan_json_for_strings(obj):
        low = s.lower()
        if any(k in low for k in APPEARANCE_KEYWORDS):
            found.append(s)
    return found


def debug_pubchem_appearance(cid):
    print(f"=== DEBUG: CID {cid} ===")

    # 1. Get all SIDs for the compound
    sids = get_sids(cid)
    print(f"Found {len(sids)} SIDs: {sids[:10]}{'...' if len(sids)>10 else ''}")

    if not sids:
        print("No substance records available.")
        return

    # 2. Analyze each SID
    for sid in sids[:10]:  # limit to first 10 to avoid spam
        print(f"\n--- Checking SID {sid} ---")
        record = get_sid_record(sid)

        # Pretty-print top-level keys
        print("Top-level keys:", list(record.keys()))

        # 3. Find all appearance strings
        hits = scan_json_for_appearance(record)

        if hits:
            print(">>> Appearance-related strings found:")
            for h in hits:
                print("   ", h)
        else:
            print("No appearance text found in this SID.")

        print("-"*50)


# Run example: POCl3 = CID 24813
if __name__ == "__main__":
    debug_pubchem_appearance(24813)

