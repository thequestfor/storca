"""Preserve exact historical transfer sources before implementing a continuation.

This is a zero-calculation archive utility, not a way to authorize execution.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def binding(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def archive_sources(campaign, destination):
    from storca.rate_transfer_execution import reopen_execution
    campaign, destination = Path(campaign).resolve(), Path(destination).resolve()
    if destination.exists():
        manifest = json.loads((destination / "source-capsule.json").read_text())
        if manifest["execution_contract"] != binding(campaign / "execution-contract.json"):
            raise ValueError("Source capsule belongs to a different contract")
        for item in manifest["files"]:
            if binding(destination / item["relative_path"])["sha256"] != item["sha256"]:
                raise ValueError("Archived source bytes changed")
        return manifest
    contract = reopen_execution(campaign)
    destination.mkdir(parents=True)
    paths = set()
    for folder in ("storca", "src", "tests"):
        paths.update(p.resolve() for p in (ROOT / folder).rglob("*.py")
                     if "__pycache__" not in p.parts)
    paths.add(ROOT / "pyproject.toml")
    required = {Path(item["path"]).resolve(): item["sha256"] for item in contract["implementation"]}
    paths.update(required)
    files = []
    for source in sorted(paths):
        relative = source.relative_to(ROOT)
        data = source.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        if source in required and digest != required[source]:
            raise ValueError("Original bound implementation source differs")
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
        files.append({"original_path": str(source), "relative_path": str(relative),
                      "sha256": digest, "originally_bound": source in required})
    manifest = {"schema_version": 1, "kind": "rate_transfer_historical_source_capsule",
                "execution_contract": binding(campaign / "execution-contract.json"),
                "files": files, "native_jobs_launched": 0}
    with (destination / "source-capsule.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    result = archive_sources(args.campaign, args.destination)
    print(json.dumps({"status": "historical_sources_archived", "files": len(result["files"]),
                      "originally_bound": sum(x["originally_bound"] for x in result["files"]),
                      "native_jobs_launched": 0}, indent=2))
