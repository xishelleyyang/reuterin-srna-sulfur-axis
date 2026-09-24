#!/usr/bin/env python
"""Fetch every raw input this pipeline needs, straight from GEO / eLife / NCBI RefSeq.

Nothing under data/ is redistributed with this repository: config/input_manifest.json
lists, for every file, its public source, its exact download URL, and the SHA-256 this
pipeline's reported numbers were last verified against. Each file is downloaded to a
temporary path, hashed, and only moved into place on a match; a mismatch (or a network
failure) is a hard error, never a silent skip or a "close enough" substitution.

Usage:
    python download_inputs.py --manifest config/input_manifest.json --dest data
    python download_inputs.py --dest data --only genomes/S.Typhimurium_LT2.fna.gz
    python download_inputs.py --dest data --skip-existing   # default; re-verifies hash of files already present
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

CHUNK = 1 << 20  # 1 MiB


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, timeout: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "reuterin-srna-sulfur-axis/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)


def fetch_one(entry: dict, dest_root: Path, skip_existing: bool) -> str:
    final_path = dest_root / entry["relpath"]
    if skip_existing and final_path.exists():
        have = sha256_of(final_path)
        if have == entry["sha256"]:
            return "already-present (hash verified)"
        print(f"  existing file hash mismatch, re-downloading: {final_path}", file=sys.stderr)

    tmp_path = final_path.with_suffix(final_path.suffix + ".part")
    download(entry["url"], tmp_path)
    got = sha256_of(tmp_path)
    got_bytes = tmp_path.stat().st_size
    if got != entry["sha256"]:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"SHA-256 mismatch for {entry['relpath']}: expected {entry['sha256']}, got {got} "
            f"({got_bytes} bytes vs {entry['bytes']} expected). Source may have changed "
            f"upstream ({entry['source']}); do not proceed without manually re-checking "
            f"this file's content and updating the manifest deliberately.")
    tmp_path.rename(final_path)
    return f"downloaded and hash-verified ({got_bytes} bytes)"


def main(manifest_path: str, dest: str, only: list[str] | None, skip_existing: bool) -> None:
    manifest = json.loads(Path(manifest_path).read_text())
    entries = manifest["files"]
    if only:
        entries = [e for e in entries if e["relpath"] in only]
        missing = set(only) - {e["relpath"] for e in entries}
        if missing:
            raise SystemExit(f"--only referenced unknown relpath(s): {sorted(missing)}")

    dest_root = Path(dest)
    dest_root.mkdir(parents=True, exist_ok=True)
    print(f"Fetching {len(entries)} file(s) into {dest_root}/ ...")
    failures = []
    for i, entry in enumerate(entries, 1):
        label = f"[{i}/{len(entries)}] {entry['relpath']}"
        try:
            status = fetch_one(entry, dest_root, skip_existing)
            print(f"{label}: {status}")
        except Exception as exc:  # noqa: BLE001 -- summarized and re-raised at the end
            print(f"{label}: FAILED -- {exc}", file=sys.stderr)
            failures.append((entry["relpath"], str(exc)))

    if failures:
        print(f"\n{len(failures)} of {len(entries)} file(s) failed:", file=sys.stderr)
        for relpath, msg in failures:
            print(f"  - {relpath}: {msg}", file=sys.stderr)
        raise SystemExit(1)
    print(f"\nAll {len(entries)} input files verified in {dest_root}/.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", default="config/input_manifest.json")
    p.add_argument("--dest", default="data")
    p.add_argument("--only", nargs="*", default=None, help="fetch only these relpath(s) from the manifest")
    p.add_argument("--force", action="store_true", help="re-download even if a hash-matching file already exists")
    a = p.parse_args()
    main(a.manifest, a.dest, a.only, skip_existing=not a.force)
