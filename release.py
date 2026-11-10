"""Build and check a portable, deterministic source-and-paper release ZIP."""

import hashlib
from pathlib import Path
import zipfile

from fragment import ROOT


VERSION = "0.2.0"


def release_files():
    names = ["README.md", "结课论文.md", "LICENSE", ".gitignore", ".gitattributes",
             "requirements-lock.txt", "pytest.ini", "study.json", "extension.json",
             "fragment.py", "experiment.py", "analyze.py", "compare_endpoints.py",
             "extension.py", "verify_project.py", "release.py", "data/study_manifest.csv",
             "analysis/global_local_discordance.csv", "analysis/baseline_reproduction.json"]
    for directory in ("tests", "examples", "data/raw", "results", "analysis/extension", "docs", ".github/workflows"):
        names.extend(p.relative_to(ROOT).as_posix() for p in (ROOT / directory).rglob("*")
                     if p.is_file() and "__pycache__" not in p.parts)
    return sorted(set(names))


def build():
    from verify_project import verify

    names = release_files()
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}
    (ROOT / "SHA256SUMS.txt").write_text("".join(f"{value}  {name}\n" for name, value in hashes.items()),
                                        encoding="utf-8", newline="\n")
    print(verify())
    destination = ROOT / f"蛋白质片段定位_源码与论文_v{VERSION}.zip"
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in [*names, "SHA256SUMS.txt"]:
            entry = zipfile.ZipInfo(name, date_time=(2026, 11, 10, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, (ROOT / name).read_bytes(), compresslevel=9)
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:
            raise ValueError("Corrupt release ZIP")
        for name, digest in hashes.items():
            if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError(f"Archive checksum mismatch: {name}")
    print(f"Built {destination.name}: {destination.stat().st_size} bytes, {len(names)+1} files")


if __name__ == "__main__":
    build()
