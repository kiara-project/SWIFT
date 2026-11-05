"""Build the stand-alone KIARA dust checks for several bin counts without
Grackle or GSL, with address and undefined-behaviour sanitizers.

Run from anywhere: python3 tests/testKIARADustGrid.py
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="kiara-dust-grid-") as scratch:
    for source in sorted((root / "tests").glob("testKIARADust*.c")):
        for count in (2, 3, 6, 12, 41):
            binary = Path(scratch) / f"{source.stem}-{count}"
            subprocess.run(
                shlex.split(os.environ.get("CC", "cc"))
                + ["-std=c11", "-Wall", "-Wextra", "-Werror",
                   "-fsanitize=address,undefined", f"-DKIARA_DUST_N_BINS={count}",
                   "-I", str(root / "src"), str(source), "-lm", "-o", str(binary)],
                check=True,
            )
            subprocess.run([str(binary)], check=True)

