#!/usr/bin/env python3
"""
build_all.py -- regenerate every per-script .docx in technical_for_me/.

Run this any time you change one of the chapter builders to refresh the
entire set in one command:

    cd mbXPro/technical_for_me
    python3 build_all.py
"""

from pathlib import Path
import sys

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from _build_technical import main as part1
from _chapters_part2 import main as part2
from _chapters_part3 import main as part3


if __name__ == "__main__":
    print("=" * 70)
    print(" Building technical_for_me/ -- internal per-script reference docs")
    print("=" * 70)
    part1()
    part2()
    part3()
    print()
    print("All 20 documents generated. Output:")
    print(f"  {HERE.resolve()}")
