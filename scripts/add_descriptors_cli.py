#!/usr/bin/env python
"""CLI tool to enrich a CSV with molecular descriptors."""
import sys
from pathlib import Path

# Ensure src is on path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.clean.add_descriptors import add_descriptors_to_csv


def main():
    if len(sys.argv) < 3:
        print("Usage: python scripts/add_descriptors_cli.py <input.csv> <output.csv>")
        sys.exit(1)
    
    in_csv = sys.argv[1]
    out_csv = sys.argv[2]
    
    print(f"Enriching {in_csv} with molecular descriptors...")
    add_descriptors_to_csv(in_csv, out_csv)
    print(f"✓ Descriptors added. Saved to {out_csv}")


if __name__ == "__main__":
    main()
