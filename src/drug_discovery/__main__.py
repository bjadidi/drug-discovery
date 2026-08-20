"""CLI: python -m drug_discovery part1 [--skip-fetch]"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Computational drug discovery pipeline"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    part1 = sub.add_parser("part1", help="Download and curate ChEMBL bioactivity data")
    part1.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Reuse data/raw CSV if it already exists",
    )
    part1.add_argument(
        "--no-search",
        action="store_true",
        help="Skip the target-search step and fetch CHEMBL220 directly",
    )

    for name, hint in [
        ("part2", "Exploratory data analysis (not implemented yet)"),
        ("part3", "Descriptor calculation (not implemented yet)"),
        ("part4", "Regression models (not implemented yet)"),
        ("part5", "Algorithm comparison (not implemented yet)"),
    ]:
        sub.add_parser(name, help=hint)

    args = parser.parse_args(argv)

    if args.command == "part1":
        from .part1 import run_part1

        run_part1(skip_fetch=args.skip_fetch, search=not args.no_search)
        return 0

    print(
        f"{args.command} is scaffolded in notebooks/ but not implemented yet. "
        "Finish Part 1 first.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
