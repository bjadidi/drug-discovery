"""CLI: python -m drug_discovery {part1,part2,part3,part4,part5,all}"""

from __future__ import annotations

import argparse


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
        ("part2", "Lipinski descriptors, pIC50, Mann-Whitney tests, EDA figures"),
        ("part3", "Morgan fingerprints and Bemis-Murcko scaffolds"),
        ("part4", "Random forest: random split vs. scaffold split"),
        ("part5", "Compare algorithms on a scaffold split"),
    ]:
        sub.add_parser(name, help=hint)
    run_all = sub.add_parser("all", help="Run parts 1-5 in order")
    run_all.add_argument("--skip-fetch", action="store_true")

    args = parser.parse_args(argv)

    from . import parts

    if args.command in ("part1", "all"):
        from .part1 import run_part1

        run_part1(skip_fetch=args.skip_fetch, search=not getattr(args, "no_search", False))
    if args.command in ("part2", "all"):
        parts.run_part2()
    if args.command in ("part3", "all"):
        parts.run_part3()
    if args.command in ("part4", "all"):
        parts.run_part4()
    if args.command in ("part5", "all"):
        parts.run_part5()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
