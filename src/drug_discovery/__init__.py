"""Computational drug discovery pipeline for acetylcholinesterase inhibitors."""

from typing import Any

__all__ = ["TARGET_CHEMBL_ID", "TARGET_NAME", "run_part1"]


def __getattr__(name: str) -> Any:
    if name in {"TARGET_CHEMBL_ID", "TARGET_NAME"}:
        from . import config

        return getattr(config, name)
    if name == "run_part1":
        from .part1 import run_part1

        return run_part1
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
