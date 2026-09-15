#!/usr/bin/env python3
"""Run both analysis versions on the two small demo series."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tifffile

from src.trafficking_core import (
    AnalysisError,
    positive_float,
    run_analysis,
    threshold_value,
    write_json,
)


ROOT = Path(__file__).resolve().parent
DEMO_DIRECTORIES = {
    "control": ROOT / "demo_data" / "control",
    "experiment": ROOT / "demo_data" / "experiment",
}


def validate_output_root(output_root: Path, input_directories: list[Path]) -> None:
    output = output_root.expanduser().resolve()
    for input_directory in input_directories:
        source = input_directory.expanduser().resolve()
        if output == source or source in output.parents:
            raise AnalysisError(
                f"demo output must be separate from every input directory: {output}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output-directory",
        type=Path,
        default=ROOT / "demo_output",
    )
    parser.add_argument("-t", "--threshold", type=threshold_value, default=0)
    parser.add_argument("-roi", "--roi")
    parser.add_argument("--frame-interval-seconds", type=positive_float)
    args = parser.parse_args()

    output_root = args.output_directory.expanduser().resolve()
    try:
        validate_output_root(output_root, list(DEMO_DIRECTORIES.values()))
    except AnalysisError as exc:
        parser.error(str(exc))
    if output_root.exists() and (not output_root.is_dir() or any(output_root.iterdir())):
        parser.error(f"output directory is not empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    summary: list[dict[str, object]] = []
    try:
        for dataset, input_directory in DEMO_DIRECTORIES.items():
            results: dict[str, dict[str, object]] = {}
            for version in ("v1", "v2"):
                output_directory = output_root / f"{dataset}_{version}"
                results[version] = run_analysis(
                    version=version,
                    input_directory=input_directory,
                    output_directory=output_directory,
                    threshold=args.threshold,
                    roi_text=args.roi,
                    parameters_confirmed=False,
                    dataset_label=dataset,
                    frame_interval_seconds=args.frame_interval_seconds,
                )

            v1_mask = tifffile.imread(
                output_root / f"{dataset}_v1" / "absolute_change_union.tif"
            )
            v2_mask = tifffile.imread(
                output_root / f"{dataset}_v2" / "absolute_change_union.tif"
            )
            mismatch = int(np.count_nonzero(v1_mask != v2_mask))
            signatures_match = (
                results["v1"]["comparison_signature"]
                == results["v2"]["comparison_signature"]
            )
            if mismatch or not signatures_match:
                raise AnalysisError(f"{dataset}: v1 and v2 absolute results do not match")
            summary.append(
                {
                    "dataset": dataset,
                    "frames": results["v1"]["frame_count"],
                    "roi_pixels": results["v1"]["roi_pixels"],
                    "absolute_change_pixels": results["v1"]["absolute_change_pixels"],
                    "absolute_change_percent": results["v1"]["absolute_change_percent"],
                    "increase_pixels": results["v2"]["increase_pixels"],
                    "increase_percent": results["v2"]["increase_percent"],
                    "decrease_pixels": results["v2"]["decrease_pixels"],
                    "decrease_percent": results["v2"]["decrease_percent"],
                    "overlap_pixels": results["v2"]["overlap_pixels"],
                    "overlap_percent": results["v2"]["overlap_percent"],
                    "v1_processing_seconds": results["v1"]["processing_seconds"],
                    "v2_processing_seconds": results["v2"]["processing_seconds"],
                    "v1_v2_mismatch_pixels": mismatch,
                    "v2_or_identity": results["v2"][
                        "absolute_equals_increase_or_decrease"
                    ],
                }
            )
    except AnalysisError as exc:
        parser.exit(2, f"Error: {exc}\n")

    write_json(output_root / "summary.json", summary)
    print(
        "dataset     frames   ROI   abs    abs %   inc   inc %   dec   dec % "
        " overlap  overlap %    v1 s    v2 s  mismatch  OR"
    )
    for row in summary:
        print(
            f"{row['dataset']:<11} {row['frames']:>6} "
            f"{row['roi_pixels']:>5} "
            f"{row['absolute_change_pixels']:>5} "
            f"{row['absolute_change_percent']:>8.4f} "
            f"{row['increase_pixels']:>5} "
            f"{row['increase_percent']:>7.4f} "
            f"{row['decrease_pixels']:>5} "
            f"{row['decrease_percent']:>7.4f} "
            f"{row['overlap_pixels']:>7} "
            f"{row['overlap_percent']:>10.4f} "
            f"{row['v1_processing_seconds']:>7.4f} "
            f"{row['v2_processing_seconds']:>7.4f} "
            f"{row['v1_v2_mismatch_pixels']:>9} "
            f"{str(row['v2_or_identity']):>3}"
        )
    print(f"Outputs: {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
