#!/usr/bin/env python3
"""Shared implementation for the revised trafficking analyses."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
import tifffile
from PIL import Image, __version__ as pillow_version


FRAME_PATTERN = r"^.+_T(?P<time>\d+)_RGB_TRITC\.(?:tif|tiff)$"
GAUSSIAN_WEIGHTS = (1, 4, 6, 4, 1)
ALGORITHM_VERSION = "1.0.0"
ALGORITHM = "full-frame blur -> ROI -> red channel -> difference -> threshold -> temporal OR"


class AnalysisError(RuntimeError):
    """A clear, expected analysis failure."""


@dataclass(frozen=True)
class Frame:
    index: int
    timepoint: int
    filename: str
    path: Path
    sha256: str
    width: int
    height: int
    dtype: str
    bits_per_sample: tuple[int, ...]
    channels: int
    pages: int
    photometric: str
    planar_configuration: str
    orientation: int


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except OSError as exc:
        raise AnalysisError(f"Could not write {path}: {exc}") from exc


def write_json(path: Path, value: object) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def threshold_value(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("threshold must be an integer") from exc
    if not 0 <= value <= 255:
        raise argparse.ArgumentTypeError("threshold must be between 0 and 255")
    return value


def positive_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a number") from exc
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return value


def discover_frames(
    input_directory: Path,
    frame_pattern: str = FRAME_PATTERN,
    allow_time_gaps: bool = False,
) -> list[tuple[int, Path]]:
    directory = input_directory.expanduser().resolve()
    if not directory.is_dir():
        raise AnalysisError(f"Input directory does not exist: {directory}")
    try:
        pattern = re.compile(frame_pattern, re.IGNORECASE)
    except re.error as exc:
        raise AnalysisError(f"Invalid frame pattern: {exc}") from exc
    if "time" not in pattern.groupindex:
        raise AnalysisError("The frame pattern must have a named 'time' group")

    found: list[tuple[int, Path]] = []
    for entry in os.scandir(directory):
        if entry.name.startswith("._") or not entry.is_file(follow_symlinks=False):
            continue
        match = pattern.fullmatch(entry.name)
        if match:
            found.append((int(match.group("time")), Path(entry.path).resolve()))

    if len(found) < 2:
        raise AnalysisError("At least two matching TIFF frames are required")
    found.sort(key=lambda item: (item[0], item[1].name.casefold()))

    timepoints = [item[0] for item in found]
    duplicates = sorted({value for value in timepoints if timepoints.count(value) > 1})
    if duplicates:
        raise AnalysisError(f"Duplicate timepoint numbers: {duplicates}")
    missing = sorted(set(range(timepoints[0], timepoints[-1] + 1)) - set(timepoints))
    if missing and not allow_time_gaps:
        raise AnalysisError(
            f"Missing timepoint numbers: {missing}. Confirm the gap and use "
            "--allow-time-gaps if it is intentional."
        )
    return found


def _bits_per_sample(value: object, channels: int) -> tuple[int, ...]:
    if isinstance(value, (tuple, list)):
        return tuple(int(item) for item in value)
    return (int(value),) * channels


def inspect_frame(index: int, timepoint: int, path: Path) -> Frame:
    try:
        before_hash = sha256_file(path)
        with tifffile.TiffFile(path) as tif:
            pages = len(tif.pages)
            if pages != 1:
                raise AnalysisError(f"{path.name}: expected one TIFF page, found {pages}")
            page = tif.pages[0]
            shape = tuple(int(value) for value in page.shape)
            dtype = np.dtype(page.dtype)
            channels = int(page.samplesperpixel or 1)
            bits = _bits_per_sample(page.bitspersample, channels)
            photometric = str(getattr(page.photometric, "name", page.photometric))
            planar = str(getattr(page.planarconfig, "name", page.planarconfig))
            orientation_tag = page.tags.get("Orientation")
            orientation = int(orientation_tag.value) if orientation_tag else 1
            if page.subifds:
                raise AnalysisError(f"{path.name}: SubIFDs are not supported")
        after_hash = sha256_file(path)
    except AnalysisError:
        raise
    except Exception as exc:
        raise AnalysisError(f"Could not read TIFF metadata from {path}: {exc}") from exc
    if before_hash != after_hash:
        raise AnalysisError(f"Input changed while metadata were being read: {path}")

    problems: list[str] = []
    if len(shape) != 3 or shape[2] != 3:
        problems.append(f"shape is {shape}, not H x W x 3")
    if dtype != np.dtype("uint8"):
        problems.append(f"dtype is {dtype}, not uint8")
    if channels != 3:
        problems.append(f"channels is {channels}, not 3")
    if bits != (8, 8, 8):
        problems.append(f"bits per sample is {bits}, not (8, 8, 8)")
    if photometric.upper() != "RGB":
        problems.append(f"photometric is {photometric}, not RGB")
    if planar.upper() != "CONTIG":
        problems.append(f"planar configuration is {planar}, not CONTIG")
    if orientation != 1:
        problems.append(f"orientation is {orientation}, not 1")
    if problems:
        raise AnalysisError(f"Unsupported input {path.name}: " + "; ".join(problems))

    return Frame(
        index=index,
        timepoint=timepoint,
        filename=path.name,
        path=path,
        sha256=after_hash,
        width=shape[1],
        height=shape[0],
        dtype=dtype.name,
        bits_per_sample=bits,
        channels=channels,
        pages=pages,
        photometric=photometric,
        planar_configuration=planar,
        orientation=orientation,
    )


def inspect_frames(discovered: list[tuple[int, Path]]) -> list[Frame]:
    frames = [
        inspect_frame(index, timepoint, path)
        for index, (timepoint, path) in enumerate(discovered, start=1)
    ]
    sizes = {(frame.width, frame.height) for frame in frames}
    if len(sizes) != 1:
        details = ", ".join(
            f"{frame.filename}={frame.width}x{frame.height}" for frame in frames
        )
        raise AnalysisError(f"Frame dimensions do not match: {details}")
    return frames


def parse_roi(text: str | None, width: int, height: int) -> tuple[int, int, int, int]:
    if text is None:
        return 0, 0, width, height
    pieces = [piece.strip() for piece in text.split(",")]
    if len(pieces) != 4:
        raise AnalysisError("ROI must contain four integers: x1,y1,x2,y2")
    try:
        x1, y1, x2, y2 = (int(piece) for piece in pieces)
    except ValueError as exc:
        raise AnalysisError("ROI must contain four integers: x1,y1,x2,y2") from exc
    if not (0 <= x1 < x2 <= width and 0 <= y1 < y2 <= height):
        raise AnalysisError(
            f"ROI {(x1, y1, x2, y2)} is outside the {width}x{height} image"
        )
    return x1, y1, x2, y2


def gaussian_blur_5x5(image: np.ndarray) -> np.ndarray:
    """Match OpenCV GaussianBlur(image, (5, 5), 0) for uint8 RGB data."""
    if image.dtype != np.uint8 or image.ndim != 3 or image.shape[2] != 3:
        raise AnalysisError(f"Expected uint8 RGB image, found {image.dtype} {image.shape}")
    if image.shape[0] < 3 or image.shape[1] < 3:
        raise AnalysisError("Images must be at least 3 by 3 pixels")

    work = image.astype(np.uint32, copy=False)
    padded_x = np.pad(work, ((0, 0), (2, 2), (0, 0)), mode="reflect")
    horizontal = sum(
        weight * padded_x[:, offset : offset + image.shape[1], :]
        for offset, weight in enumerate(GAUSSIAN_WEIGHTS)
    )
    padded_y = np.pad(horizontal, ((2, 2), (0, 0), (0, 0)), mode="reflect")
    total = sum(
        weight * padded_y[offset : offset + image.shape[0], :, :]
        for offset, weight in enumerate(GAUSSIAN_WEIGHTS)
    )
    return ((total + 128) // 256).astype(np.uint8)


def load_frame(frame: Frame) -> np.ndarray:
    try:
        if sha256_file(frame.path) != frame.sha256:
            raise AnalysisError(f"Input changed after inspection: {frame.path}")
        with tifffile.TiffFile(frame.path) as tif:
            image = tif.pages[0].asarray()
        after_hash = sha256_file(frame.path)
    except AnalysisError:
        raise
    except Exception as exc:
        raise AnalysisError(f"Could not read {frame.path}: {exc}") from exc
    if image.shape != (frame.height, frame.width, 3) or image.dtype != np.uint8:
        raise AnalysisError(f"Input changed after inspection: {frame.path}")
    if after_hash != frame.sha256:
        raise AnalysisError(f"Input changed while it was being read: {frame.path}")
    return image


def compute_unions(
    frames: list[Frame],
    roi: tuple[int, int, int, int],
    threshold: int,
    include_direction: bool,
) -> dict[str, np.ndarray]:
    x1, y1, x2, y2 = roi
    absolute = np.zeros((y2 - y1, x2 - x1), dtype=bool)
    increase = np.zeros_like(absolute)
    decrease = np.zeros_like(absolute)
    previous: np.ndarray | None = None

    for frame in frames:
        # The order below is part of the method and must not be rearranged.
        blurred_full_frame = gaussian_blur_5x5(load_frame(frame))
        red_roi = blurred_full_frame[y1:y2, x1:x2, 0]
        if previous is not None:
            difference = red_roi.astype(np.int16) - previous.astype(np.int16)
            absolute |= np.abs(difference) > threshold
            if include_direction:
                increase |= difference > threshold
                decrease |= difference < -threshold
        previous = red_roi

    result = {"absolute": absolute}
    if include_direction:
        result.update({"increase": increase, "decrease": decrease})
    return result


def prepare_output(input_directory: Path, output_directory: Path) -> Path:
    source = input_directory.expanduser().resolve()
    output = output_directory.expanduser().resolve()
    if output == source or source in output.parents:
        raise AnalysisError("Output must be separate from the input directory")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise AnalysisError(f"Output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def save_mask(path: Path, mask: np.ndarray) -> None:
    data = np.where(mask, 255, 0).astype(np.uint8)
    temporary = path.with_name(path.stem + ".tmp" + path.suffix)
    try:
        tifffile.imwrite(
            temporary,
            data,
            photometric="minisblack",
            compression=None,
            metadata=None,
            software=False,
        )
        os.replace(temporary, path)
    except OSError as exc:
        raise AnalysisError(f"Could not write {path}: {exc}") from exc
    restored = tifffile.imread(path)
    if not np.array_equal(data, restored):
        raise AnalysisError(f"Saved mask failed verification: {path}")


def save_pseudocolor(
    path: Path,
    increase: np.ndarray,
    decrease: np.ndarray,
    absolute: np.ndarray,
) -> None:
    # RGB: magenta=increase only, cyan=decrease only, white=both, black=none.
    image = np.stack((increase, decrease, absolute), axis=-1).astype(np.uint8) * 255
    temporary = path.with_name(path.stem + ".tmp" + path.suffix)
    try:
        Image.fromarray(image).save(temporary)
        os.replace(temporary, path)
    except OSError as exc:
        raise AnalysisError(f"Could not write {path}: {exc}") from exc
    with Image.open(path) as restored:
        restored_array = np.asarray(restored.convert("RGB"))
    if not np.array_equal(image, restored_array):
        raise AnalysisError(f"Saved pseudocolor image failed verification: {path}")


def percent(count: int, total: int) -> float:
    return 100.0 * float(count) / float(total)


def run_analysis(
    version: str,
    input_directory: Path,
    output_directory: Path,
    threshold: int = 0,
    roi_text: str | None = None,
    frame_pattern: str = FRAME_PATTERN,
    allow_time_gaps: bool = False,
    parameters_confirmed: bool = False,
    dataset_label: str | None = None,
    frame_interval_seconds: float | None = None,
) -> dict[str, object]:
    if version not in {"v1", "v2"}:
        raise AnalysisError("Version must be v1 or v2")
    if not 0 <= threshold <= 255:
        raise AnalysisError("Threshold must be between 0 and 255")
    if frame_interval_seconds is not None and (
        not math.isfinite(frame_interval_seconds) or frame_interval_seconds <= 0
    ):
        raise AnalysisError("Frame interval must be greater than zero")

    started = time.perf_counter()
    started_at_utc = datetime.now(timezone.utc).isoformat()
    discovered = discover_frames(input_directory, frame_pattern, allow_time_gaps)
    frames = inspect_frames(discovered)
    roi = parse_roi(roi_text, frames[0].width, frames[0].height)
    output = prepare_output(input_directory, output_directory)
    x1, y1, x2, y2 = roi
    roi_pixels = (x2 - x1) * (y2 - y1)

    processing_started = time.perf_counter()
    unions = compute_unions(frames, roi, threshold, version == "v2")
    processing_seconds = time.perf_counter() - processing_started

    absolute = unions["absolute"]
    save_mask(output / "absolute_change_union.tif", absolute)
    absolute_count = int(np.count_nonzero(absolute))
    metrics: dict[str, object] = {
        "version": version,
        "dataset": dataset_label or input_directory.name,
        "frame_count": len(frames),
        "transition_count": len(frames) - 1,
        "frame_interval_seconds": frame_interval_seconds,
        "roi_pixels": roi_pixels,
        "absolute_change_pixels": absolute_count,
        "absolute_change_percent": percent(absolute_count, roi_pixels),
        "processing_seconds": processing_seconds,
    }

    if version == "v2":
        increase = unions["increase"]
        decrease = unions["decrease"]
        overlap = increase & decrease
        direction_or = increase | decrease
        mismatch = int(np.count_nonzero(absolute ^ direction_or))
        increase_count = int(np.count_nonzero(increase))
        decrease_count = int(np.count_nonzero(decrease))
        overlap_count = int(np.count_nonzero(overlap))
        save_mask(output / "increase_union.tif", increase)
        save_mask(output / "decrease_union.tif", decrease)
        save_pseudocolor(
            output / "pseudo_colormap_diff.png", increase, decrease, absolute
        )
        metrics.update(
            {
                "increase_pixels": increase_count,
                "increase_percent": percent(increase_count, roi_pixels),
                "decrease_pixels": decrease_count,
                "decrease_percent": percent(decrease_count, roi_pixels),
                "overlap_pixels": overlap_count,
                "overlap_percent": percent(overlap_count, roi_pixels),
                "increase_only_pixels": int(np.count_nonzero(increase & ~decrease)),
                "decrease_only_pixels": int(np.count_nonzero(decrease & ~increase)),
                "absolute_equals_increase_or_decrease": mismatch == 0,
                "identity_mismatch_pixels": mismatch,
                "inclusion_exclusion_residual_pixels": absolute_count
                - (increase_count + decrease_count - overlap_count),
            }
        )
        if mismatch or metrics["inclusion_exclusion_residual_pixels"] != 0:
            raise AnalysisError("Internal identity failed: absolute != increase OR decrease")

    order_text = "index\ttimepoint\tfilename\tsha256\n" + "".join(
        f"{frame.index}\t{frame.timepoint}\t{frame.filename}\t{frame.sha256}\n"
        for frame in frames
    )
    atomic_write_text(output / "frame_order.tsv", order_text)
    signature_source = {
        "algorithm_version": ALGORITHM_VERSION,
        "frames": [(frame.timepoint, frame.filename, frame.sha256) for frame in frames],
        "threshold": threshold,
        "roi": roi,
        "frame_interval_seconds": frame_interval_seconds,
        "algorithm": ALGORITHM,
        "gaussian_weights": GAUSSIAN_WEIGHTS,
        "gaussian_border": "OpenCV BORDER_REFLECT_101 equivalent",
        "source_channel": "RGB index 0",
        "input_contract": "single-page uint8 contiguous RGB, orientation 1",
    }
    signature = hashlib.sha256(
        json.dumps(signature_source, sort_keys=True).encode("utf-8")
    ).hexdigest()
    parameters = {
        "version": version,
        "dataset": dataset_label or input_directory.name,
        "input_directory": str(input_directory.expanduser().resolve()),
        "output_directory": str(output),
        "parameter_status": (
            "confirmed" if parameters_confirmed else "runtime_test_only_unconfirmed"
        ),
        "algorithm_version": ALGORITHM_VERSION,
        "frame_count": len(frames),
        "frame_interval_seconds": frame_interval_seconds,
        "image_width": frames[0].width,
        "image_height": frames[0].height,
        "image_dtype": frames[0].dtype,
        "bits_per_sample": frames[0].bits_per_sample,
        "channels": frames[0].channels,
        "pages_per_file": frames[0].pages,
        "photometric": frames[0].photometric,
        "planar_configuration": frames[0].planar_configuration,
        "orientation": frames[0].orientation,
        "threshold": threshold,
        "threshold_rule": "strict post-blur red-channel difference > threshold",
        "roi_half_open_xyxy": roi,
        "frame_pattern": frame_pattern,
        "frame_sort": "integer value of the named time group",
        "algorithm": ALGORITHM,
        "gaussian_kernel_1d": [value / 16 for value in GAUSSIAN_WEIGHTS],
        "gaussian_border": "OpenCV BORDER_REFLECT_101 equivalent",
        "source_channel": "RGB red, array index 0; equivalent to OpenCV BGR index 2",
        "comparison_signature": signature,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "tifffile": tifffile.__version__,
            "Pillow": pillow_version,
        },
    }
    metrics["comparison_signature"] = signature
    metrics["total_seconds"] = time.perf_counter() - started
    write_json(output / "parameters.json", parameters)
    write_json(output / "metrics.json", metrics)
    atomic_write_text(
        output / "run.log",
        "\n".join(
            [
                f"started_utc={started_at_utc}",
                f"version={version}",
                f"frames={len(frames)}",
                f"frame_order={','.join(frame.filename for frame in frames)}",
                f"threshold={threshold}",
                f"roi={roi}",
                f"frame_interval_seconds={frame_interval_seconds}",
                f"absolute_change_pixels={absolute_count}",
                f"processing_seconds={processing_seconds:.6f}",
                "status=success",
            ]
        )
        + "\n",
    )
    return metrics


def build_parser(version: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=f"Trafficking analysis {version}")
    parser.add_argument(
        "-i", "--input-directory", "--input_directory", type=Path, required=True
    )
    parser.add_argument("-o", "--output-directory", type=Path, required=True)
    parser.add_argument(
        "-t",
        "--threshold",
        "--threadshold",
        type=threshold_value,
        default=0,
        help="integer threshold from 0 to 255; differences must be greater than it",
    )
    parser.add_argument(
        "-roi",
        "--roi",
        "--region-of-interest",
        "--region_of_interest",
        dest="roi_text",
        help="half-open coordinates x1,y1,x2,y2; default is the full image",
    )
    parser.add_argument("--frame-pattern", default=FRAME_PATTERN)
    parser.add_argument("--allow-time-gaps", action="store_true")
    parser.add_argument(
        "--frame-interval-seconds",
        type=positive_float,
        help="optional acquisition interval recorded as metadata; it does not change masks",
    )
    parser.add_argument("--parameters-confirmed", action="store_true")
    parser.add_argument("--dataset-label")
    return parser


def cli_main(version: str, argv: Sequence[str] | None = None) -> int:
    parser = build_parser(version)
    args = parser.parse_args(argv)
    try:
        metrics = run_analysis(
            version=version,
            input_directory=args.input_directory,
            output_directory=args.output_directory,
            threshold=args.threshold,
            roi_text=args.roi_text,
            frame_pattern=args.frame_pattern,
            allow_time_gaps=args.allow_time_gaps,
            parameters_confirmed=args.parameters_confirmed,
            dataset_label=args.dataset_label,
            frame_interval_seconds=args.frame_interval_seconds,
        )
    except AnalysisError as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit("Run trafficking_v1.py or trafficking_v2.py")
