from __future__ import annotations

import argparse
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from run_demo import validate_output_root
from src.trafficking_core import (
    AnalysisError,
    cli_main,
    compute_unions,
    discover_frames,
    gaussian_blur_5x5,
    inspect_frames,
    parse_roi,
    positive_float,
    prepare_output,
    run_analysis,
    sha256_file,
    threshold_value,
)


def write_rgb(path: Path, red: int, green: int = 0, blue: int = 0) -> None:
    image = np.zeros((7, 9, 3), dtype=np.uint8)
    image[:, :, 0] = red
    image[:, :, 1] = green
    image[:, :, 2] = blue
    tifffile.imwrite(path, image, photometric="rgb", metadata=None, software=False)


def write_array(path: Path, image: np.ndarray) -> None:
    tifffile.imwrite(path, image, photometric="rgb", metadata=None, software=False)


class TraffickingTests(unittest.TestCase):
    def test_committed_demo_frames_are_metadata_stripped(self) -> None:
        root = Path(__file__).resolve().parents[1] / "demo_data"
        expected = {
            Path("control/control_T01_RGB_TRITC.tif"): (
                "8663589d1efd55418fbcde11dde080f006c77fc1c56eed6529b0e303a053ae5d"
            ),
            Path("control/control_T02_RGB_TRITC.tif"): (
                "b5e80c80746a16ef93ccaa85eeb5417b39bbc1cb6a85c5c45a0b3a05174bf9ad"
            ),
            Path("experiment/experiment_T01_RGB_TRITC.tif"): (
                "86ab7815a82bb62353e4b20ae21e01e868853389cde9ea93e26e5f6958e89b1e"
            ),
            Path("experiment/experiment_T02_RGB_TRITC.tif"): (
                "6bd4c0a3f9cbf6e93f75946b23f036b3c91e058451cbea5e362582158ba7a55b"
            ),
        }
        allowed_tags = {
            256,
            257,
            258,
            259,
            262,
            273,
            274,
            277,
            278,
            279,
            282,
            283,
            284,
            296,
        }
        for relative_path, expected_hash in expected.items():
            path = root / relative_path
            with self.subTest(frame=str(relative_path)):
                self.assertEqual(sha256_file(path), expected_hash)
                with tifffile.TiffFile(path) as tif:
                    self.assertEqual(len(tif.pages), 1)
                    page = tif.pages[0]
                    self.assertEqual(page.shape, (512, 512, 3))
                    self.assertEqual(page.dtype, np.dtype("uint8"))
                    self.assertEqual(page.photometric.name, "RGB")
                    self.assertEqual(page.planarconfig.name, "CONTIG")
                    self.assertEqual(
                        {tag.code for tag in page.tags.values()}, allowed_tags
                    )
                    self.assertEqual(page.tags[274].value, 1)

        expected_counts = {
            "control": (123344, 64493, 58851, 0),
            "experiment": (185812, 94471, 91341, 0),
        }
        for series, counts in expected_counts.items():
            frames = inspect_frames(discover_frames(root / series))
            roi = parse_roi(None, width=512, height=512)
            v1 = compute_unions(frames, roi, threshold=0, include_direction=False)
            v2 = compute_unions(frames, roi, threshold=0, include_direction=True)
            actual = (
                int(np.count_nonzero(v1["absolute"])),
                int(np.count_nonzero(v2["increase"])),
                int(np.count_nonzero(v2["decrease"])),
                int(np.count_nonzero(v2["increase"] & v2["decrease"])),
            )
            with self.subTest(series=series):
                self.assertEqual(actual, counts)
                self.assertTrue(np.array_equal(v1["absolute"], v2["absolute"]))
                self.assertTrue(
                    np.array_equal(
                        v2["absolute"], v2["increase"] | v2["decrease"]
                    )
                )

    def test_numeric_time_sort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "sample_T10_RGB_TRITC.tif", 10)
            write_rgb(directory / "sample_T1_RGB_TRITC.tif", 1)
            write_rgb(directory / "sample_T2_RGB_TRITC.tif", 2)
            found = discover_frames(directory, allow_time_gaps=True)
            self.assertEqual([item[0] for item in found], [1, 2, 10])

    def test_v1_v2_and_direction_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "sample_T01_RGB_TRITC.tif", 0)
            write_rgb(directory / "sample_T02_RGB_TRITC.tif", 100)
            write_rgb(directory / "sample_T03_RGB_TRITC.tif", 0)
            frames = inspect_frames(discover_frames(directory))
            roi = parse_roi(None, 9, 7)
            v1 = compute_unions(frames, roi, threshold=0, include_direction=False)
            v2 = compute_unions(frames, roi, threshold=0, include_direction=True)
            self.assertTrue(np.array_equal(v1["absolute"], v2["absolute"]))
            self.assertTrue(np.all(v2["increase"]))
            self.assertTrue(np.all(v2["decrease"]))
            self.assertTrue(
                np.array_equal(v2["absolute"], v2["increase"] | v2["decrease"])
            )

    def test_only_source_red_is_analyzed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "sample_T01_RGB_TRITC.tif", 20, 0, 0)
            write_rgb(directory / "sample_T02_RGB_TRITC.tif", 20, 200, 200)
            frames = inspect_frames(discover_frames(directory))
            result = compute_unions(
                frames, parse_roi(None, 9, 7), threshold=0, include_direction=False
            )
            self.assertFalse(np.any(result["absolute"]))

    def test_invalid_roi_is_rejected(self) -> None:
        with self.assertRaises(AnalysisError):
            parse_roi("-1,0,5,5", width=9, height=7)
        with self.assertRaises(AnalysisError):
            parse_roi("4,4,4,5", width=9, height=7)

    def test_uint16_is_rejected_without_conversion(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            image = np.zeros((7, 9, 3), dtype=np.uint16)
            for timepoint in (1, 2):
                tifffile.imwrite(
                    directory / f"sample_T{timepoint:02d}_RGB_TRITC.tif",
                    image,
                    photometric="rgb",
                    metadata=None,
                    software=False,
                )
            with self.assertRaisesRegex(AnalysisError, "uint8"):
                inspect_frames(discover_frames(directory))

    def test_integration_outputs_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_directory = root / "input"
            input_directory.mkdir()
            write_rgb(input_directory / "sample_T01_RGB_TRITC.tif", 0)
            write_rgb(input_directory / "sample_T02_RGB_TRITC.tif", 10)
            metrics = run_analysis(
                "v2", input_directory, root / "output", threshold=0
            )
            self.assertTrue(metrics["absolute_equals_increase_or_decrease"])
            self.assertEqual(metrics["absolute_change_pixels"], 63)
            self.assertTrue((root / "output" / "frame_order.tsv").is_file())
            self.assertFalse((input_directory / "absolute_change_union.tif").exists())

    def test_demo_output_cannot_enter_either_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            control = root / "control"
            experiment = root / "experiment"
            control.mkdir()
            experiment.mkdir()
            with self.assertRaises(AnalysisError):
                validate_output_root(experiment / "new_output", [control, experiment])

    def test_duplicate_and_missing_timepoints_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "a_T1_RGB_TRITC.tif", 0)
            write_rgb(directory / "b_T1_RGB_TRITC.tif", 0)
            with self.assertRaisesRegex(AnalysisError, "Duplicate"):
                discover_frames(directory)

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "sample_T1_RGB_TRITC.tif", 0)
            write_rgb(directory / "sample_T3_RGB_TRITC.tif", 0)
            with self.assertRaisesRegex(AnalysisError, "Missing"):
                discover_frames(directory)

    def test_threshold_rule_is_strict(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_rgb(directory / "sample_T1_RGB_TRITC.tif", 0)
            write_rgb(directory / "sample_T2_RGB_TRITC.tif", 10)
            frames = inspect_frames(discover_frames(directory))
            equal = compute_unions(
                frames, parse_roi(None, 9, 7), threshold=10, include_direction=False
            )
            above = compute_unions(
                frames, parse_roi(None, 9, 7), threshold=9, include_direction=False
            )
            self.assertFalse(np.any(equal["absolute"]))
            self.assertTrue(np.all(above["absolute"]))

    def test_blur_occurs_before_roi_crop(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            frame_1 = np.zeros((7, 7, 3), dtype=np.uint8)
            frame_2 = frame_1.copy()
            frame_2[3, 2, 0] = 255
            write_array(directory / "sample_T1_RGB_TRITC.tif", frame_1)
            write_array(directory / "sample_T2_RGB_TRITC.tif", frame_2)
            result = compute_unions(
                inspect_frames(discover_frames(directory)),
                (3, 3, 4, 4),
                threshold=0,
                include_direction=False,
            )
            self.assertEqual(result["absolute"].shape, (1, 1))
            self.assertTrue(bool(result["absolute"][0, 0]))

    def test_gaussian_blur_preserves_constants_and_is_deterministic(self) -> None:
        image = np.full((9, 11, 3), 137, dtype=np.uint8)
        first = gaussian_blur_5x5(image)
        second = gaussian_blur_5x5(image)
        self.assertTrue(np.array_equal(first, image))
        self.assertTrue(np.array_equal(first, second))

    def test_inconsistent_dimensions_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_array(
                directory / "sample_T1_RGB_TRITC.tif",
                np.zeros((8, 8, 3), dtype=np.uint8),
            )
            write_array(
                directory / "sample_T2_RGB_TRITC.tif",
                np.zeros((9, 8, 3), dtype=np.uint8),
            )
            with self.assertRaisesRegex(AnalysisError, "dimensions do not match"):
                inspect_frames(discover_frames(directory))

    def test_grayscale_and_multipage_inputs_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            grayscale = root / "grayscale"
            grayscale.mkdir()
            for timepoint in (1, 2):
                tifffile.imwrite(
                    grayscale / f"sample_T{timepoint}_RGB_TRITC.tif",
                    np.zeros((8, 8), dtype=np.uint8),
                    metadata=None,
                    software=False,
                )
            with self.assertRaisesRegex(AnalysisError, "not H x W x 3"):
                inspect_frames(discover_frames(grayscale))

            multipage = root / "multipage"
            multipage.mkdir()
            for timepoint in (1, 2):
                path = multipage / f"sample_T{timepoint}_RGB_TRITC.tif"
                with tifffile.TiffWriter(path) as writer:
                    writer.write(np.zeros((8, 8, 3), dtype=np.uint8), photometric="rgb")
                    writer.write(np.zeros((8, 8, 3), dtype=np.uint8), photometric="rgb")
            with self.assertRaisesRegex(AnalysisError, "expected one TIFF page"):
                inspect_frames(discover_frames(multipage))

    def test_corrupt_matching_tiff_reports_read_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "sample_T1_RGB_TRITC.tif").write_bytes(b"not a TIFF")
            (directory / "sample_T2_RGB_TRITC.tif").write_bytes(b"not a TIFF")
            with self.assertRaisesRegex(AnalysisError, "Could not read TIFF metadata"):
                inspect_frames(discover_frames(directory))

    def test_frozen_input_hash_is_rechecked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            first = directory / "sample_T1_RGB_TRITC.tif"
            second = directory / "sample_T2_RGB_TRITC.tif"
            write_rgb(first, 0)
            write_rgb(second, 10)
            frames = inspect_frames(discover_frames(directory))
            write_rgb(second, 20)
            with self.assertRaisesRegex(AnalysisError, "changed after inspection"):
                compute_unions(
                    frames,
                    parse_roi(None, 9, 7),
                    threshold=0,
                    include_direction=False,
                )

    def test_output_directory_safety(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "input"
            source.mkdir()
            with self.assertRaises(AnalysisError):
                prepare_output(source, source)
            with self.assertRaises(AnalysisError):
                prepare_output(source, source / "inside")
            occupied = Path(temporary) / "occupied"
            occupied.mkdir()
            (occupied / "existing.txt").write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(AnalysisError, "not empty"):
                prepare_output(source, occupied)

    def test_cli_value_parsers(self) -> None:
        self.assertEqual(threshold_value("0"), 0)
        self.assertEqual(threshold_value("255"), 255)
        self.assertEqual(positive_float("30"), 30.0)
        for value in ("-1", "256", "abc", "1.5"):
            with self.subTest(threshold=value), self.assertRaises(
                argparse.ArgumentTypeError
            ):
                threshold_value(value)
        for value in ("0", "-1", "nan", "inf", "abc"):
            with self.subTest(interval=value), self.assertRaises(
                argparse.ArgumentTypeError
            ):
                positive_float(value)

    def test_directional_outputs_and_pseudocolor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input"
            source.mkdir()
            frame_1 = np.zeros((32, 32, 3), dtype=np.uint8)
            frame_2 = frame_1.copy()
            frame_3 = frame_1.copy()
            frame_1[4:12, 4:12, 0] = 100
            frame_2[14:22, 4:12, 0] = 100
            frame_3[14:22, 4:12, 0] = 100
            frame_2[14:22, 20:28, 0] = 100
            for index, frame in enumerate((frame_1, frame_2, frame_3), start=1):
                write_array(source / f"sample_T{index}_RGB_TRITC.tif", frame)

            v1_output = root / "v1"
            v2_output = root / "v2"
            v1_metrics = run_analysis("v1", source, v1_output)
            v2_metrics = run_analysis("v2", source, v2_output)
            v1_absolute = tifffile.imread(v1_output / "absolute_change_union.tif")
            v2_absolute = tifffile.imread(v2_output / "absolute_change_union.tif")
            increase = tifffile.imread(v2_output / "increase_union.tif") == 255
            decrease = tifffile.imread(v2_output / "decrease_union.tif") == 255
            self.assertTrue(np.array_equal(v1_absolute, v2_absolute))
            self.assertTrue(np.array_equal(v2_absolute == 255, increase | decrease))
            self.assertEqual(
                v1_metrics["absolute_change_pixels"],
                v2_metrics["absolute_change_pixels"],
            )
            self.assertGreater(v2_metrics["overlap_pixels"], 0)
            pseudocolor = np.asarray(
                Image.open(v2_output / "pseudo_colormap_diff.png").convert("RGB")
            )
            self.assertEqual(tuple(pseudocolor[8, 8]), (0, 255, 255))
            self.assertEqual(tuple(pseudocolor[18, 8]), (255, 0, 255))
            self.assertEqual(tuple(pseudocolor[18, 24]), (255, 255, 255))
            self.assertEqual(tuple(pseudocolor[28, 16]), (0, 0, 0))

    def test_cli_uses_supplied_arguments_and_records_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input"
            source.mkdir()
            write_rgb(source / "sample_T1_RGB_TRITC.tif", 0)
            write_rgb(source / "sample_T2_RGB_TRITC.tif", 10)
            output = root / "output"
            result = cli_main(
                "v1",
                [
                    "-i",
                    str(source),
                    "-o",
                    str(output),
                    "--threshold",
                    "0",
                    "--frame-interval-seconds",
                    "30",
                ],
            )
            self.assertEqual(result, 0)
            parameters = json.loads(
                (output / "parameters.json").read_text(encoding="utf-8")
            )
            self.assertEqual(parameters["input_directory"], str(source.resolve()))
            self.assertEqual(parameters["threshold"], 0)
            self.assertEqual(parameters["frame_interval_seconds"], 30.0)


if __name__ == "__main__":
    unittest.main()
