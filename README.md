# Prostate Cancer: WNT- Mediated communication, Population Adaptation and Coordinated Therapy Resistance

For the NLN project.

This private repository contains two separate collections of research code: the existing NLN single-cell R analysis and a Python package for AI assisted spatial temporal live cell imaging analysis.

## Repository layout

- `R/` contains the original NLN single-cell analysis scripts.
- `trafficking_analysis/` contains the maintained Python image-analysis code, sanitized legacy scripts, tests, and four author-selected demo images.

## Existing R analysis

The `R/` directory is preserved unchanged from repository commit `05456f5`. It contains 34 scripts in four groups:

- `R/LNCaP_AR/` covers count import, Seurat-object creation, quality control, doublet filtering, object reduction and merging, Harmony integration, clustering, annotation, enrichment, trajectory analysis, reference mapping, and plotting.
- `R/patient/` covers gene-signature preparation, tumor-epithelial reclustering and annotation, GSEA, and pathway plots.
- `R/pseudotime/` contains Monocle, Slingshot, and experimental Palantir trajectory analyses and plots.
- `R/utilities/` contains helpers for gene signatures, GMT files, and GSEA plots.

All 34 files pass syntax parsing with R 4.4.2. They are historical interactive analysis scripts rather than a validated end-to-end pipeline. They depend on external study data and intermediate objects that are not included, use study-specific HPC paths, and do not record a complete set of package versions. The numeric filename prefixes reflect historical working order and alternative branches; they do not define one portable run sequence. No claim is made that the R analysis can be reproduced from this repository alone.

## Python trafficking analysis

Version 1 measures the temporal union of absolute red-channel changes. Version 2 also separates increases and decreases and creates a pseudocolor image.

The calculation order is:

1. Apply a 5 by 5 Gaussian blur to each complete RGB frame.
2. Crop the region of interest, if supplied.
3. Select the red channel.
4. Compare adjacent timepoints.
5. Keep differences strictly greater than the threshold.
6. Combine binary masks across time with a logical OR.

For Version 2, `absolute = increase OR decrease`. A pixel may occur in both directional temporal unions if it changes in opposite directions during different transitions. Therefore, `absolute count = increase count + decrease count - overlap count`.

### System requirements

- Python 3.12 or newer.
- NumPy 2.3.3.
- tifffile 2026.5.15.
- Pillow 11.3.0.

The package was tested with Python 3.12.11 on macOS 15.6.1 on Apple silicon. Windows and Linux have not been tested. No GPU or non-standard hardware is required. A normal desktop with at least 2 GB of free memory is sufficient for the included demo.

### Installation and tests

```bash
cd trafficking_analysis
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

On Windows, use `.venv\Scripts\activate`. Installation normally takes less than five minutes when Python is already installed. The expected test result is 20 passing tests, normally in less than one second on the tested computer.

### Demo

The demo contains Control and Experiment series with two 512 by 512 pixel, 8-bit RGB TIFF frames each. The images are pixel-identical, metadata-stripped copies of four experimental frames selected by the author. The original files were not modified. The committed copies do not contain OME, vendor, user, source-path, acquisition-time, instrument, or position metadata.

Run all four demo analyses:

```bash
cd trafficking_analysis
python run_demo.py
```

The command creates `demo_output/` and should finish in less than one minute on a normal desktop. With threshold 0 and the full image, the expected summary is:

| Series | Frames | ROI pixels | Absolute pixels | Absolute % | Increase | Increase % | Decrease | Decrease % | Overlap | Overlap % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Control | 2 | 262,144 | 123,344 | 47.0520% | 64,493 | 24.6021% | 58,851 | 22.4499% | 0 | 0.0000% |
| Experiment | 2 | 262,144 | 185,812 | 70.8817% | 94,471 | 36.0378% | 91,341 | 34.8438% | 0 | 0.0000% |

For both series, the Version 1 and Version 2 absolute masks must match pixel by pixel, and the Version 2 absolute mask must equal `increase OR decrease`. The expected mismatch count is zero. These demo parameters are for software verification and are not confirmed manuscript parameters or a reproduction of a paper result.

### Use with other data

Input frames must be single-page, contiguous, three-channel, 8-bit RGB TIFF files of equal dimensions. Filenames must contain numeric timepoints such as `sample_T1_RGB_TRITC.tif`, `sample_T2_RGB_TRITC.tif`, and `sample_T10_RGB_TRITC.tif`. Timepoint numbers are sorted numerically. Grayscale, 16-bit, other multichannel, multipage, duplicate-timepoint, and inconsistent-size inputs are rejected rather than silently converted.

Run Version 1:

```bash
python trafficking_analysis/src/trafficking_v1.py \
  -i /path/to/input_frames \
  -o /path/to/new_v1_output \
  -t 0
```

Run Version 2 with the same frames and parameters:

```bash
python trafficking_analysis/src/trafficking_v2.py \
  -i /path/to/input_frames \
  -o /path/to/new_v2_output \
  -t 0
```

The threshold is an integer from 0 through 255. The optional ROI is a half-open rectangle supplied as `-roi x1,y1,x2,y2`; omitting it uses the full image. The output directory must be new or empty and must not be inside the input directory.

Version 1 writes `absolute_change_union.tif`, `frame_order.tsv`, `parameters.json`, `metrics.json`, and `run.log`. Version 2 also writes `increase_union.tif`, `decrease_union.tif`, and `pseudo_colormap_diff.png`.

Real-data outputs can expose paths, filenames, hashes, and study measurements. Keep all such outputs confidential and outside this repository unless they receive a separate review and authorization.

The manuscript threshold, ROI, frame interval, frame order, and final figure-generating run still require confirmation. Results from the two demo folders alone must not be presented as a statistical significance test.

## Licensing

- The source code under `R/` is licensed under the MIT License in `LICENSE-R`. No R script was changed when this license file was added.
- Everything under `trafficking_analysis/`, including the software, documentation, tests, legacy scripts, and demo images, is covered by the Restricted Research Review License in `trafficking_analysis/LICENSE`. Redistribution is prohibited without prior written permission from Zirui Fu.

The trafficking package is proprietary and is not open-source software. This repository must remain private. Editors and reviewers require explicit repository access.
