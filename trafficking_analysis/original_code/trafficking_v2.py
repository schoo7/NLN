import cv2
import argparse
import sys
import os
import numpy as np


def chunk_image(image, p1, p2):
    x1, y1 = p1
    x2, y2 = p2

    return image[y1:y2, x1:x2]


def diff_red_channel_all(
    image1,
    image2,
    threshold
):

    red_c1 = image1[:, :, 2]
    red_c2 = image2[:, :, 2]

    diff = (
        red_c2.astype(np.int16)
        -
        red_c1.astype(np.int16)
    )

    increased = np.where(
        diff > int(threshold),
        255,
        0
    ).astype(np.uint8)

    decreased = np.where(
        diff < -int(threshold),
        255,
        0
    ).astype(np.uint8)

    absdiff = cv2.absdiff(
        red_c1,
        red_c2
    )

    _, absdiff_binary = cv2.threshold(
        absdiff,
        int(threshold),
        255,
        cv2.THRESH_BINARY
    )

    return (
        increased,
        decreased,
        absdiff_binary
    )


def compute_diff_percent(image):

    white_pixel_count = np.sum(
        image == 255
    )

    total_pixel_count = image.size

    return (
        white_pixel_count,
        total_pixel_count
    )


def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input_directory",
        required=True,
        help="directory containing tif images"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        required=False,
        default=0,
        help="threshold for ignoring small changes"
    )

    parser.add_argument(
        "-roi",
        "--region_of_interest",
        required=False,
        help="ROI: x1,y1,x2,y2"
    )

    args = parser.parse_args()

    input_directory = (
        args.input_directory
    )

    threshold = int(
        args.threshold
    )

    if not os.path.exists(
        input_directory
    ):

        print("input directory does not exist")

        exit(1)

    image_files = [
        f
        for f in os.listdir(
            input_directory
        )
        if f.endswith("tif")
    ]

    image_files.sort()

    if len(image_files) < 2:

        print(
            "at least two tif images are required"
        )

        exit(1)

    print(f"detected {len(image_files)} tif files")

    images = []

    for f in image_files:

        path = os.path.join(
            input_directory,
            f
        )

        image = cv2.imread(
            path
        )

        image_blur = (
            cv2.GaussianBlur(
                image,
                (5, 5),
                0
            )
        )

        images.append(
            image_blur
        )

    if args.region_of_interest:

        coords = list(
            map(
                int,
                args.region_of_interest.split(",")
            )
        )

        assert len(coords) == 4

        p1 = (
            coords[0],
            coords[1]
        )

        p2 = (
            coords[2],
            coords[3]
        )

        images = [
            chunk_image(
                img,
                p1,
                p2
            )
            for img in images
        ]

    union_inc = None
    union_dec = None
    union_abs = None

    for i in range(
        len(images) - 1
    ):

        increased, decreased, absdiff = (
            diff_red_channel_all(
                images[i],
                images[i + 1],
                threshold
            )
        )

        if union_abs is None:

            union_inc = (
                increased.copy()
            )

            union_dec = (
                decreased.copy()
            )

            union_abs = (
                absdiff.copy()
            )

        else:

            union_inc = (
                cv2.bitwise_or(
                    union_inc,
                    increased
                )
            )

            union_dec = (
                cv2.bitwise_or(
                    union_dec,
                    decreased
                )
            )

            union_abs = (
                cv2.bitwise_or(
                    union_abs,
                    absdiff
                )
            )

    white_count, total_count = (
        compute_diff_percent(
            union_abs
        )
    )

    white_ratio = (
        white_count
        /
        total_count
        *
        100
    )

    print(
        f"\033[32m"
        f"changed pixels: "
        f"{white_count}, "
        f"percentage: "
        f"{white_ratio:.2f}%"
        f"\033[0m"
    )

    cv2.imwrite(
        os.path.join(
            input_directory,
            "only_increased.tif"
        ),
        union_inc
    )

    cv2.imwrite(
        os.path.join(
            input_directory,
            "only_decreased.tif"
        ),
        union_dec
    )

    cv2.imwrite(
        os.path.join(
            input_directory,
            "total_absdiff.tif"
        ),
        union_abs
    )

    pseudo_color = cv2.merge(
        [
            union_abs,
            union_dec,
            union_inc
        ]
    )

    cv2.imwrite(
        os.path.join(
            input_directory,
            "pseudo_colormap_diff.png"
        ),
        pseudo_color
    )


if __name__ == "__main__":
    main(
        sys.argv[1:]
    )
