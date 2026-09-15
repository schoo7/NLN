import cv2
import argparse
import sys
import os
import numpy as np


def chunck_image(image, p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return image[y1:y2, x1:x2]


def dif_red_channle(image1, image2, threadshold):
    assert image1.shape == image2.shape

    red_c1 = image1[:, :, 2]
    red_c2 = image2[:, :, 2]

    red_dif = cv2.absdiff(red_c1, red_c2)

    _, red_dif_t = cv2.threshold(
        red_dif,
        int(threadshold),
        255,
        cv2.THRESH_BINARY
    )

    return red_dif_t


def comput_diff_percent(image):
    white_pixel_count = np.sum(image == 255)
    total_pixel_count = image.size

    return white_pixel_count, total_pixel_count


def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input_directory",
        required=True,
        help=(
            "the input directory which contains your input images. "
            "All images should have the same shape"
        )
    )

    parser.add_argument(
        "-t",
        "--threadshold",
        required=False,
        help=(
            "the threadshold that you want to specify "
            "to ignore the difference"
        )
    )

    parser.add_argument(
        "-roi",
        "--region_of_interest",
        required=False,
        help=(
            "specify the region of interests in the picture "
            "instead of the entire picture"
        )
    )

    args = parser.parse_args()

    input_directory = args.input_directory

    if not os.path.exists(input_directory):
        print("input directory does not exist")
        exit(1)

    image_files = [
        f for f in os.listdir(input_directory)
        if f.endswith(("tif"))
    ]

    if len(image_files) < 2:
        print("require at least 2 pictures")
        exit(1)

    print(f"detected {len(image_files)} tif files")

    images = list()

    for f in image_files:

        f = os.path.join(input_directory, f)

        image = cv2.imread(f)

        image_blur = cv2.GaussianBlur(
            image,
            (5, 5),
            0
        )

        images.append(image_blur)

    assert len(images) == len(image_files)

    threadshold = 0

    if args.threadshold:
        threadshold = args.threadshold

    if args.region_of_interest:

        region_of_interest = (
            args.region_of_interest.split(",")
        )

        assert len(region_of_interest) == 4

        rio_p = [
            int(p)
            for p in region_of_interest
        ]

        p1 = (
            rio_p[0],
            rio_p[1]
        )

        p2 = (
            rio_p[2],
            rio_p[3]
        )

        for i in range(len(images)):
            images[i] = chunck_image(
                images[i],
                p1,
                p2
            )

    union_image = None

    for i in range(len(images) - 1):

        red_dif = dif_red_channle(
            images[i],
            images[i + 1],
            threadshold
        )

        output_file = os.path.join(
            input_directory,
            f"zirui{i}.tif"
        )

        if union_image is None:
            union_image = red_dif.copy()

        else:

            assert union_image.shape == red_dif.shape

            union_image = cv2.bitwise_or(
                union_image,
                red_dif
            )

    white_pixel_count, total_pixel_count = (
        comput_diff_percent(union_image)
    )

    white_pixel_ratio = (
        float(white_pixel_count)
        / total_pixel_count
        * 100
    )

    print(
        f"\033[32m "
        f"the number of differeent pixel: "
        f"{white_pixel_count}, "
        f"different ratio: "
        f"{white_pixel_ratio} %"
        f"\033[0m"
    )

    output_file = os.path.join(
        input_directory,
        "zirui.png"
    )

    cv2.imwrite(
        output_file,
        union_image
    )


if __name__ == "__main__":
    main(sys.argv[1:])
