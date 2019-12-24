# -*- coding: utf-8 -*-
"""Miscellaneous utilities."""
from difflib import SequenceMatcher
import math
from pathlib import Path
import subprocess


def common_wildcard(s1, s2):
    """Find a common substring as a wild card."""
    match = SequenceMatcher(lambda x: x == "%", s1, s2).find_longest_match(
        0, len(s1), 0, len(s2)
    )
    if match.a > 0:
        wcard = "*"
    else:
        wcard = ""
    wcard += s1[match.a:match.size]
    if match.size == len(s1):
        wcard += ""
    else:
        wcard += "*"
    return wcard


def unit_format(value, unit="1"):
    """Pretty-print unit string."""
    if value == 1:
        if unit == "1":
            string = ""
        else:
            string = f'${str(unit).replace(" ", "$ $")}$'
            if "%" in string:
                string = string.replace("%", r"\%")
    else:
        exp = math.floor(math.log10(value))
        base = value / 10 ** exp
        if exp == 0 or exp == 1:
            string = r"${0}$".format(value)
        elif exp == -1:
            string = r"${0:0.1f}$".format(value)
        else:
            if int(base) == 1:
                string = r"$10^{{{0:d}}}$".format(int(exp))
            else:
                string = r"${0:d}\times10^{{{1:d}}}$".format(int(base), int(exp))
        if not unit == "1":
            string += r" ${}$".format(unit)
    return string


def merge_to_gif(input_images, gifname, delay=20, resize=100, output_dir=None):
    """Use list of image files to create GIF using imagemagick convert."""
    if isinstance(input_images, str):
        input_path = Path(input_images)
    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir)
    subp_args = [
        "convert",
        "-resize",
        str(resize) + "%",
        "-delay",
        str(delay),
        "-dispose",
        "Background",
        "+page",
        input_images,
        "-loop",
        "0",
        output_dir / gifname,
    ]
    print("Generating .gif: {0}".format(" ".join([str(i) for i in subp_args])))
    subprocess.check_call(subp_args)
