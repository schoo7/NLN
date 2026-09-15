#!/usr/bin/env python3
"""Safe revision of trafficking_v2.py."""

import sys

if __package__:
    from .trafficking_core import cli_main
else:
    from trafficking_core import cli_main


if __name__ == "__main__":
    raise SystemExit(cli_main("v2", sys.argv[1:]))
