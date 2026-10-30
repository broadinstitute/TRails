"""TRails: one command from your tables to a running results server.

    python3 trails.py \\
        --repeat-copy-numbers-tsv  your_cohort.tsv.gz \\
        --sample-metadata-tsv      your_samples.tsv.gz \\
        [--phenotypes-table        your_phenotypes.tsv.gz]

This file only puts the trails/ code directory on the import path and hands over to
trails/trails_launcher.py, which does the work. Run with --help for every option.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "trails"))

import trails_launcher  # noqa: E402 - needs the path set up above

if __name__ == "__main__":
    trails_launcher.main()
