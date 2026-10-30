"""TRails tests. Run from the TRails directory with: .venv/bin/python -m unittest discover -p '*tests.py'

The code lives in ../trails/ and its modules import each other by bare name, so it goes on the
import path here, before any test module is loaded.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "trails"))
