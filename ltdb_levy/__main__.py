"""Allow ``python -m ltdb_levy`` to run the command-line interface."""

from __future__ import annotations

import sys

from .cli import main

sys.exit(main())

