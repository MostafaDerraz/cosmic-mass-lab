"""`python -m cosmicmass ...`: the command line of the engine (see `runner.main`)."""

import sys

from .runner import main

sys.exit(main())
