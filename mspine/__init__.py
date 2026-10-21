"""mspine -- measured-table experimental design, with an audit trail.

The unit of cost in every campaign here is a *real experiment*: a row somebody
already ran and published. No surrogate in this package ever sees a simulated
objective, so a claimed gain is a statement about data the method did not get to
invent.

Four layers, each usable alone:

    mspine.data      pinned measured tables (SHA-256 verified on load)
    mspine.campaign  budgeted screening under a metered API + held-out check
    mspine.laws      engineering laws tested against permuted controls
    mspine.mechanism reaction-network reduction, with a gate that can only reject
    mspine.ledger    hash-chained record of every number this package emits

Every result is a dict with the same three things: the number, the control that
could have destroyed it, and the provenance of the bytes it came from.
"""
import os

__version__ = "1.0.0"

# This host caps one allocation near 10 MB and BLAS grabs threads per call, so an
# unrestricted OpenBLAS dies with "Memory allocation still failed after 10
# retries" on a 300x300 solve that has no business failing. Pinning the thread
# pools before numpy loads is the difference between a benchmark that runs and a
# benchmark that cannot be reproduced on the machine it was written on. Set
# MSPINE_ALLOW_MT=1 to opt out on a bigger box.
if os.environ.get("MSPINE_ALLOW_MT") != "1":
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(_v, "1")

from .space import DesignSpace, Factor
from .gp import MixedGP

__all__ = ["DesignSpace", "Factor", "MixedGP", "__version__"]
