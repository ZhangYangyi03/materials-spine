"""Pin BLAS thread pools before anything imports numpy.

This host caps a single allocation near 10 MB and OpenBLAS spawns a thread pool
per call, so an unrestricted build dies with "Memory allocation still failed
after 10 retries" on a 300x300 solve that has no business failing. mspine/__init__
does the same pinning, but a test module that does `import numpy` on its first
line gets there first -- which is exactly how the whole suite aborted with no
output at all. A conftest is imported by pytest before any test module, so this
is the only place the pin is early enough. Set MSPINE_ALLOW_MT=1 to opt out.
"""
import os

if os.environ.get("MSPINE_ALLOW_MT") != "1":
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[_v] = "1"
