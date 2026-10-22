# Contributing

This repository has one author and a very specific defect it is built to avoid: a
number stated in prose that no longer matches the run that produced it. Every
contribution is judged against that first.

## The invariant to preserve

`results/ledger.jsonl` is append-only. A row is never edited or removed; a re-run
appends. `python -m mspine verify` re-derives every row's hash and its parent's,
and compares the head with `results/anchor.json`. If you find yourself wanting to
delete a row because it is wrong, the correct action is to append the corrected
row and let the history show both.

## If you change a number

Any change that moves a result has to move three things together, or CI fails:

    the code            mspine/ -- and the benchmark that produced the row
    the chain           results/ledger.jsonl, appended, never rewritten
    the anchor          python -m mspine anchor, committed with the chain
    the prose           README.md / docs/paper.md, checked by tests/test_docs.py

`tests/test_docs.py` parses the tables in the README and compares them with the
committed chain, so editing a number in the prose to look better than the run
fails in under a second. `tests/test_published.py` re-derives the numbers by
re-running the campaigns, and is marked `slow` and run by CI.

## What gets merged

    fixes to a measurement      always -- especially when the old number flattered us
    a new acquisition arm       yes, if it runs inside the same budget and control
    a new table                 yes, with its upstream licence, its SHA-256 pinned
                                in mspine/tables/PROVENANCE.md, and a campaign on it
    a new claim                 only with the control that could kill it

The last line is the whole review. "It works" is not a result; "it works and here
is the run in which it should have failed and did not" is.

## Running it

    python -m pip install numpy scipy pytest
    python -m pytest tests -q          # ~20 s, the invariants
    python -m pytest tests -q -m slow  # ~30 s, re-derives the published numbers
    python -m mspine bench --tables concrete --seeds 1 --budget 32
    python -m mspine verify

If the fast suite fails on a machine with more than one core, set
`OPENBLAS_NUM_THREADS=1`. This host caps a single allocation near 10 MB and an
unrestricted BLAS thread pool dies on a 300x300 solve that has no business
failing; `conftest.py` pins it for the suite, but a script you run by hand does
not get that.

## Licence

MIT. Contributions are accepted under it.
