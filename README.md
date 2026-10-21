# materials-spine

A surrogate and an audit trail for materials chemistry where every evaluation is
a **real experiment somebody already ran**. No simulator sits between the
optimiser and the objective, so a claimed gain is a statement about data the
method did not get to invent.

This is the merge of three earlier repositories, and it exists because splitting
them was the mistake. They were one pipeline cut into three pieces:

    process-spine     the surrogate, the metered campaign, four measured tables
    formula-screen    the engineering laws, the carbon-fibre table, the audit layer
    kinetics-spine    the reaction-network reduction, with a gate

The same Buchwald table shipped in two of them. The same held-out check was
written twice. A gain quoted in one repo and a gate quoted in another could not be
read against each other, because there was no single place where the provenance of
the bytes and the number that came out of them met. Here they meet: one ledger,
one hash chain, one report.

    python -m mspine tables      # what ships, with coverage and provenance per table
    python -m mspine bench       # run the campaigns, the laws and the mechanism
    python -m mspine report      # print what the ledger holds, and nothing else
    python -m mspine verify      # re-derive every pinned digest and the chain
    python -m pytest tests -q    # the invariants

## The four layers

    data        four+ measured tables, SHA-256 pinned, verified on load
    campaign    budgeted screening under a metered API, held-out scoring
    laws        engineering laws tested against controls that could kill them
    mechanism   reaction-network reduction, with a gate that can only reject
    ledger      append-only hash chain; every number this package emits

Anything not on the chain does not appear in the report. The report has no
memory of its own, which is the whole difference between a result and a story.

## The fuel, and why coverage is printed first

    table        rows   grid    cells  cover   target
    buchwald     4599   4608     4599  99.80%  isolated yield (%)
    ccpp         9568  65536     4367   6.66%  net hourly output
    gasturbine  36733 390625     3381   0.87%  NOx emission
    concrete     1030  65536      269   0.41%  28-day strength (MPa)

Coverage is measured cells over the product grid. It is the honest measure of how
much of a design space a table explores, and it ranges over three orders of
magnitude -- which means a campaign's difficulty is not the same thing on two
rows of that table, and any comparison that does not print it is comparing
apples to a plate.

## What the screening actually found

64 experiments, batches of 8, 3 seeds, held-out gain (mean measured response of
the top decile the surrogate picked, minus the held-out mean) with the permuted
control beside it:

    buchwald   ofat +28.6   random +27.9   fill +23.9   ucb +20.6   ei +19.7
    ccpp       fill +28.7   pi +27.5      ofat +27.4   ei +26.9   random +26.3
    concrete   fill +21.5   random +18.9   ucb +16.9    ofat +16.3  ei +12.4

Every control collapses to about zero, which is the part that makes the numbers
worth reading: the gain is coming from the measured outcomes, not from the
features' correlation with themselves.

The uncomfortable part, kept in the table rather than in a footnote: on none of
the three tables does a Bayesian acquisition beat the two baselines a chemist
would use without a model. One-factor-at-a-time wins Buchwald. Greedy
space-filling wins the other two. EI is last or second-last on all three. If the
project had reported only its best arm it would have reported the opposite of
its own result.

## The laws, and the one that failed

Hydration (concrete, 136 mixes, leave-one-age-out, anchored on each mix's own
S(28)): the bounded CEB-FIP form halves the error of a power law -- 3.45 vs 6.85
MPa mean absolute, and at 90 days 2.39 vs 10.61. Its one free parameter k is
carried by composition (held-out R2 +0.25, permuted control at +0.05, p=0.000)
and *not* by strength level (R2 -0.06, p=0.32). Both claims run against
permuted-target controls.

Weakest-link (carbon fibre, 255 filaments, 3 gauge lengths): the law requires one
Weibull modulus to govern both the scatter within a length and the scale shift
across lengths. They disagree by a factor 1.55 -- and the shuffle null, which
contains no length effect by construction, produces a *larger* ratio in 99% of
its draws. So this design cannot support the claim either way, and the reusable
product is the instrument plus that null, not a verdict.

## The mechanism, and the honesty about its window

Four stages over a 9-species / 14-reaction hydrogen mechanism: propose
reductions, keep only those an equivalence gate passes, re-prove an invariant on
the reduction, break a reaction and require the gate to say no, then measure what
moved.

    t_end     active rxns   H2 converted   reductions accepted
    1e-03          2            0.0%           3 of 4
    1e-02          2            0.1%           3 of 4
    1e-01          4            1.1%           0 of 4
    1e+00          7            8.5%           0 of 4
    1e+01          9           24.9%           0 of 4

That table is the finding. Agreement between a reduced and a full mechanism is a
property of an *observation horizon*, not of the mechanism: over a window short
enough that the network has barely moved, every reduction passes because there is
nothing to see. At the ignited window the gate refuses all four candidates and
the run reports NO_WIN -- and the one candidate that is accepted at 1 ms drops a
step carrying 8e-08 of the flux, so it is labelled vacuous. A report that says
"the reduction is equivalent" without naming the window has not said anything.

## What is deliberately not here

No experimental validation, and no claim to have discovered a material. The
tables are finished data; nothing here predicts a synthesis nobody has run. What
it can do is say, with a measured held-out number and a control behind it,
whether a loop ranks untried recipes better than chance -- a much smaller claim,
and the only one the data supports.

The mechanism module is one mechanism and four hand-written reduction
candidates. It is a measuring instrument, not a mechanism-reduction algorithm.

## Reproducing

    python -m mspine bench --tables buchwald,ccpp,concrete --seeds 3 --budget 64
    python -m mspine report

Every number above is a row on `results/ledger.jsonl`, each hashed with its
predecessor, each carrying the SHA-256 of the table it was computed from. The
report reads the chain and prints the last row per (table, arm, control), so a
re-run adds rows instead of overwriting them -- the history is the audit.

## Provenance and licences

This repository merges three MIT/Apache-2.0 projects of my own; the licences are
carried, not dropped (see `mspine/tables/PROVENANCE.md` and
`LICENSE-kinetics-spine`). Every dataset is redistributed under its upstream
licence with its digest pinned, and the digests are checked on load: a mismatch
stops the run rather than producing a plausible number from unknown bytes.
