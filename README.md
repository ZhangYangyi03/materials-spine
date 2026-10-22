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
    python -m mspine bench --refits 0,8   # both fitting policies, side by side
    python -m mspine report      # print what the ledger holds, and nothing else
    python -m mspine verify      # re-derive every pinned digest and the chain,
                                 # and check the head against results/anchor.json
    python -m mspine anchor      # publish the head digest outside the ledger
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
    gasturbine  36733 390625     3381   0.87%  NOx emission (mg/m3)
    concrete     1030  65536      269   0.41%  28-day strength (MPa)

Coverage is measured cells over the product grid. It is the honest measure of how
much of a design space a table explores, and it ranges over three orders of
magnitude -- which means a campaign's difficulty is not the same thing on two
rows of that table, and any comparison that does not print it is comparing
apples to a plate.

## What the screening actually found

64 experiments per run, batches of 8, five seeds, held-out gain (mean measured
response of the top decile the surrogate ranked among cells the campaign never
spent, minus the mean of all unspent cells), with the permuted control beside it:

    table      arm        fit   gain       sd    rank   top10   control   %of range
    buchwald   random     0     +29.12    4.15  +0.581   0.02    +0.22      29.1%
               ofat       0     +27.61    4.34  +0.509   0.06    +0.23      27.6%
               fill       0     +25.38    4.06  +0.559   0.00    +0.92      25.4%
               qei        8     +26.06    8.94  +0.295   0.06    +0.39      26.1%
               qei        0     +22.20    7.10  +0.341   0.04    +0.89      22.2%
               ei         0     +21.97    4.71  +0.352   0.04    +1.27      22.0%
    ccpp       fill       0     +28.76    0.11  +0.958   0.28    -0.53      38.1%
               pi         0     +27.61    0.71  +0.826   0.18    +0.05      36.6%
               qei        0     +27.06    1.58  +0.803   0.22    -0.07      35.8%
               ei         0     +26.39    1.93  +0.771   0.20    -0.13      35.0%
               random     0     +25.38    1.93  +0.926   0.02    +0.21      33.6%
    gasturbine fill       0     +24.63    0.46  +0.719   0.06    -0.02      26.2%
               random     0     +19.96    5.48  +0.677   0.00    -0.10      21.2%
               ofat       0     +16.47    3.57  +0.512   0.02    +0.18      17.5%
               qei        0     +15.44    4.78  +0.211   0.02    -0.14      16.4%
               ei         0     +13.61    5.60  +0.163   0.02    +0.27      14.5%
    concrete   fill       0     +22.24    1.35  +0.797   0.32    +0.10      27.7%
               random     0     +20.23    3.18  +0.769   0.32    +0.06      25.2%
               ucb        0     +15.95    4.24  +0.364   0.28    +0.23      19.9%
               pi         0     +12.85    3.58  +0.351   0.32    -0.39      16.0%
               ei         0     +12.67    3.26  +0.350   0.32    -0.21      15.8%

(abridged; `python -m mspine report` prints every arm on the chain)

The controls are small and they are not zero -- the largest reads +1.27, about 1.3%
of the yield range, against real gains of 22-29 units. Read that column as a scale
and not as a threshold: with five seeds, a control of +1.27 with a standard
deviation of 1.9 is indistinguishable from zero, while a control of -0.53 with a
standard deviation of 0.2 is a small real negative, which is what a surrogate
fitted to noise should give you.

The uncomfortable part, kept in the table rather than in a footnote: on none of
the four tables does any Bayesian acquisition beat the best baseline a chemist
would use without a model. On Buchwald the ordering is blunt -- random search
+29.1, one-factor-at-a-time +27.6, space-filling +25.4, and all eight Bayesian
arms below all three. On gasturbine space-filling wins by 4.7 units over random
search with a standard deviation of 0.46, which is the tightest result in the
whole study and the least ambiguous: there is no seed noise left to explain it
away. That is not a tuning problem a better acquisition would fix; eight were run.
If the project had reported only its best arm it would have reported the opposite
of its own result.

**And the number depends on how the surrogate is fitted, which is why that is a
column.** `qei` reaches +26.1 on Buchwald only when the marginal-likelihood fit is
re-run after every batch; with the default (fit once on the initial design, then
condition) it gets +22.2. Both policies are rows on the chain, because a single
number without the policy is not a result. The per-seed table in
[docs/paper.md](docs/paper.md) shows what that +3.9 actually is: one seed out of
five recovering from a stall by +30, against three of four seeds moving slightly
the wrong way. The effect is real on the mean and its spread is wider than the
effect, and both facts are in the paper.

## The laws, and the one that failed

Hydration (concrete: 427 mixes, 416 of them with a measured 28-day strength, 529
leave-one-age-out rows, each anchored on its own mix's S(28)). The bounded CEB-FIP
form roughly halves the error of a power law:

    form      mean |err|   median    p90      n
    ceb         3.45 MPa   2.51     7.51    529     closer than power on 339 of 529
    power       6.85 MPa   4.62    15.91    529
    aci        15.96 MPa   8.05    24.00    420

and it holds at the ages where the two forms diverge: at 90 days, 2.39 vs 10.61
MPa; at 365 days, 2.73 vs 4.67. Its one free parameter k is carried by composition
(held-out R2 +0.25, permuted-target 95th percentile +0.05, p=0.000) and not by
strength level (R2 -0.06, p=0.32). Both claims run against permuted-target
controls, so "the parameter is a material property" is a measured statement rather
than a fit that happened to work.

Weakest-link (carbon fibre: 255 single filaments at 3 gauge lengths). The law
requires one Weibull modulus to govern both the scatter within a length and the
scale shift across lengths. They disagree by a factor of 1.55 -- scatter modulus
6.58, scaling modulus 10.19 -- and the pool-and-redistribute null, which by
construction contains no length effect at all, produces a ratio at or above the
observed one in 99.3% of 153 draws (median 37.4). So this design cannot support the
claim in either direction, and the reusable product is the instrument plus that
null, not a verdict.

## The mechanism, and the honesty about its window

Four stages over a 9-species / 14-reaction hydrogen mechanism: propose reductions,
keep only those an equivalence gate passes, re-prove an invariant on the reduction,
break a reaction and require the gate to say no, then measure what moved.

    t_end     active rxns   H2 converted   reductions accepted
    1e-03          2            0.0%           3 of 4
    1e-02          2            0.1%           3 of 4
    1e-01          4            1.1%           0 of 4
    1e+00          7            8.5%           0 of 4
    1e+01          9           24.9%           0 of 4

That table is the finding. Agreement between a reduced and a full mechanism is a
property of an observation horizon, not of the mechanism: over a window short
enough that the network has barely moved, every reduction passes because there is
nothing to see. At the ignited window the gate refuses all four candidates and the
run reports NO_WIN. The one candidate that is accepted at 1 ms drops a step
carrying 8e-08 of the flux, so it is labelled vacuous -- a pass that says the test
could not look, not that the reduction is safe.

Two more numbers are reported because the gate cannot be finer than its own
integrator: at the configured 400 RK4 steps the solver's self-error against a 16x
refined run is 9e-4, three orders coarser than the gate's 1e-6 tolerance, and the
second integrator (LSODA) is compared at 16x for the same reason (agreement 5.9e-15
there, 4.8e-4 at the coarse step). A report that says "equivalent to 1e-6" without
that floor next to it has claimed a resolution the numerics do not have.

## What is deliberately not here

No experimental validation, and no claim to have discovered a material. The
tables are finished data; nothing here predicts a synthesis nobody has run. What
it can do is say, with a measured held-out number and a control behind it,
whether a loop ranks untried recipes better than chance -- a much smaller claim,
and the only one the data supports.

The mechanism module is one mechanism and four hand-written reduction
candidates. It is a measuring instrument, not a mechanism-reduction algorithm.

## Reproducing

    python -m mspine bench --tables buchwald,ccpp,concrete --seeds 5 --budget 64 \
        --refits 0,8
    python -m mspine report

Every number above is a row on `results/ledger.jsonl`, each hashed with its
predecessor, each carrying the SHA-256 of the table it was computed from. The
report reads the chain and prints the last row per (table, arm, policy, control),
so a re-run adds rows instead of overwriting them -- the history is the audit.

Two properties of that chain are worth stating plainly, because both were found by
trying to break it rather than by design:

  * it detects an edit to any interior row, and it cannot detect an edit to the
    last row -- that is a property of hash chains, not a defect, and the answer is
    `Ledger.anchor()`, which publishes the head digest somewhere the ledger cannot
    reach (a commit message will do). `verify(anchor_head=...)` checks it.
  * it needed a canonical form for keys, not just for values. `json.dumps` writes
    an integer key as `"56"` in numeric order and reloads it in lexical order, so
    a `{56: ..., 180: ...}` payload produced bytes that did not match the bytes
    that were hashed, and an honest ledger reported itself as tampered. Two tests
    now hold that down (`tests/test_ledger.py`).

The campaigns themselves are run in pieces -- eight arms x five seeds x three
tables is more than one process on this host is allowed to hold -- and
`bench.summary_from_ledger()` rebuilds the whole table from the chain, so a
benchmark split across ten restarts is still one result. The part files are a
convenience; the chain is the record.

## Provenance and licences

This repository merges three MIT/Apache-2.0 projects of my own; the licences are
carried, not dropped (see `mspine/tables/PROVENANCE.md` and
`LICENSE-kinetics-spine`). Every dataset is redistributed under its upstream
licence with its digest pinned, and the digests are checked on load: a mismatch
stops the run rather than producing a plausible number from unknown bytes.
