# Two numbers, one campaign, and which of them survives

*This is the study a merge makes possible and a split makes impossible. The three
repositories this one replaces each published a number, and the numbers
disagreed. Nothing was wrong with either measurement: they were answers to
different questions, taken on the same table.*

## 1. What disagreed

`process-spine` reported that batch EI reaches 0.9896 of the Buchwald ceiling
against one-factor-at-a-time's 0.9785, and concluded that the Bayesian advantage
on that table "lives entirely in the batch rule". `formula-screen` reported, on
the same 4599 rows, that expected improvement *lost* to one-factor-at-a-time. Both
printed multi-seed means. Both were reproducible. They could not both be the
headline.

Two differences turn out to be load-bearing, and neither was written down:

    metric    best-so-far / ceiling      rows the campaign itself chose to run
              held-out screening gain    rows it never saw

    fit       refit=0   hyperparameters fitted once on the initial design, then
                        held while the posterior is conditioned
              refit=8   the marginal-likelihood fit re-run after every batch

`process-spine` used best-so-far with refit=8. This package defaulted to held-out
gain with refit=0. So the disagreement was not about the data or the seeds.

## 2. The fitting policy, isolated

Same table, same 64-experiment budget, batches of 8, same 8-point initial design,
same seeds. One variable:

    buchwald, 5 seeds, best-so-far / ceiling
      arm      refit=0             refit=8
      qei      0.9422 +- 0.101     0.9896 +- 0.010
      ei       0.9429 +- 0.078     0.9760 +- 0.033
      ofat     0.9785 +- 0.026     0.9785 +- 0.026    no surrogate: identical by construction
      random   0.9394 +- 0.061     0.9394 +- 0.061    no surrogate

Re-tuning is worth +0.047 to batch EI here, and the arms without a surrogate do
not move by a single digit, which is the check that the variable is the surrogate
and not the seed. `process-spine`'s published number reproduces exactly at
refit=8.

## 3. And then it stops being a result

The mean is carried by one seed out of five:

    buchwald, per seed, best-so-far / ceiling
      arm   seed    0       1       2       3       4
      qei   r0   0.9729  0.9969  0.7414  1.0000  1.0000
      qei   r8   0.9729  0.9969  0.9818  1.0000  0.9962
      ei    r0   0.9165  0.9969  0.8014  1.0000  1.0000
      ei    r8   0.9729  0.9969  0.9139  1.0000  0.9962

Three of five seeds are bit-identical under both policies; a fourth moves by
-0.0038. The entire +0.047 is seed 2 recovering from a stall, which alone
contributes +0.048 of the mean. "Re-tuning is what makes batch EI work" is not
supported by five seeds of which three are untouched and a fourth moves the wrong
way by a thousandth. What is supported is "one seed in five stalls without
re-tuning, and re-tuning unsticks it". Different claims; the second one is worth
knowing, and a test now fails if the pattern changes
(`tests/test_published.py::test_one_seed_carries_the_refit_gap`).

Re-tuning also has no consistent sign. On ccpp, same budget, same seeds, held-out
gain:

    ccpp, held-out gain, 3 seeds
      ei       refit=0  +26.90      refit=8  +21.96      worse
      qei      refit=0  +27.32      refit=8  +28.19      better
      fill     no surrogate, +28.74 either way
      ofat     no surrogate, +27.35 either way

So the policy is a source of variance to be reported, not a dial to be set to
whatever wins.

## 4. Which metric is the measurement

Every `best-so-far` row above is measured on rows the campaign chose. That number
cannot separate a method that finds good recipes from one that merely looked at
the ones it found, and it is the number both old papers quoted. The held-out gain
asks about the cells the campaign did not spend:

    fit on what the campaign spent, rank the unspent cells, take the top decile,
    compare their measured mean with the measured mean of all unspent cells

Under that metric, three seeds, both policies, the merged table:

    table      best arm        best Bayesian arm         permuted control
    buchwald   ofat +28.58     qei +27.82 (refit=8)      collapses to ~0
    ccpp       fill +28.74     qei +28.19 (refit=8)      collapses to ~0
    concrete   fill +21.51     ucb +16.92 (refit=0)      collapses to ~0

No Bayesian acquisition beats the best model-free baseline on any of the three
tables, under either fitting policy or metric. The controls collapse, so the gains
are real; they are just not the surrogate's wins. On concrete the ordering is
blunt: greedy space-filling +21.5, random +18.9, UCB +16.9,
one-factor-at-a-time +16.3. Random search is second best on a table 0.41% covered
by its own grid, which says more about the table than about random search.

## 5. What this changes in the repository

Three things, all code rather than prose:

  * `refit` is a field on every campaign record and part of the aggregation key,
    so one arm under two fitting policies is two rows and neither can silently
    overwrite the other (`mspine/campaign.py`);
  * `--refits 0,8` runs both policies side by side, so the cost of re-tuning is a
    number on the chain instead of an assumption in a docstring
    (`mspine/bench.py`);
  * the report prints the policy as a column, because a table of gains that hides
    it is a table that made the mistake above (`mspine/report.py`).

`process-spine`'s paper is not reproduced here and is not withdrawn there: its
numbers stand at refit=8 on best-so-far, and that repository says so. What is new
is that the two conclusions can be read against each other, and the reading is
that the interesting quantity was the one neither repository reported -- how much
of a published gain is one seed.

## 6. What this does not show

Three seeds for the per-table means, five for the per-seed table. One table
decides the refit question, and it decides it on one seed. A sweep over budgets
and batch sizes would say whether the stall is a property of refit=0 at batch 8
or a property of seed 2. It has not been run, and until it is, the honest
statement is the narrow one: on Buchwald at budget 64, batch 8, these five seeds,
refit=8 removes one stall. Nothing here licenses a general claim about re-tuning.
