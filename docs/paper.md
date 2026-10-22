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

    metric    best-so-far / ceiling      rows the campaign itself chose to run (n=5)
              held-out screening gain    rows it never saw (n=5)

    fit       refit=0   hyperparameters fitted once on the initial design, then
                        held while the posterior is conditioned
              refit=8   the marginal-likelihood fit re-run after every batch

`process-spine` used best-so-far with refit=8. This package defaulted to held-out
gain with refit=0. So the disagreement was not about the data or the seeds.

## 2. The fitting policy, isolated

Same table, same 64-experiment budget, batches of 8, same 8-point initial design,
five seeds, one variable. Held-out gain, which is the metric this package reports
for the reasons in section 4:

    buchwald, held-out gain, 5 seeds
      arm      refit=0              refit=8
      qei      +22.20 +- 7.10       +26.06 +- 8.94
      ei       +21.97 +- 4.71       not re-run
      fill     +25.38 +- 4.06       no surrogate: identical by construction
      ofat     +27.61 +- 4.34       identical by construction
      random   +29.12 +- 4.15       identical by construction

Re-tuning is worth +3.9 to `qei` on this table, and the three arms without a
surrogate do not move by a digit -- which is the check that the variable is the
surrogate and not the seed. (`ei` at refit=8 was not re-run at five seeds; the
three-seed study that did run it is on the chain as `buchwald/ei/refit8/real`.) The first version of this section quoted a +0.047 on
best-so-far/ceiling, reproducing `process-spine`'s published number at refit=8;
that metric is the one section 4 argues against, and the effect is +3.9 units of
yield on the metric that asks about rows the campaign did not choose.

## 3. And then most of it stops being a result

Per seed, held-out gain on Buchwald:

    arm   refit   seed0    seed1    seed2    seed3    seed4     mean
    qei   r0     +24.52   +22.69   +10.42   +29.59   +23.80    +22.20
    qei   r8     +22.66   +19.84   +40.96   +19.32   +27.53    +26.06

Seed 2 is the whole difference: +10.42 without re-tuning, +40.96 with it, a
+30.5 swing on that seed alone, contributing +6.1 to a five-seed mean. On the
other four seeds re-tuning is -1.9, -2.9, -10.3 and +3.7 -- three of four worse.
The mean is positive because one seed recovered from a stall by a large margin,
and the standard deviation is 7-9 units, which is larger than the effect.

So "re-tuning is what makes batch EI work" is not supported, and neither is
"re-tuning is a dial to be set to whatever wins". What is supported, and it is
narrower: on Buchwald at budget 64, batch 8, refit=0 lets batch EI stall on at
least one seed, and refit=8 unsticks it at the cost of variance elsewhere.
Different claims; the second is worth knowing and the first was what both
repositories' prose implied.

## 4. Which metric is the measurement

Every `best-so-far` row above is measured on rows the campaign chose. That number
cannot separate a method that finds good recipes from one that merely looked at
the ones it found, and it is the number both old papers quoted. The held-out gain
asks about the cells the campaign did not spend:

    fit on what the campaign spent, rank the unspent cells, take the top decile,
    compare their measured mean with the measured mean of all unspent cells

Under that metric, five seeds, the merged table (`python -m mspine report`, all
rows on the chain):

    table      best arm            best Bayesian arm       permuted control
    buchwald   random +29.12       qei +26.06 (refit=8)     +0.05 .. +1.27
    ccpp       fill +28.76         pi +27.61                -0.53 .. +0.21
    concrete   fill +22.24         pi +12.85                -0.39 .. +0.45

No Bayesian acquisition beats the best model-free baseline on any of the three
tables, under either fitting policy. On Buchwald the ordering is blunt: random
search +29.1, one-factor-at-a-time +27.6, greedy space-filling +25.4, and every
Bayesian arm below all three. That is not a tuning problem to be fixed by choosing
a better acquisition -- eight arms were run.

The controls are small and not zero: the largest is +1.27, about 1.3% of the yield
range, against real gains of 22-29 units. Read the control column as a scale
rather than a threshold. With five seeds a control row that reads +1.27 with a
standard deviation of 1.9 is not distinguishable from zero; a row that reads -0.53
with a standard deviation of 0.2 is, and it is a small negative, which is what a
surrogate fitted to noise should produce.

On concrete, random search is second best on a grid its own table covers 0.41% of,
which says more about the table than about random search: with 269 measured cells
out of 65536, a surrogate has very little to learn from and a space-filling
baseline is already near-optimal.

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

Five seeds for every per-table mean, which is enough to show the ordering between
random search and every Bayesian arm is stable on Buchwald and not enough to
resolve a 3.9-unit effect whose per-seed spread is 7-9. One table decides the
refit question, and it is decided by one seed. A sweep over budgets
and batch sizes would say whether the stall is a property of refit=0 at batch 8
or a property of seed 2. It has not been run, and until it is, the honest
statement is the narrow one: on Buchwald at budget 64, batch 8, these five seeds,
refit=8 removes one stall. Nothing here licenses a general claim about re-tuning.
