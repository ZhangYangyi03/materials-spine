# Security policy

## Reporting

Open a private security advisory on this repository, or email
zhangyangyi03@outlook.com. There is no bounty and no SLA; there is a person who
reads the mail.

## What this repository does and does not touch

It reads CSV files that ship inside the package and writes JSON to `results/`.
It makes no network requests at run time, holds no credentials, and executes no
data from the tables it loads -- `mspine/data.py` parses with `csv` and `numpy`,
not with `eval`, and the tables are verified by SHA-256 before they are used.

The one thing worth knowing: `results/ledger.jsonl` is data, not code, and it is
loaded with `json.loads`. A tampered ledger is detected by `mspine verify` rather
than trusted, which is the point of the whole layer.

## Supported versions

The default branch. This is a research repository with one author; there are no
backported branches.
