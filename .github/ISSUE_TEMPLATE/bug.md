---
name: a number is wrong
about: A figure in the README, the paper or the report does not match the chain
title: ""
labels: correctness
---

## The number

Which document, which line, what it says.

## What the chain says

    python -m mspine report
    python -m mspine verify

Paste the relevant rows. If `verify` is not OK, say so first -- that is a different
and more serious report.

## Which of the three the fix belongs in

    the code             the benchmark computed the wrong thing
    the prose            the number was right when written and the run changed
    the chain            a row is missing or was edited (this should be impossible)

Please say which, because the three have different fixes and only one of them is
editing a document.
