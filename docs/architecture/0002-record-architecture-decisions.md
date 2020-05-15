# 2. Calculate a prior of root distributions based on user input

Date: 2020-05-15

## Status

Accepted

## Context

A prior distribution needs to be specified to accurately calculate probabilities. The prior can
be calculated in many different ways. Hopefully this decision reflects the least surprising
results for the user.

## Decision

The prior will be calculated as follows:

* if -p specified on command line: calculate a Poisson distribution with the specified lambda
* if -f specified on command line: The user has specified a root distribution. Use that. Issue a warning if the user has ALSO specified a Poisson lambda
* if -i specified on command line: Estimate a Poisson distribution from the families provided.
* Otherwise, use a uniform distribution. Issue a warning as this is not a very reasonable prior.

## Consequences

Some users may be confused by the ordering. There may be some unforeseen circumstance where the combination of input flags is still surprising.


