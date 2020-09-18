# 3. Calculate a prior of root distributions based on user input

Date: 2020-09-18

## Status

Accepted

## Context

Simulations are generated based on the lambda value provided by the user. It seemed that the simulations
we generated were not realistic enough. It was thought that slightly modifying the lambda value every
so often would create more realistic simulations, and this code was implemented.

## Decision

The effect on the simulated families doesn't seem to be large. We will simply base all simulated families
on the given lambda.

## Consequences

The simulated families should have a slightly smaller variance in sizes.


