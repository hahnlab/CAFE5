# 4. Use median rather than mean value when discretizing gamma values

Date: 2020-05-28

## Status

Accepted

## Context

In order to simulate a smooth gamma curve, we discretize values from the gamma curve and
set categories based on them.

## Decision

It was decided to use the median value rather than the mean for calculating the discrete values.

## Consequences

Slightly more accurate.



