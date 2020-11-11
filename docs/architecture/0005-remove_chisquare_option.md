# 5. Remove chi-squared option

Date: 2020-11-11

## Status

Accepted

## Context

The original intent was to run both the base and gamma models and compare them with a chi-squared test.

## Decision

Currently it seems more effective to allow the user to run base and gamma models separately and compare them
using their own favorite method.

## Consequences

Remove -r option from the code as well as chi-square classes.
