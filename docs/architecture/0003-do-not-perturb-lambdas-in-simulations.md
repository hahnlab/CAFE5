# 3. Do not perturb lambdas when simulating families

Date: 2020-04-14

## Status

Accepted

## Context

Simulations are generated based on the lambda value provided by the user. It seemed that the simulations
we generated were not realistic enough. It was thought that slightly modifying the lambda value every
so often would create more realistic simulations, and this code was implemented.

## Decision

We will simply base all simulated families on the given lambda and not perturb the lambda.

## Consequences

The simulated families should have a slightly smaller variance in sizes. If one wants to recreate simulated 
conditions, I think they have to not use perturb, but if one wants to simulate realistic conditions and 
recover them with the same number of rate categories the perturb works better (although I made the sttdev on the normal 
distribution around the multipliers much smaller and it got rid of the negative lambdas and still worked well 
(dividing by like 50 instead of 5).




