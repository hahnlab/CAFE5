# 6. The calculation of P-values

Date: 2020-11-12

## Status

Accepted

## Context

We want to provide users with the pvalue of each individual family under the lambda provided (either
calculated or provided by the user)

## Decision

The procedure is as follows: 

1. For every possible root family size, randomly generate 1000 families based on the lambda
2. Compute the probability of each family generated, sort the results smallest to largest. The result is the
   conditional distribution.
3. Compute the probability of each of the user's families and calculate the pvalue from this likelihood and the 
   conditional distribution

## Consequences

Note that although the family generation routine can accept an error model, we do not use the model when 
calculating p-values. It is not clear if using a lambda estimated with the error model will produce the 
correct conditional distribution at the tips.

The error model attempts to account for error in the lambda due to poor annotation, genome, etc. Once this has 
been accounted for in the rate estimation, if this same rate is used to calculate probabilities, we should
generate an expected conditional distribution that reflects the "true", corrected, distribution. (However, when 
families are generated in the simulator, we take an error model into account by adjusting family sizes at the tip,
adding or subtracting a single copy based on the model. This is at the discretion of the user. ) 