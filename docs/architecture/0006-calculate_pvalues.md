# 6. The calculation of P-values

Date: 2020-11-12

## Status

Accepted

## Context

We want to provide users with the pvalue of each individual family under the lambda provided (either
calculated or provided by the user)

## Decision

The procedure is as follows: 

1. For every possible root family size (starting with 1, not 0), randomly generate 1000 families based on the lambda
2. Compute the likelihood of each family generated, sort the results smallest to largest. The result is the
   conditional distribution.
3. Compute the likelihood of each of the user's families for every root family size 
4. For every root family size less than 125% of the largest species size in the family, compute a pvalue for each 
   of the user's families, based on the conditional distribution at that family size
5. Take the maximum pvalue calculated

## Consequences

To correctly calculate the likelihood of the generated families, the calculation at the root has to take
into account that the root size is known. This is reflected in the matrix multiplication, which uses a 
single row at the root rather than all possible rows.

Note that although the family generation routine can accept an error model, we do not use the model when 
calculating p-values. It is not clear if using a lambda estimated with the error model will produce the 
correct conditional distribution at the tips.

The error model attempts to account for error in the lambda due to poor annotation, genome, etc. Once this has 
been accounted for in the rate estimation, if this same rate is used to calculate probabilities, we should
generate an expected conditional distribution that reflects the "true", corrected, distribution. (However, when 
families are generated in the simulator, we take an error model into account by adjusting family sizes at the tip,
adding or subtracting a single copy based on the model. This is at the discretion of the user. ) 
