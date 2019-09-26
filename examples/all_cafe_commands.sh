# Estimating a single lambda for the whole tree

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -o singlelambda

# Estimating a single lambda for the whole tree using a Poisson distribution for the root frequency

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -o singlelambda_poisson

# Estimating a separate lambda for the chimp/human branch of the tree

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -y chimphuman_separate_lambda.txt -o doublelambda

# Reconstruct a phylogeny using a given lambda and taking into account an error model

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -l 0.01 -e errormodel_0.1.txt -o errormodel

# Estimate a lambda and error model

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -e errormodel_0.1.txt -o lambda_epsilon

# Estimate a lambda along with a gamma distribution using three rate categories

../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -k 3 -o gamma_dist

# Reconstruct a phylogeny with two lambda values 

./bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -m 0.01,0.05 -y chimphuman_separate_lambda.txt -o lambdas01_05

# Reconstruct a phylogeny with a given lambda and shaped gamma distribution
../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -l 0.002 -k 3 -a 0.425 -o alpha425

# Simulate 100 families with randomly selected root sizes
../bin/cafexp -s100 -l 0.002 -t mammals_tree.txt -o sim100

# Simulate 1000 families with a Poisson distribution of root sizes
 ../bin/cafexp -s -f poisson_root_dist_1000.txt -l 0.002 -t mammals_tree.txt -o simpoisson1000

# Simulate 1000 families with a shaped gamma distribution
../bin/cafexp -s1000 -l 0.002 -k 4 -a .4 -t mammals_tree.txt -o simalpha4 

# Estimate a separate lambda for each family
../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -b -o lambdaperfamily

