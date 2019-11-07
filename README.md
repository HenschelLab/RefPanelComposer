# RefPanelComposer
Composing a reference panel for genotype imputation using a pool of haplotypes (local and public)

Algorithm:

Input 
1. target study samples, phased, genotype array
2. Distance matrix between all samples, ie. target and reference(in later stages, this will be calculated automatically with a separate script)

Output:
imputed genotypes

Process:
the algorithm chooses the most suitable references using a K-NN approach:
for every target genotype array, k nearest neighbors from the pool of potential references are identified
then the union of those neighbors is formed and converted into a reference panel compatible with impute2.
Note if your query samples are too different, you might want to run this separately.

