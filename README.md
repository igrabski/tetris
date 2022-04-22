# Tetris
Tetris performs Bayesian combinatorial multi-study factor analysis (see our preprint [here](https://arxiv.org/pdf/2007.12616.pdf)) in R.

# Usage
To run Tetris, you should first prepare your data as a list of matrices, where each matrix corresponds to data from a different study/group and is in samples by features format. You can then run the sampler with the following command:

```
# X is a list of matrices
# alpha, beta are hyperparameters 
run <- tetris(X,alpha=5,beta=1)
```

To recover point estimates of the factor indicator matrix and loadings matrix, you can run:

```
# alpha_IBP is the hyperparameter; S is the number of studies
A <- choose.A(run,alpha_IBP=5,S=3)
run_fixed <- tetris(X,alpha=5,beta=1,fixed=TRUE,A.fixed=A)
Lambda <- getLambda(run_fixed,A)
```

Tetris also has a clustering extension, which simultaneously estimates group structure alongside factor analysis parameters. To run this extension, you should prepare your data as a single matrix in samples by features format, and run:

```
# X is a matrix
# S is the number of clusters you want to find
# alpha, beta are hyperparameters
run_clustering <- tetris_clustering(X,S=3,alpha=5,beta=1)
```
