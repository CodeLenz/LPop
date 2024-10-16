# LPop

Optimization using Populational Methods (PSO and Genetics)

<p align="center">
<img src="./docs/population.png" alt="LPop" width="500">
</p>

## Install
```julia
using Pkg
Pkg.add("url=https://github.com/CodeLenz/LPop")
```

## Examples

First, define a function of $x \in \mathbb{R}^n$ only, such that $f:\mathbb{R}^n \to \mathbb{R}$

```julia
f(x) = (x[1]-3)^2 + (x[2]-5)^2
```
with minimum at $x_1=3$ and $x_2=5$.


The solution using PSO is 
```julia
# Number of particles
Np = 20

# Number of iterations (epochs)
niter = 100

# Lower bounds
xl = zeros(2)

# Upper bounds
xu = 10*ones(2)

# Call the optimizer
x_opt = PSO(f,Np,niter,xl,xu)

2-element Vector{Float64}:
 2.999999999994289
 5.000000000002956

```

Solution using Genetic Algorithms
```julia
# Number of individuals
Np = 20

# Number of iterations (epochs)
ngera = 100

# Lower bounds
xl = zeros(2)

# Upper bounds
xu = 10*ones(2)

# Number of bits
Nb = 6

# Mutation rate
mutrate = 2/100

# Type of replacement
replacement = "Elite"

# Call the optimizer
x_opt = Genetic(f,Np,ngera,Nb,mutrate,xl,xu,replacement)


```
