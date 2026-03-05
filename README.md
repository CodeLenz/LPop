# LPop

Unconstrained optimization using Populational Methods (PSO and Genetics)

For constrained optimization there is a PSO implementation (see the examples)

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

# Number of bits (per variable)
Nb = 6

# Mutation rate
mutrate = 2/100

# Type of replacement
replacement = "Elite"

# Call the optimizer
x_opt = Genetic(f,Np,ngera,Nb,mutrate,xl,xu,replacement)

```

## Constrained problems using PSO. 

The approach uses the procedure proposed by Deb (Comput. Methods Appl. Mech. Engrg. 186 (2000) 311±338). Lets solve the first example in this manuscript (converting the inequalities to <=0).
This is a hard problem, since the feasible region is only $0.7\%$ of the design space.

```julia

# Objective function and constraints (<=0)
function f(x) 
        obj = (x[1]^2+x[2]-11)^2 + (x[1] + x[2]^2 -7)^2 
        g1 = -(4.84 - (x[1]-0.05)^2 -(x[2]-2.5)^2)
        g2 = -(x[1]^2 + (x[2]-2.5)^2 -4.84 )
        return [obj; g1; g2]
end

# Lower and upper side constraints
xl = [0.0 ; 0.0]
xu = [6.0 ; 6.0]

# Number of particles
NP = 100

# Number of iterations 
niter = 100

# Solve the problem
x_opt, fg_opt = PSO_Deb(f,NP,niter,xl,xu)


```

The optimal constrained solution is $\mathbf{x} \approx (2.24, 2.38)$ with $f\approx 13.59$. The first constraint is active, $g_1 \approx 0$ and the second constraint is inactive, $g_2 \approx -0.22$.

Using the standart coefficients (no tunning to this problem) and running the program 1000 times with 20 particles and 100 iterations gives 753 correct runs (75\%). Using 100 particles gives 1000 correct runs (100\%)


# TODO

Paralelization

Other Selection and Replacement operators for GA