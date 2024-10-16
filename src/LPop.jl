module LPop

    using Statistics

    include("pso.jl")
    include("aux_gene.jl")
    include("genetic.jl")

    export PSO, Genetic

end
