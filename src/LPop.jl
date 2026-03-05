module LPop

    using Statistics

    include("pso.jl")
    include("pso_deb.jl")
    include("aux_gene.jl")
    include("genetic.jl")

    export PSO, PSO_Deb, Genetic

end
