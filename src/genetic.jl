"""
Genetic algorithm (classic method)

     Genetic(f::Function,Np::Int,ngera::Int,Nb::Int,mutrate::Float,xl::Vector,xu::Vector,replacement;verbose)

Inputs: 

     f(x):        Objective function. Must return a single scalar 

     Np:          Number of individuals 

     ngera:       Number of generations (epochs) 

     Nb:          Number of bits 

     mutrate:     Mutation rate

     xl:          Vector with lower bounds (side constraints) 

     xu:          Vector with upper bounds (side constraints) 

     replacement: type of replacement (Full, Elite)

     verbose: show the best result of each iteration (if true)

Outputs:

     vector with the best solution

"""
function Genetic(f::Function, Np::Int64, ngera::Int64,Nb::Int64,
                 mutrate::Float64, xl::Vector, xu::Vector,
                 replacement="Elite"; verbose=false)

    # Verifica se as opções de replacement são válidas
    @assert replacement in ["Full","Elite"] "Geneticos::tipo de replacement não implementado"

    # Garante que as dimensões são iguais
    @assert length(xl)==length(xu)

    # Número de variáveis de projeto
    n = length(xl)    

    # Gera a população inicial
    # Matriz com Np × (n * Nb)
    P = rand(Bool,Np,n*Nb)

    # Calcula o vetor de variáveis de projeto
    # para essa população
    valores =  Calcula_Objetivos(P,f,n,Nb,xl,xu)

    # Loop das gerações do algoritmo genético
    for gera=1:ngera

        # Agora Vamos pegar a melhor solução 
        val_opt, pos_opt = findmin(valores)

        # Melhor solução é 
        x_opt = Extrai_Individuo(P[pos_opt,:],n,Nb,xl,xu)

        # Mostra a melhor solução dessa geração
        verbose &&  println("Generation ", gera, " best ",val_opt)
        
        # Faz o torneio simples
        papis, mamis = Torneio_Simples(valores,Np)

        # Crossover de um ponto
        Pbar = CrossOver(papis,mamis,P,n,Nb,Np)

        # Mutação
        Phat = Mutacao(Pbar,Np,n,Nb,mutrate)

        # Calcula os objetivos da nova população
        valores_novos = Calcula_Objetivos(Phat,f,n,Nb,xl,xu)

        # Replacement
        Replacement!(replacement,Np,P,Phat,valores,valores_novos)


    end #gera

    # Agora Vamos pegar a melhor solução 
    val_opt, pos_opt = findmin(valores)

    # Melhor solução é 
    x_opt = Extrai_Individuo(P[pos_opt,:],n,Nb,xl,xu)

    # Mostra a melhor solução
    println("Solução é ", val_opt)
    print(x_opt)

    return x_opt
               
end 
