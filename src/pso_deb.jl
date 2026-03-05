#
# Rotina principal. 
#
#  min f(x)
#
#  T.q
#         g_j(x) <= 0 
#      xl_i <= x_i <= xu_i
#
#
# As entradas são
#
# f     -> função que retorna um vetor com função objetivo e restrições de desigualdade (deve depender somente de x) 
# NP    -> número de partículas
# niter -> número de iterações do método PSO
# xl    -> restrições laterais inferiores
# xu    -> restrições laterais superiores
# w     -> fator de inércia
# c1    -> coeficiente para o termo local
# c2    -> coeficiente para o termo global
# verbose -> true para mostrar a evolução da otimização
#
"""
Particle swarm optimization (classic method)

     PSO_Deb(f::Function, NP::Int,niter::Int,xl::Vector,xu::Vector;w=0.5,c1=1.5,c2=1.5,verbose=false)

Inputs: 

     f(x):    Vector with objective function and constraints

     NP:      Number of particles in the swarm 

     niter:   Number of iterations 

     xl:      Vector with lower bounds (side constraints) 

     xu:      Vector with upper bounds (side constraints) 

     w:       Weight factor (inertia) 

     c1:      Coefficient for local direction 

     c2:      Coefficient for global direction 

     verbose: show the best result of each iteration (if true)
     

Outputs:

     x_global_best: vector with the best solution
     s1: objective and constrains at the best solution

"""
function PSO_Deb(f::Function, NP::Int,niter::Int,xl::Vector,xu::Vector;w=0.5,c1=1.5,c2=1.5,verbose=false)

    # Verificações de consistência
    length(xl)==length(xu) || error("PSO_Deb:: dimensão dos vetores xl e xu deve ser igual")

    # Número de variáveis de projeto
    n = length(xl)

    # Verifica se o número de partículas é positivo
    NP > 1 || error("PSO_Deb:: número de partículas deve ser maior do que 1")
    
    # Verifica se temos ao menos NP = 10*n 
    NP >= 10*n || println("PSO_Deb:: sugerimos utilizar NP>10n se possível. Continuando...")

    # Verifica se temos niter > 1
    niter > 1 || error("PSO_Deb:: número de iterações deve ser maior do que um")

    #
    # Inicializa o array de partículas, respeitando as restrições laterais
    #
    X = zeros(NP,n)
    for p=1:NP
        X[p,:] .= xl .+ (xu.-xl).*rand(n)
    end
 
    #
    # Inicializa o array de velocidades, utilizando uma distribuição normal 
    # com média nula.
    #
    V = zeros(NP,n)
    for p=1:NP
    	V[p,:] = randn(n)
    end
 
    # Chama f(x) para uma partícula, só para sabermos a dimensão 
    s1 = f(vec(X[1,:]))

    # Dimensão que f(x) retorna  
    m = length(s1)

    # Teste de consistência
    if m==1 
       error("PSO_Deb:: use PSO para problemas sem restrições")
    end

    # Array que contem o valor do objetivo e restrições para cada partícula
    F = zeros(NP,m)

    # Aproveita que já chamamos o f(x) para a primeira partícula
    F[1,:] .= vec(s1)

    # Melhor valor global das partículas e sua posição
    global_best = maxintfloat(1.0)
    x_global_best = zeros(n)

    # Melhor posição e objetivo de cada partícula
    personal_best = zeros(NP)
    X_personal_best = zeros(NP,n)

    # Indica se a partícula é viável ou não, calculando o valor da violação
    viola = zeros(NP)

    # Armazena os dados da primeira partícula 
    viola[1] = Violacao(F[1,:])

    # Inicializa pior viável da iteração
    pior_viavel = F[1,1]
    if viola[1]>0
        pior_viavel = -maxintfloat(1.0)
    end
        

    #
    # Inicializa F, lembrando que já calculamos para a primeira partícula
    #
    for p=2:NP

        # Objetivo e restrições para esta partícula
        F[p,:] .= vec(f(vec(X[p,:])))

        # Violação 
        viola[p] = Violacao(F[p,:])

        # Se for viável e se objetivo for o pior, guarda
        if viola[p]==0 && F[p,1]>pior_viavel
            pior_viavel = F[p,1]
        end

    end #p
 
    if verbose
       println("Número de partículas violadas: ",length(viola[viola.>0]))
       println("Pior viável ", pior_viavel)
    end

    # Todo mundo é inviável 
    if length(viola[viola.>0])==NP
        verbose && println("Não temos partículas viáveis")
        pior_viavel = maxintfloat(1.0)
    end

    # Processa os 
    objetivo = zeros(NP)

    # Se inviável, soma violação ao pior objetivo viável. Do contrário, soma ao 
    # valor do objetivo
    for p=1:NP

        if viola[p]>0 
           objetivo[p] = pior_viavel + viola[p]
        else
            objetivo[p] = F[p,1]
        end 

    end

    # Inicializa global_best, x_global_best, personal_best e X_personal_best
    for p=1:NP

          # Aqui o objetivo é o melhor da partícula
          personal_best[p] = objetivo[p]
          X_personal_best[p,:] .= vec(X[p,:])

          # Se também for o melhor de todos, armazena
          if objetivo[p] < global_best
               global_best = objetivo[p]
               x_global_best .= vec(X[p,:])
          end

    end #p


    #
    # Loop principal de otimização
    # 
    for iter=1:niter

        # Atualiza a velocidade de cada partícula baseado na lógica do PSO. Aproveitamos
        # para atualizar as posições também
        for p=1:NP
         	
         	# Atualiza a velocidade 
            V[p,:] .= w*V[p,:] .+ c1*rand(n).*(X_personal_best[p,:] .- X[p,:]) .+ 
                    c2*rand(n).*(x_global_best .- X[p,:])

            # Atualiza posição
            X[p,:] .= X[p,:] .+ V[p,:]

            # Aplica as restrições laterais 
            X[p,:] = max.(xl, min.(xu, X[p,:]))

        end #p

        # Inicializa o pior viável 
        pior_viavel = -maxintfloat(1.0)

        # Calcula objetivo e restrições para cada partícula
        for p=1:NP

            # Calcula o novo objetivo para esta partícula
            F[p,:] .=  vec(f(vec(X[p,:])))

            # Violação 
            viola[p] = Violacao(F[p,:])

            # Se for viável e se objetivo for o pior, guarda
            if viola[p]==0 && F[p,1]>pior_viavel
               pior_viavel = F[p,1]
            end

        end #p      


        # Todo mundo é inviável 
        if length(viola[viola.>0])==NP
            verbose && println("Não temos partículas viáveis")
            pior_viavel = maxintfloat(1.0)
        end

        # Se inviável, soma violação ao pior objetivo viável. Do contrário, soma ao 
        # valor do objetivo
        for p=1:NP

            if viola[p]>0 
               objetivo[p] = pior_viavel + viola[p]
            else
               objetivo[p] = F[p,1]
            end 

        end

        # Processa as informações de melhor valor para cada partícula 
        # e também para o enxame
        for p=1:NP

            # Verifica se ela melhorou
            if objetivo[p] <  personal_best[p] 
                personal_best[p]  = objetivo[p]
                X_personal_best[p,:] .= vec(X[p,:])
            end 

            # Verifica se ela é melhor do que o global best
            if objetivo[p]<global_best
                global_best = objetivo[p]
                x_global_best .= vec(X[p,:])
            end 

        end #p

        # Mostra a evolução (verbose=true)
        verbose && println("Melhor objetivo até agora $(global_best)")


    end #iter

   # Faz um teste com a mehor partícula 
   s1 = vec(f(x_global_best))

   #
   # Retorna a melhor solução e f(x_global_best)
   #
   return  x_global_best, s1

end


#
# Devolve a violação total de uma partícula, isto é, 
# a soma de todas as restrições violadas (>=0)
#
function Violacao(v::Vector)

    # Utilizamos somente as posições 2:end
    v1 = v[2:end]

    # Soma das posições violadas
    return sum(v1[v1.>0])

end