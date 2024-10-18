#
# Rotina principal. As entradas são
#
# f     -> função objetivo (deve depende somente de x) 
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

     PSO(f::Function, NP::Int,niter::Int,xl::Vector,xu::Vector;w=0.5,c1=1.5,c2=1.5,verbose=false)

Inputs: 

     f(x):    Objective function. Must return a single scalar 

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

"""
function PSO(f::Function, NP::Int,niter::Int,xl::Vector,xu::Vector;w=0.5,c1=1.5,c2=1.5,verbose=false)

    # Verificações de consistência
    length(xl)==length(xu) || error("PSO:: dimensão dos vetores xl e xu deve ser igual")

    # Número de variáveis de projeto
    n = length(xl)

    # Verifica se o número de partículas é positivo
    NP > 1 || error("PSO:: número de partículas deve ser maior do que 1")
    
    # Verifica se temos ao menos NP = 10*n 
    NP >= 10*n || println("PSO:: sugerimos utilizar NP>10n se possível. Continuando...")

    # Verifica se temos niter > 1
    niter > 1 || error("PSO:: número de iterações deve ser maior do que um")

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
 
    # Array que contem o valor o objetivo  para cada partícula
    F = zeros(NP)

    # Melhor valor global das pariculas e sua posição
    global_best = maxintfloat(1.0)
    x_global_best = zeros(n)

    # Melhor posição e objetivo de cada partícula
    personal_best = zeros(NP)
    X_personal_best = zeros(NP,n)

    #
    # Inicializa F
    #
    Threads.@threads for p=1:NP

          # Objetivo para esta partícula
          F[p] = f(vec(X[p,:]))

    end #p
 
    # Inicializa global_best, x_global_best, personal_best e X_personal_best
    for p=1:NP

          # Aqui o objetivo é o melhor da partícula
          personal_best[p] = F[p]
          X_personal_best[p,:] .= vec(X[p,:])

          # Se também for o melhor de todos, armazena
          if F[p] < global_best
               global_best = F[p]
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

         # Calcula o objetivo para cada partícula
         Threads.@threads  for p=1:NP

              # Calcula o novo objetivo para esta partícula
              F[p] =  f(vec(X[p,:]))

         end #p      

         # Processa as informações de melhor valor para cada partícula 
         # e também para o enxame
         for p=1:NP

               # Verifica se ela melhorou
               if F[p] <  personal_best[p] 
                    personal_best[p]  = F[p]
                    X_personal_best[p,:] .= vec(X[p,:])
               end 

               # Verifica se ela é melhor do que o gloal best
               if F[p]<global_best
                    global_best = F[p]
                    x_global_best .= vec(X[p,:])
               end 

         end #p

         # Mostra a evolução (verbose=true)
         verbose && println("Melhor objetivo até agora $(global_best)")

    end #iter

   #
   # Retorna a melhor solução
   #
   return  x_global_best

end
