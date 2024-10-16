#
# Converte um vetor x para um número natural
#
# Eq 1 do texto
function Decodifica(x::Vector)

    # Número de bits do Vetor
    Nb = length(x)

    # Acumla o valor do número 
    acumulado = 0

    for j=Nb:-1:1
        bit = x[j]
        acumulado = acumulado + bit*2^(Nb-j)
    end

    return acumulado
end

#
# Converte um vetor x para um número real
# na faixa [l,u]
#
# Eq 4 do texto
function Decodifica(x::Vector,l::Float64,u::Float64)

   # Teste de consistência
   @assert u>l "Decodifica:: u deve ser maior do que l"

   # Converte x para um número natural
   nat = Decodifica(x)

   # Intervalo da variável
   intervalo = u-l

   # Resolução da nossa representação binária
   R = 2^length(x) - 1

   # Converte para real na faixa [a,b]
   l + (intervalo/R)*nat

end

#
# Dado um membro da população (uma linha da matriz P)
# extrai as variáveis de projeto e converte cada uma
# para um número real
#
# X -> uma linha da matriz P (indivíduo)
# n -> número de variáveis de projeto
# Nb -> número de bits por variável de projeto
# xmin -> vetor com os valores inferiores de cada x
# xmax -> vetor com os valores superiores de cada x
#
# Eq 17 do texto
function Extrai_Individuo(X::Vector,n::Int64,Nb::Int64,
                          xmin::Vector,xmax::Vector)
 
    # Aloca um vetor com valores reais
    saida = zeros(n)

    # Para cada variável de projeto, "captura" 
    # os bits e converte para um número real
    faixa = 1:Nb
    for j=1:n
 
        # Pega os bits
        x = X[faixa] 
 
        # Converte para real
        valor = Decodifica(x,xmin[j],xmax[j])

        # Armazena
        saida[j] = valor

        # Atualiza as posições inicial e final
        faixa = faixa .+ Nb

    end
 
   # Devolve o vetor com as variáveis de projeto
   # em representação real 
   return saida
end

#
# Torneio simples (Operação de Seleção)
#
function Torneio_Simples(valores::Vector,Np::Int)

    # Converte Np/2 para um número inteiro
    np2 = round(Int64,Np/2)

    # Aloca o vetor com os "pais"
    papis = zeros(Int64,np2)

    # Aloca o vetor com as "mamis"
    mamis = zeros(Int64,np2)

    # Loop pelos pares
    for j=1:np2

        # Sorteia os pais
        p1 = rand(1:Np)
        p2 = rand(1:Np)
        if valores[p1]<valores[p2]
            papis[j] = p1
        else
            papis[j] = p2
        end
        
        # Sorteia as mães
        m1 = rand(1:Np)
        m2 = rand(1:Np)
        if valores[m1]<valores[m2]
            mamis[j] = m1
        else
            mamis[j] = m2
        end

    end #j

    # Devolve papis e mamis
    return papis,mamis

end

#
# Recombinação - Crossover de um ponto
#
function CrossOver(papis::Vector{Int64},mamis::Vector{Int64},
                   P::Array{Bool},n::Int64,Nb::Int64,Np::Int64)

    # Cria a matriz Pbar
    Pbar = zeros(Bool,Np,n*Nb)

    # Converte Np/2 para um número inteiro
    np2 = round(Int64,Np/2)

    # Loop pelos pares
    linha_1 = 1
    linha_2 = 2
    for j=1:np2

        # Pai e mãe
        pai = papis[j]
        mae = mamis[j]

        # "Código genético"
        xpai = P[pai,:]
        xmae = P[mae,:]

        # Sorteia um ponto de crossover
        pco = rand(2:n*Nb-1)

        # Código genético dos filhos
        f1 = vcat(xpai[1:pco],xmae[pco+1:end])
        f2 = vcat(xmae[1:pco],xpai[pco+1:end])
        
        # Armazena na matriz Ptil
        Pbar[linha_1,:] .= f1
        Pbar[linha_2,:] .= f2
        
        # Incrementa as linhas
        linha_1 += 2
        linha_2 += 2

    end

    # Retorna a Pbar
    return Pbar

end

#
# Mutação
#
function Mutacao(P::Array{Bool},Np::Int64,n::Int64,Nb::Int64,taxa_mut::Float64)

    # Dado o número total de bits em P
    ntot = Np*(n*Nb)

    # Número de bits a modificar
    nmod = ceil(Int64,taxa_mut*ntot)

    # Visitamos a matriz P nmod vezes
    for j=1:nmod

        linha = rand(1:Np)
        coluna = rand(1:n*Nb)
        P[linha,coluna] = !(P[linha,coluna])

    end #j

    # Retorna a matriz modificada
    return P

end

#
# Calcula os objetivos para a população
#
# Eq 7
function Calcula_Objetivos(P::Array{Bool},obj::Function,n::Int64,Nb::Int64,xmin::Vector,xmax::Vector)

    # Número de linhas em Pai
    nlp = size(P,1)

    # Aloca a saída com a dimensão da população
    valores = zeros(nlp)

    # Loop por todos os indivíduos
    for j=1:nlp
    
        # Extrai as variáveis de projeto reais
        x = Extrai_Individuo(P[j,:],n,Nb,xmin,xmax)

        # Calcula o valor da função nesse ponto
        valores[j] = obj(x)

    end    

    return valores
end


#
# Modifica  valores e P, de acordo com a estratégia de replacement
#
function Replacement!(tipo_replacement,Np,P,Phat,valores,valores_novos)

    # No full replacement trocamos a população atual pela nova
    if tipo_replacement=="Full"

        #  População
        P .= Phat

        # objetivos
        valores .= valores_novos

    elseif tipo_replacement=="Elitismo"
        
        # ELISTISMO COMPLETO
        
        # Concatena objetivos atuais com os novos
        obj_conc = [valores ; valores_novos] 

        # Coloca em ordem crescente
        pos_ord = sortperm(obj_conc)

        # Só queremos os Np primeiros
        lista = pos_ord[1:Np]

        # Concatena as populações, só para faciltar o acesso
        P_conc = [P ; Phat] 

        # Pega os indivíduos mais aptos das duas populações
        P .= P_conc[lista,:]

        # Valores do objetivo desses caras
        valores .= obj_conc[lista]

    end

end