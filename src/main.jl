using Random 
using DelimitedFiles
using TOML
using Plots 


#lire le CSV et renvoyer la matrice X 
function charger_donnees(chemin::AbstractString)
    n=0
    d=0
    for ligne in eachline(chemin)
        resultat = match(r"n=(\d+), d=(\d+)", ligne)
        if resultat!==nothing
            n=parse(Int,resultat.captures[1]) #récupère les valeurs et convertit en entier
            d=parse(Int,resultat.captures[2])
        end
    end
    if n==0 || d==0 
        error("dimension pas trouvé")
    end
    X=zeros(n,d)
     
    i=0 #nombre de points lu

    for ligne in eachline(chemin)
        ligne = strip(ligne)
        if isempty(ligne) || startswith(ligne,"#")
            continue #on passe à la ligne suivante

        end
        valeurs = split(ligne,',')
        i+=1
        for j in 1:d
            X[i,j]=parse(Float64,valeurs[j])
        end
    end 
    
    return X

end


#normalise les vecteurs de la matrice X selon la stratégie vu au CM1 AED
function normaliser_donnees(X::Matrix{Float64})
    for j in 1:size(X, 2)
        mn = minimum(X[:, j]); mx = maximum(X[:, j])
        X[:, j] = (X[:, j] .- mn) ./ max(mx - mn, 1e-8)
    end
    return X 
end 

function taille_cluster(z::Vector{Int},k::Int,n::Int)
    cpt_cluster=zeros(Int,k)
    for i in 1:n
        cluster=z[i]
        cpt_cluster[cluster]+=1
    end
    return cpt_cluster
end 
        


#calcul le centroides de chaque cluster 
function calculer_centroides(X::Matrix{Float64}, z::Vector{Int}, k::Int, d::Int, n::Int,cpt_cluster::Vector{Int})
    centroides = zeros(k,d)
    for i in 1:n
        cluster = z[i]
        for j in 1:d
            centroides[cluster,j] = centroides[cluster,j] + X[i,j]
        end
    end
    for i in 1:k
        for j in 1:d
            centroides[i,j] = centroides[i,j]/cpt_cluster[i]
        end
    end
    return centroides 
end


#calcul l'inertie d'une solution via X, z et la matrice des centres 
function calculer_inertie(X::Matrix{Float64}, z::Vector{Int}, centroides::Matrix{Float64}, n::Int, d::Int)
    I = 0
    for i in 1:n
        cluster = z[i]
        for j in 1:d
            I = I + (X[i,j]-centroides[cluster,j])^2
        end
    end
    return I 
end



#construction de la solution initiale baselyne
function baseline_cyclique(k::Int,n::Int)
    z_base=zeros(Int,n)
    for i in 1:n
        z_base[i] = mod(i-1,k)+1
    end
    return z_base
end 


#calibrage de T0 si T0_auto est TRUE sinon utilise T_init 
function calibrer_T0(X::Matrix{Float64}, z::Vector{Int}, k::Int, taux_cible::Float64; rng::Random.AbstractRNG = Random.default_rng())
    #TODO 
end


#fonction du recuit simulé et retournes : centroides, inertie, hist_best, hist_cur et hist_T
function sa_kmeans(X::Matrix{Float64}, k::Int, cfg::Dict; rng::Random.AbstractRNG = Random.default_rng())
    #TODO 
end


#trace les historique de best, cur et baseline
function tracer_convergence(hist_best::Vector{Float64}, hist_cur::Vector{Float64},inertie_baseline::Float64, chemin::AbstractString)
    #TODO
end


#trace l'historique de température
function tracer_temperature(hist_T::Vector{Float64}, chemin::AbstractString)
    #TODO 
end


#trace les cluster si en dimension 2 
function tracer_clusters_2d(X::Matrix{Float64}, z::Vector{Int}, centroides::Matrix{Float64}, chemin::AbstractString, inertie::Float64, k::Int)
    graphique = plot(
        xlabel = "x1",
        ylabel = "x2",
        title = "Affectation cyclique (déterministe)  (inertie=$(round(inertie, digits=3)))",
    )

    for cluster in 1:k
        lignes = (z .== cluster)
        scatter!(graphique,
            X[lignes, 1],             
            X[lignes, 2],              
            color = cluster,           
            label = "Cluster $cluster"
        )
        scatter!(graphique,
            [centroides[cluster, 1]],  
            [centroides[cluster, 2]],  
            color = cluster,           
            markershape = :star5,
            markersize = 12,
            markerstrokecolor = :black,
            label = ""                 
        )
    end
    savefig(graphique, chemin)  
    display(graphique)         
end



function main()
    dossier_projet = dirname(@__DIR__)
    chemin_config = joinpath(dossier_projet,"config.toml")
    config = TOML.parsefile(chemin_config) #lit le fichier config 

    data_file = config["data"]["data_file"]
    k=config["data"]["k"]

    n_iter = config["sa"]["n_iter"]
    T_init = config["sa"]["T_init"]
    T_min = config["sa"]["T_min"]
    T0_auto = config["sa"]["T0_auto"]
    taux_cible = config["sa"]["taux_cible"]
    n_runs = config["sa"]["n_runs"]
    trace = config["sa"]["trace"]
    log_interval = config["sa"]["log_interval"]
    L = config["sa"]["L"]

    arret_par_stagnation = config["convergence"]["enabled"]
    patience = config["convergence"]["patience"]

    chemin_csv = joinpath(dossier_projet,data_file)
    X = charger_donnees(chemin_csv)
    n,d=size(X)
    X_norm = normaliser_donnees(X)
    z_base=baseline_cyclique(k,n)
    cpt_cluster = taille_cluster(z_base,k,n)
    println(z_base)
    centroides_base = calculer_centroides(X,z_base,k,d,n,cpt_cluster)
    println(centroides_base)
    I_base=calculer_inertie(X,z_base,centroides_base,n,d) 
    println(I_base)
    if d==2 
        dossier_sortie = joinpath(dossier_projet, "output")
        chemin_sortie = joinpath(dossier_sortie, "clusters_baseline.png")
        tracer_clusters_2d(X,z_base,centroides_base,chemin_sortie,I_base,k)
    end
end 

