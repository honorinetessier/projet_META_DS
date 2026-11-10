using Random 


#lire le CSV et renvoyer la matrice X 
function charger_donnees(chemin::AbstractString)
    #TODO 
end


#normalise les vecteurs de la matrice X selon la stratégie vu au CM1 AED
function normaliser_donnees(X::Matrix{Float64})
    #TODO 
end 


#calcul le centroides de chaque cluster 
function calculer_centres(X::Matrix{Float64}, z::Vector{Int}, k::Int)
    #TODO 
end


#calcul l'inertie d'une solution via X, z et la matrice des centres 
function calculer_inertie(X::Matrix{Float64}, z::Vector{Int}, centres::Matrix{Float64})
    #TODO 
end



#construction de la solution initiale baselyne
function baseline_cyclique(X::Matrix{Float64}, k::Int)
    #TODO 
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
function tracer_clusters_2d(X::Matrix{Float64}, z::Vector{Int}, centres::Matrix{Float64}, chemin::AbstractString)
    #TODO 
end



function main()
    #TODO 
end 

