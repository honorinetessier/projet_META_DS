using DelimitedFiles
using TOML
using Plots 
using Statistics

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


#normalise les vecteurs de la matrice X selon la norme min-max
function normaliser_donnees(X::Matrix{Float64})
    for j in 1:size(X, 2)
        mn = minimum(X[:, j]); mx = maximum(X[:, j])
        X[:, j] = (X[:, j] .- mn) ./ max(mx - mn, 1e-8)
    end
    return X 
end 

#stocke ensuite la taille des cluster, ce qui permet de ne pas la recalculer à chaque fois et juste de la mettre à jours 
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


#calcul l'inertie d'une solution  
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


#générer un voisin, récupérer aléaoirement un point et son nouveau cluser 
function voisin(z::Vector{Int},k::Int, n::Int)
    x_random = rand(1:n)
    cluster_random = rand(1:k)
    while cluster_random == z[x_random] 
        cluster_random = rand(1:k)
    end
    return x_random, cluster_random
end 


#règne d'acceptation de Metropolis, retourne la probabilité 
function metropolis(I::Float64,I_voisin::Float64,T::Float64)
    deltaI = I_voisin - I
    if deltaI <= 0 
        return 1
    else
        return exp(-deltaI / T)
    end
end 


#calibrage de T0 si T0_auto est TRUE sinon utilise T_init 
function calibrer_T0(X::Matrix{Float64},z::Vector{Int}, k::Int, taux_cible::Float64, cpt_cluster::Vector{Int},n::Int,d::Int)
    centroides = calculer_centroides(X,z,k,d,n,cpt_cluster)
    I = calculer_inertie(X,z,centroides,n,d)
    somme = 0.0 #somme et ctp pour faire la moyenne des dégradations 
    cpt = 0

    for i in 1:500
        x_random,cluster_random = voisin(z,k,n)
        ancien_cluster=z[x_random]

        #validation du voisin
        if cpt_cluster[ancien_cluster] <= 1
            continue
        end 

        z_voisin=copy(z)
        cpt_cluster_voisin = copy(cpt_cluster)
        z_voisin[x_random] = cluster_random
        cpt_cluster_voisin[ancien_cluster] -= 1
        cpt_cluster_voisin[cluster_random] += 1 
        centroides_voisin = calculer_centroides(X,z_voisin,k,d,n,cpt_cluster_voisin)
        I_voisin = calculer_inertie(X,z_voisin,centroides_voisin,n,d)
        deltaI = I_voisin - I
        if deltaI > 0
            somme += deltaI
            cpt+=1
        end
    end
    if cpt==0 #division par zéro 
        error("Pas de dégénération sur les 500 essais")
    end
    moyenne = somme/cpt
    return (-moyenne)/log(taux_cible)
end


#fonction du recuit simulé et retournes : centroides, inertie, hist_best, hist_cur et hist_T, et z_best 
#ici je fais le choix de rajouter des paramètres par rapport au sujet pour éviter de recalculer à chaque run la baseline
function sa_kmeans(X::Matrix{Float64}, config::Dict,z_base::Vector{Int},I_base::Float64,cpt_cluster_base::Vector{Int}, centroides_base::Matrix{Float64},n::Int,d::Int)

    k=config["data"]["k"]
    n_iter = config["sa"]["n_iter"]
    T_init = config["sa"]["T_init"]
    T_min = config["sa"]["T_min"]
    T0_auto = config["sa"]["T0_auto"]
    taux_cible = config["sa"]["taux_cible"]
    trace = config["sa"]["trace"]
    L = config["sa"]["L"]
    arret_par_stagnation = config["convergence"]["enabled"]
    patience = config["convergence"]["patience"]
    log_interval = config["sa"]["log_interval"]

    if T0_auto 
        T0 = calibrer_T0(X,z_base,k,taux_cible,cpt_cluster_base,n,d)
    else
        T0 = T_init
    end

    z_best = z_base
    I_best = I_base
    z_cur = z_base
    I_cur = I_base
    centroides = centroides_base
    hist_best = Float64[]
    hist_cur = Float64[]
    hist_T = Float64[]

    #calcul de alpha
    nombre_paliers = div(n_iter,L)
    alpha = (T_min/T0)^(1/nombre_paliers)

    T = T0
    A_plus = 0
    premier_A_plus = true #pour paciente, on attends le premier A++ pour la mettre en place 
    cpt_cluster_cur = cpt_cluster_base
    
    for i in 1:n_iter #première condition d'arrêt : le nombre d'itération
        push!(hist_best,I_best) #mise à jour des historique à chaque itération 
        push!(hist_cur,I_cur)
        push!(hist_T,T)

        if T<=T_min #deuxième condition d'arrêt : la température minimal 
            cpt_cluster_best = taille_cluster(z_best,k,n)
            centroides = calculer_centroides(X,z_best,k,d,n,cpt_cluster_best)
            return I_best,z_best,centroides,hist_best,hist_cur,hist_T
        end

        if arret_par_stagnation #troisième condition d'arrêt (optionnelle) : arrêt si pendant un certain temps pas de nouvelle amélioration 
            if A_plus >= patience && !premier_A_plus
                cpt_cluster_best = taille_cluster(z_best,k,n)
                centroides = calculer_centroides(X,z_best,k,d,n,cpt_cluster_best)
                return I_best,z_best,centroides,hist_best,hist_cur,hist_T
            end
        end

        #récupération du voisin 
        x_random,cluster_random = voisin(z_cur,k,n)
        ancien_cluster=z_cur[x_random]

        if cpt_cluster_cur[ancien_cluster] > 1 
            z_voisin=copy(z_cur)
            cpt_cluster_voisin = copy(cpt_cluster_cur)
            z_voisin[x_random] = cluster_random
            cpt_cluster_voisin[ancien_cluster] -= 1
            cpt_cluster_voisin[cluster_random] += 1 
            centroides_voisin = calculer_centroides(X,z_voisin,k,d,n,cpt_cluster_voisin)
            I_voisin = calculer_inertie(X,z_voisin,centroides_voisin,n,d)
            proba = metropolis(I_cur,I_voisin,T)
            p = rand()

            #validation du voisin 
            if p<proba
                z_cur = z_voisin
                I_cur = I_voisin
                cpt_cluster_cur = cpt_cluster_voisin
            end
        end

        #mise à jour de la meilleur solution 
        if I_cur < I_best 
            I_best = I_cur
            z_best = z_cur
            premier_A_plus = false 
            A_plus = 0
        elseif !premier_A_plus #on a déjà eu le premier A++ 
            A_plus+=1
        end 

        #affichage de la trace 
        if trace || (log_interval > 0 && mod(i, log_interval) == 0)
            println("Itération $i | T = $T | I_cur = $I_cur | I_best = $I_best ")
        end

        #mise à jour de la température en fonction des paliers L 
        if mod(i,L)==0
            T=T*alpha
        end
    end

    push!(hist_best,I_best)
    push!(hist_cur,I_cur)
    push!(hist_T,T)
    cpt_cluster_best = taille_cluster(z_best,k,n)
    centroides = calculer_centroides(X,z_best,k,d,n,cpt_cluster_best)
    return I_best,z_best,centroides,hist_best,hist_cur,hist_T
    
end


#trace les historique de best, cur et baseline
function tracer_convergence(hist_best::Vector{Float64}, hist_cur::Vector{Float64},inertie_baseline::Float64, chemin::AbstractString)
    iterations = 0:(length(hist_best) - 1)

    graphique = plot(
        iterations, 
        hist_cur, 
        label="Courante",
        xlabel="Itération", 
        ylabel="Inertie"
    )

    plot!(graphique, iterations, hist_best, label="Meilleure")
    hline!(graphique, [inertie_baseline], label="Baseline", linestyle=:dash)

    savefig(graphique, chemin)
    display(graphique)
end


#trace l'historique de température
function tracer_temperature(hist_T::Vector{Float64}, chemin::AbstractString)
    iterations = 0:(length(hist_T) - 1)

    graphique = plot(
        iterations, 
        hist_T, 
        seriestype = :steppost,
        label="Température",
        xlabel="Itération", 
        ylabel="Température"
    )

    savefig(graphique, chemin)
    display(graphique) 
end


#trace les cluster si en dimension 2 
function tracer_clusters_2d(X::Matrix{Float64}, z::Vector{Int}, centroides::Matrix{Float64}, chemin::AbstractString, inertie::Float64, k::Int,titre::String)
    graphique = plot(
        xlabel = "x1",
        ylabel = "x2",
        title = titre,
    )

    for cluster in 1:k
        lignes = (z .== cluster) #ligne[i] vaut true si le point i appartient au cluster courant 
        scatter!(graphique,
            X[lignes, 1],             
            X[lignes, 2],              
            color = cluster,           
            label = "Cluster $cluster"
        )
        #ajout des centroides 
        scatter!(graphique,
            [centroides[cluster, 1]],  
            [centroides[cluster, 2]],  
            color = cluster, #même couleur que les points       
            markershape = :star5, #forme étoile 
            markersize = 12,
            markerstrokecolor = :black,
            label = ""                 
        )
    end
    savefig(graphique, chemin)  
    display(graphique)         
end



function main()
    #récupération du fichier config 
    dossier_projet = dirname(@__DIR__)
    chemin_config = joinpath(dossier_projet,"config.toml")
    config = TOML.parsefile(chemin_config) 

    data_file = config["data"]["data_file"]
    k=config["data"]["k"]
    n_runs=config["sa"]["n_runs"]
    
    #récupération de l'instance 
    chemin_csv = joinpath(dossier_projet,data_file)

    X = charger_donnees(chemin_csv)
    n,d=size(X)
    X_norm = normaliser_donnees(X)

    #création de la solution de base 
    z_base=baseline_cyclique(k,n)
    cpt_cluster_base = taille_cluster(z_base,k,n)
    centroides_base = calculer_centroides(X,z_base,k,d,n,cpt_cluster_base)
    I_base=calculer_inertie(X,z_base,centroides_base,n,d) 

    #première runs de sa 
    I_best,z_best,centroides,hist_best,hist_cur,hist_T = sa_kmeans(X_norm, config, z_base, I_base, cpt_cluster_base,centroides_base,n,d)

    #runs suivantes de sa 
    for i in 1:n_runs-1
        I_best_2,z_best_2,centroides_2,hist_best_2,hist_cur_2,hist_T_2 = sa_kmeans(X_norm, config, z_base, I_base, cpt_cluster_base,centroides_base,n,d)
        #mise à jour de la solution si la nouvelle run est meilleur (plus petite inertie)
        if I_best_2 < I_best
            I_best = I_best_2
            z_best=z_best_2
            centroides=centroides_2
            hist_best=hist_best_2
            hist_cur=hist_cur_2
            hist_T=hist_T_2
        end
    end

    #tracer les graphiques 
    dossier_sortie = joinpath(dossier_projet, "output")
    mkpath(dossier_sortie)
    tracer_convergence(hist_best, hist_cur, I_base,joinpath(dossier_sortie, "convergence.png"))
    tracer_temperature(hist_T,joinpath(dossier_sortie, "temperature.png"))

    #si en dimension 2, retour visuel de la meilleur solution 
    if d == 2
        dossier_sortie = joinpath(dossier_projet, "output")
        mkpath(dossier_sortie)
        chemin_base = joinpath(dossier_sortie, "clusters_baseline.png")
        titre_base = "Baseline cyclique (inertie=$(round(I_base, digits=3)))"
        tracer_clusters_2d(X_norm, z_base, centroides_base,chemin_base, I_base, k, titre_base)
        chemin_sa = joinpath(dossier_sortie, "clusters_sa.png")
        titre_sa = "Recuit simulé (inertie=$(round(I_best, digits=3)))"
        tracer_clusters_2d(X_norm, z_best, centroides,chemin_sa, I_best, k, titre_sa)
    end
end 


#pour expérimentation, comme main précédent mais avec config qui est modifié, le parametre pour créer le dossier de sortie des graphiques , repetition pour respecter les 5 repetitions demandé des expérimentations 
function main_expe(config::Dict,parametre::String, valeur,repetition::Int)
    dossier_projet = dirname(@__DIR__)

    data_file = config["data"]["data_file"]
    k=config["data"]["k"]
    n_runs=config["sa"]["n_runs"]


    chemin_csv = joinpath(dossier_projet,data_file)
    X = charger_donnees(chemin_csv)
    n,d=size(X)
    X_norm = normaliser_donnees(X)
    z_base=baseline_cyclique(k,n)
    cpt_cluster_base = taille_cluster(z_base,k,n)
    centroides_base = calculer_centroides(X,z_base,k,d,n,cpt_cluster_base)
    I_base=calculer_inertie(X,z_base,centroides_base,n,d) 

    I_best,z_best,centroides,hist_best,hist_cur,hist_T = sa_kmeans(X_norm, config, z_base, I_base, cpt_cluster_base,centroides_base,n,d)
    inerties=[I_best]

    for i in 1:n_runs-1
        I_best_2,z_best_2,centroides_2,hist_best_2,hist_cur_2,hist_T_2 = sa_kmeans(X_norm, config, z_base, I_base, cpt_cluster_base,centroides_base,n,d)
        push!(inerties, I_best_2)
        if I_best_2 < I_best
            I_best = I_best_2
            z_best=z_best_2
            centroides=centroides_2
            hist_best=hist_best_2
            hist_cur=hist_cur_2
            hist_T=hist_T_2
        end
    end

    #calcul des statisques de l'expérimentation
    println("$parametre = $valeur | Baseline = $I_base | Meilleure inertie = $I_best")
    moyenne = mean(inerties)
    println("Inertie moyenne = $moyenne")
    if length(inerties) >= 2
        ecart_type = std(inerties)
        println("Écart-type = $ecart_type")
    end
    nb_iter = length(hist_cur)-1
    println("Nombre d'itération = $nb_iter")
    println("")

    #un dossier différent en fonction du paramètre testé, et un fichier par paramètres 
    dossier_sortie = joinpath(dossier_projet, "output", parametre, string(valeur),"repetition_$(repetition)")
    mkpath(dossier_sortie)
    tracer_convergence(hist_best, hist_cur, I_base,joinpath(dossier_sortie, "convergence.png"))
    tracer_temperature(hist_T,joinpath(dossier_sortie, "temperature.png"))

    if d == 2
        chemin_base = joinpath(dossier_sortie, "clusters_baseline.png")
        titre_base = "Baseline cyclique (inertie=$(round(I_base, digits=3)))"
        tracer_clusters_2d(X_norm, z_base, centroides_base,chemin_base, I_base, k, titre_base)
        chemin_sa = joinpath(dossier_sortie, "clusters_sa.png")
        titre_sa = "Recuit simulé (inertie=$(round(I_best, digits=3)))"
        tracer_clusters_2d(X_norm, z_best, centroides,chemin_sa, I_best, k, titre_sa)
    end

    return I_best
end 


#les différentes fonctions d'expérimentation 
function experimentation_taux_cible()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))
    config["sa"]["T0_auto"] = true #si on test taux_cible, on doit valider T0_auto, sinon utilise T_init de config 

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "taux_cible")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")

    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            for taux in [0.50, 0.65, 0.80]
                config["sa"]["taux_cible"] = taux #mise à jour de config 
                meilleures_inerties = Float64[]
                
                for repetition in 1:5
                    println("taux_cible = $taux | Répétition $repetition")
                    I_best =  main_expe(config,"taux_cible",taux,repetition)
                    push!(meilleures_inerties,I_best)
                end 

                moyenne = mean(meilleures_inerties)
                ecart_type = std(meilleures_inerties)
                best = minimum(meilleures_inerties)

                println("----Bilan pour taux_cible = $taux----")
                println("Meilleure inertie = $best")
                println("Moyenne des meilleures inerties = $moyenne")
                println("Écart-type = $ecart_type")
                println("")

            end
        end
    end
end

function experimentation_n_runs()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "n_runs")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")

    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            for runs in [1,3,5]
                config["sa"]["n_runs"] = runs
                meilleures_inerties = Float64[]

                for repetition in 1:5 
                    println("n_runs = $runs | Répétition $repetition")
                    I_best =  main_expe(config,"n_runs",runs,repetition)
                    push!(meilleures_inerties,I_best)
                end

                moyenne = mean(meilleures_inerties)
                ecart_type = std(meilleures_inerties)
                best = minimum(meilleures_inerties)

                println("----Bilan pour n_runs = $runs----")
                println("Meilleure inertie = $best")
                println("Moyenne des meilleures inerties = $moyenne")
                println("Écart-type = $ecart_type")
                println("")

            end
        end
    end
end

function experimentation_n_iter()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "n_iter")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")

    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            for iter in [20000,50000,100000]
                config["sa"]["n_iter"] = iter
                meilleures_inerties = Float64[]

                for repetition in 1:5 
                    println("n_iter = $iter | Répétition $repetition")
                    I_best =  main_expe(config,"n_iter",iter,repetition)
                    push!(meilleures_inerties,I_best)
                end

                moyenne = mean(meilleures_inerties)
                ecart_type = std(meilleures_inerties)
                best = minimum(meilleures_inerties)

                println("----Bilan pour n_iter = $iter----")
                println("Meilleure inertie = $best")
                println("Moyenne des meilleures inerties = $moyenne")
                println("Écart-type = $ecart_type")
                println("")
            end
        end
    end
end

#pour experimentation resultat dans le fichier txt resultat qui se trouve dans le dossier output/patience
function experimentation_patience()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))
    config["convergence"]["enabled"] = true #on doit bien valider enabled sinon le critère d'arrêt par stagnation n'est pas activé 

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "patience")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")
    
    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            for pat in [1000, 2000, 5000]
                config["convergence"]["patience"] = pat
                meilleures_inerties = Float64[]

                for repetition in 1:5 
                    println("patience = $pat | Répétition $repetition")
                    I_best =  main_expe(config,"patience",pat,repetition)
                    push!(meilleures_inerties,I_best)
                end

                moyenne = mean(meilleures_inerties)
                ecart_type = std(meilleures_inerties)
                best = minimum(meilleures_inerties)

                println("----Bilan pour patience = $pat----")
                println("Meilleure inertie = $best")
                println("Moyenne des meilleures inerties = $moyenne")
                println("Écart-type = $ecart_type")
                println("")
                
            end
        end
    end
end

function experimentation_L()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "L")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")
    
    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            for l in [1, 10, 100, 1000]
                config["sa"]["L"] = l
                meilleures_inerties = Float64[]

                for repetition in 1:5 
                    println("L = $l | Répétition $repetition")
                    I_best =  main_expe(config,"L",l,repetition)
                    push!(meilleures_inerties,I_best)
                end

                moyenne = mean(meilleures_inerties)
                ecart_type = std(meilleures_inerties)
                best = minimum(meilleures_inerties)

                println("----Bilan pour L = $l----")
                println("Meilleure inertie = $best")
                println("Moyenne des meilleures inerties = $moyenne")
                println("Écart-type = $ecart_type")
                println("")
                
            end
        end
    end
end

function experimentation_t_init()
    dossier_projet = dirname(@__DIR__)
    config = TOML.parsefile(joinpath(dossier_projet, "config.toml"))
    config["sa"]["T0_auto"] = false

    #affichage des sorties sur terminal dans un fichier txt 
    dossier_sortie = joinpath(dossier_projet, "output", "t_init")
    mkpath(dossier_sortie)
    chemin_resultats = joinpath(dossier_sortie, "resultats.txt")
    
    open(chemin_resultats, "w") do fichier
        redirect_stdout(fichier) do
            meilleures_inerties = Float64[]

            for repetition in 1:5 
                println("t_init | Répétition $repetition")
                I_best =  main_expe(config,"t_init",0.08,repetition)
                push!(meilleures_inerties,I_best)
            end

            moyenne = mean(meilleures_inerties)
            ecart_type = std(meilleures_inerties)
            best = minimum(meilleures_inerties)

            println("----Bilan pour t_inir = 0.08----")
            println("Meilleure inertie = $best")
            println("Moyenne des meilleures inerties = $moyenne")
            println("Écart-type = $ecart_type")
            println("")
            
        end
    end
end



