# Clustering K-means par recuit simulé

## Organisation du projet

- `src/` : code Julia du projet, dans le fichier `main.jl`.
- `config.toml` : configuration des données et des paramètres du recuit simulé.
- `data/` : instances à tester.
- `output/` : résultats des expérimentations et graphiques.


## Exécution

Depuis le dossier `src`, ouvrir un terminal et lancer Julia :

```bash
julia
```

Charger ensuite le code :

```julia
include("main.jl")
```

Puis appeler la fonction souhaitée, par exemple :

```julia
main()
```

Après une modification du code, recharger le fichier avec `include("main.jl")`.

## Fonctions disponibles

| Fonction | Action |
|---|---|
| `main()` | Lance le recuit avec les paramètres actuels de `config.toml`. |
| `experimentation_L()` | Teste les valeurs de `L` : `1`, `10`, `100`, `1000`. |
| `experimentation_patience()` | Teste les valeurs de `patience` : `1000`, `2000`, `5000`. |
| `experimentation_taux_cible()` | Teste les valeurs de `taux_cible` : `0.50`, `0.65`, `0.80`. |
| `experimentation_n_runs()` | Teste les valeurs de `n_runs` : `1`, `3`, `5`. |
| `experimentation_n_iter()` | Teste les valeurs de `n_iter` : `20000`, `50000`. |



## Protocole expérimental

Chaque valeur de paramètre est testée sur **cinq répétitions**.

Chaque répétition lance `n_runs` recuits indépendants et conserve la meilleure inertie. Par exemple, avec `n_runs = 3`, les cinq répétitions représentent au total **15 recuits**.

Le bilan indique notamment :

- la meilleure inertie observée ;
- la moyenne des meilleures inerties des cinq répétitions ;
- leur écart-type ;
- le nombre d’itérations réellement effectuées.


## Résultats

Le dossier `output/` regroupe les résultats par instance, puis par paramètre étudié.

Par exemple, pour `donnees_didactique`, les expériences sont réparties dans les dossiers `L`, `n_iter`, `n_runs`, `patience` et `taux_cible`.

Chaque dossier d’expérimentation contient :

- un compte rendu TXT des résultats ;
- un sous-dossier par valeur testée ;
- dans chaque sous-dossier, les graphiques des cinq répétitions.

Les graphiques produits sont :

- **Convergence** : inertie courante, meilleure inertie et baseline.
- **Température** : évolution de la température avec les paliers de refroidissement.
- **Clusters de la baseline** : affectation initiale et centroïdes, pour les instances en dimension 2.
- **Clusters après recuit** : meilleure affectation trouvée et centroïdes, pour les instances en dimension 2.

Pendant les expérimentations, les messages sont enregistrés dans les fichiers TXT à la place du terminal, pour pouvoir garder une trace. 
