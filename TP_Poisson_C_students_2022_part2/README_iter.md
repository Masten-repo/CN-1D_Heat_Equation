
# TP Poisson 1D - Partie 4 : Méthodes Itératives

Ce document résume les commandes nécessaires pour compiler, exécuter les benchmarks (Richardson, Jacobi, Gauss-Seidel) 

## 1. Compilation

Assurez-vous d'être dans le conteneur Docker et à la racine du projet.

```bash
# Crée le dossier bin s'il n'existe pas et compile tout
mkdir -p bin
make

```

## 2. Exécution des Méthodes et Collecte des Données

Le programme `tpPoisson1D_iter` génère un fichier `RESVEC.dat` à chaque exécution. Il faut donc renommer ce fichier entre chaque test pour ne pas l'écraser.

### Étape 2.1 : Richardson Alpha (Méthode 0)

```bash
./bin/tpPoisson1D_iter 0
mv RESVEC.dat resvec_richardson.dat

```

### Étape 2.2 : Jacobi (Méthode 1)

```bash
./bin/tpPoisson1D_iter 1
mv RESVEC.dat resvec_jacobi.dat

```

### Étape 2.3 : Gauss-Seidel (Méthode 2)

```bash
./bin/tpPoisson1D_iter 2
mv RESVEC.dat resvec_gs.dat

```

*Vérification :* Vous devez avoir 3 fichiers `.dat` dans votre dossier :

```bash
ls -l resvec_*.dat

```

---

Les courbes Rouge et Bleue se superposent parfaitement, c'est normal pour la matrice du Laplacien 1D (D = 2I et Alpha_opt ≈ 0.5).*

```

```