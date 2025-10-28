#  Résolution Numérique de l’Équation d’Advection-Diffusion 2D

##  Aperçu du projet

Ce projet implémente la **résolution numérique de l’équation d’advection-diffusion 2D** en utilisant la **méthode des volumes finis (MVF)**.  
Le but est d’analyser la compétition entre **advection** et **diffusion** dans différents régimes physiques, à travers plusieurs cas d’étude.

Le code, écrit en **C++**, repose sur la discrétisation spatiale et temporelle de l’équation de la chaleur appliquée à un fluide (l’eau), et met en œuvre plusieurs schémas numériques pour évaluer leur **stabilité**, **convergence** et **précision**.

 Code source : [github.com/mousdid/Advection](https://github.com/mousdid/Advection)

---

##  Objectifs

- Étudier la stabilité et la convergence des schémas numériques :
  - **Explicite centré (instable)**
  - **Explicite upwind (condition CFL)**
  - **Implicite centré et upwind (stables)**
- Simuler le **réchauffement d’un fluide** dans une bouilloire chauffée par des résistances cylindriques.
- Étudier l’impact :
  - Du **maillage** sur la précision et le temps de calcul.
  - De la **puissance thermique** et de la **vitesse du fluide** sur la distribution de température.
- Identifier les **conditions de stabilité** (CFL, Péclet) et la **vitesse optimale** évitant l’ébullition.

---

##  Motivation scientifique

Le transfert de chaleur dans les fluides est un phénomène fondamental en mécanique des fluides et en thermique.  
Cette étude illustre comment les méthodes numériques (ici les volumes finis) permettent de :
- Approcher des solutions analytiques impossibles à obtenir autrement.  
- Étudier la stabilité numérique et la précision selon le choix du schéma.  
- Comprendre le compromis entre **temps de calcul**, **raffinement spatial**, et **qualité des résultats**.

---

##  Méthodologie

### 1. Validation du schéma
- Étude sur un **cas carré de référence**.
- Analyse de l’instabilité du schéma explicite centré.
- Vérification de la stabilité du schéma upwind selon la **condition CFL** :
 
  $$\(\text{CFL} = \frac{u \, \Delta t}{\Delta x}
  \)$$
- Convergence testée en **temps** et en **espace** (erreur $$\( L_\infty \)$$).

### 2. Réchauffement d’un fluide par résistances
- Simulation du chauffage de l’eau dans une bouilloire sans convection.
- Conditions limites :
  - **Dirichlet** : température imposée au bord externe.
  - **Neumann** : flux de chaleur sur les résistances.
- Étude de la température moyenne et maximale en fonction :
  - Du **temps**
  - De la **taille du maillage**
  - De la **puissance de chauffe**

### 3. Réchauffement avec écoulement potentiel
- Introduction d’un champ de vitesse potentiel :
  

  $$\(V_r = (1 - R^2/r^2) U \cos(\theta)\)$$ ;
  $$\(V_\theta = -(1 + R^2/r^2) U \sin(\theta)\)$$
  
  
- Étude de la compétition **advection–diffusion** selon la vitesse U :
  - Calcul des **nombres de Reynolds et de Péclet**.
  - Analyse de la température moyenne selon différentes vitesses.
  - Détermination de la **vitesse optimale (~1.77×10⁻⁴ m/s)** évitant l’ébullition.

---

##  Principaux résultats

| Étude | Observation clé |
|--------|------------------|
| Schéma explicite centré | Instable, divergence rapide |
| Schéma upwind explicite | Stable sous CFL |
| Schéma upwind implicite | Stable et précis |
| Maillage fin | Précision ↑ mais temps de calcul ↑ |
| Vitesse fluide ↑ | Homogénéisation, ébullition évitée |
| Puissance 2200 W | Température moyenne de 70 °C atteinte en ~603 s |
| Vitesse optimale | ~1.77×10⁻⁴ m/s → 23.6 °C moyenne, sans ébullition |

---

##  Structure du projet

```
Advection/
│
├── src/                 # Code C++ principal
│   ├── main.cpp
│   ├── Function.cpp
│   ├── DataFile.cpp
│   └── ...
│
├── data/                # Fichiers de paramètres (.toml)
│   ├── data_ref.toml
│   ├── data_real_case_1.toml
│   └── data_real_case_2.toml
│
├── mesh/                # Maillages GMSH (.geo, .mesh)
│   ├── square.mesh
│   ├── kettle2D.geo
│   └── circle_with_hole.geo
│
├── output/              # Résultats (Paraview, CSV, logs)
│
└── README.md
```

---

##  Exécution

### 1. Compiler le projet
```bash
mkdir build && cd build
cmake ..
make
```

### 2. Lancer une simulation
```bash
./advection data/data_real_case_1.toml
```

### 3. Visualiser les résultats
Les sorties peuvent être ouvertes avec **Paraview** :
```bash
paraview output/*.vtu
```

---

##  Outils utilisés

| Type | Outil / Bibliothèque |
|------|----------------------|
| Langage principal | **C++** |
| Visualisation | **Paraview**, **Gmsh** |
| Maillage | **GMSH** (.geo, .mesh) |
| Compilation | **CMake**, **g++** |
| Analyse et tracés | **Python (Matplotlib, NumPy)**|

---

##  Références

- Cours de **Volumes Finis et Équations aux Dérivées Partielles (MATMECA – ENSEIRB)**  
- Théorie de l’advection-diffusion :  
 $$\(\rho c_p \frac{\partial T}{\partial t} + (u \cdot \nabla)T - \lambda \Delta T = 0\)$$
- Logiciels : [Gmsh](https://gmsh.info), [Paraview](https://www.paraview.org)

---

##  Auteur

**OUSDID Mohamed Yassir**  
_Projet réalisé dans le cadre du cours de POO en C++ ( Janvier 2024)._
