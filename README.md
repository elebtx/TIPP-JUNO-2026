# JUNO Neutrino Oscillation Simulation

Ce projet propose une chaîne d'analyse modulaire et statistique pour simuler le flux d'antineutrinos électroniques ($\bar{\nu}_e$) provenant de réacteurs nucléaires et leur détection par l'expérience **JUNO**.
L'objectif est de quantifier la sensibilité de JUNO à la **Hiérarchie de Masse des Neutrinos (NMO)** en comparant les spectres d'énergie pour l'ordonnancement Normal (NO) et Inversé (IO), tout en intégrant les effets de résolution expérimentale et les incertitudes systématiques.


## Fonctionnalités du Code
La simulation est capable de modéliser :
1. **Source** : Flux d'antineutrinos issus de la fission ($^{235}U, ^{238}U, ^{239}Pu, ^{241}Pu$) avec les paramétrisations de Vogel-Engel.
2. **Oscillations** : Probabilité de survie $P_{ee}$ complète à 3 saveurs, incluant les phases d'interférence sensibles à la hiérarchie.
3. **Détection** : Section efficace IBD, efficacité de sélection (plateau à 73%) et résolution en énergie ($3\%/\sqrt{E}$).
4. **Statistiques & Fitting** :
   -  Génération de données par Monte Carlo (Accept-Reject).
   -  Analyse de sensibilité via un test de $\chi^2$ (Asimov dataset).
   -  Prise en compte des incertitudes de forme (shape) et de normalisation (pull terms).


## Structure du Projet
Le projet est divisé en modules Python pour une flexibilité maximale :

**Modules de Physique de base**
* **`JunoPhysics.py`** : Constantes fondamentales ($N_p$, $E_{fiss}$), paramètres du détecteur et fonctions de résolution.
* **`NeutrinoFlux.py`** : Modélisation des spectres de fission et calcul du flux incident par réacteur.
* **`Oscillation.py`** : Calcul des probabilités de survie $P_{ee}$ pour les configurations NO et IO.
* **`main.py`** : Point d'entrée pour visualiser les spectres théoriques vs spectres avec résolution.
  
**Modules d'Analyse Avancée**
* **`ResolutionMonteCarlo.py`** : Simule une acquisition de 6 ans via une méthode Monte Carlo pour observer les fluctuations statistiques sur le spectre reconstruit.
* **`Chi.py`** : Effectue un scan de $\Delta m^2_{ee}$ et calcule le $\Delta\chi^2$ pour évaluer la capacité de discrimination de la hiérarchie de masse. Compare également le cas "Idéal" (1 réacteur) au cas "Réel" (10 cœurs distants).


## Résultats & Analyse
La simulation permet de générer deux types de diagnostics critiques : 
1. **Spectres d'Énergie** : Visualisation de l'écrasement des oscillations dû à la résolution énergétique et à la distribution spatiale des réacteurs (effet de baseline smearing).
2. **Sensibilité (Scan $\chi^2$)** : Graphique du $\Delta\chi^2$ en fonction de la valeur de $|\Delta m^2_{ee}|$, permettant d'identifier le minimum de la "bonne" hiérarchie et le rejet de la "mauvaise".


## Installation et Utilisation

**Prérequis**

* Python 3.x
* NumPy, SciPy, Matplotlib

**Exécution**

Pour voir les spectres comparatifs :

python main.py

Pour lancer l'analyse de sensibilité ($\chi^2$) :

python Chi.py

## Références

[1] F. An et al. [JUNO Collaboration], J. Phys. G 43 (2016) no.3, 030401 [arXiv:1507.05613].

[2] Vogel, P., & Beacom, J. F., [Angular distribution of neutrinos from inverse beta decay].
