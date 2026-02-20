
```markdown
# Simulation des Oscillations de Neutrinos – Expérience JUNO

[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ce projet propose une chaîne d'analyse modulaire en Python pour simuler le flux d'antineutrinos électroniques ($\bar{\nu}_e$) provenant de réacteurs nucléaires et leur détection par l'expérience **JUNO** (Jiangmen Underground Neutrino Observatory). 

L'objectif principal est d'évaluer la sensibilité du détecteur à la **Hiérarchie de Masse des Neutrinos (NMO)** en comparant de manière statistique les spectres attendus pour l'ordonnancement normal (NO) et l'ordonnancement inversé (IO).

## 🚀 Fonctionnalités Principales

* **Flux de réacteurs** : Modélisation du flux de $\bar{\nu}_e$ issus de la fission des isotopes nucléaires majeurs (U-235, U-238, Pu-239, Pu-241) via la paramétrisation de Vogel-Engel.
* **Oscillations à 3 saveurs** : Calcul de la probabilité de survie des neutrinos intégrant les termes solaires et atmosphériques.
* **Détection IBD** : Modélisation de la détection basée sur la cinématique de la désintégration bêta inverse (Inverse Beta Decay).
* **Simulation Monte-Carlo** : Génération d'événements stochastiques sur 6 ans intégrant la résolution en énergie nominale du détecteur ($3\%/\sqrt{E}$).
* **Analyse Statistique** : Calcul du $\Delta\chi^2$ permettant de quantifier la capacité de JUNO à discriminer les modèles NO et IO, avec comparaison entre des configurations de cœurs idéales et réalistes.

## 🗂️ Structure du Projet

Le projet est divisé en modules Python spécialisés pour faciliter l'évolution vers des modèles d'analyse de données plus complexes :

* **`main.py`** : Script principal orchestrant la simulation sur les réacteurs (complexes de Yangjiang et Taishan), appliquant les oscillations et générant les spectres finaux.
* **`JunoPhysics.py`** : Regroupe les constantes physiques fondamentales, les paramètres du détecteur (masse cible, efficacité de 73%) et le calcul de la section efficace IBD.
* **`NeutrinoFlux.py`** : Gère la source en calculant le spectre de fission par isotope, le taux de fission par cœur et le flux géométrique incident.
* **`Oscillation.py`** : Implémente les calculs de la probabilité de survie $\bar{\nu}_e \rightarrow \bar{\nu}_e$ en fonction de l'énergie et de la distance propre à chaque réacteur.
* **`ResolutionMonteCarlo.py`** : Construit une approche statistique par acceptation/rejet pour tirer les événements reconstruits et appliquer un lissage simulant la réponse du détecteur sur 6 ans.
* **`Chi.py`** : Calcule le $\chi^2$ asimovien en intégrant des incertitudes de forme et de normalisation, puis génère les courbes de sensibilité $\Delta\chi^2$ en fonction de $|\Delta m^2_{ee}|$.

## ⚙️ Installation et Utilisation

### Prérequis
Assurez-vous d'avoir installé **Python 3.x**. Les bibliothèques scientifiques standards sont requises.

```bash
# Cloner le dépôt
git clone [https://github.com/votre-nom-utilisateur/votre-depot-juno.git](https://github.com/votre-nom-utilisateur/votre-depot-juno.git)
cd votre-depot-juno

# Installer les dépendances
pip install numpy scipy matplotlib

```

### Exécution

Vous pouvez lancer les différents modules de manière indépendante selon l'analyse souhaitée :

```bash
# 1. Lancer la simulation complète (spectres théoriques et effectifs)
python main.py

# 2. Tester la génération de spectres stochastiques (Monte-Carlo)
python ResolutionMonteCarlo.py

# 3. Générer les courbes de sensibilité statistiques (Analyse Chi2)
python Chi.py

```

## 📊 Résultats Attendus

L'exécution des scripts génère plusieurs analyses visuelles clés :

1. **Spectres Théoriques** : Comparaison des spectres énergétiques sans oscillation, avec NO, et avec IO.
2. **Spectres Reconstruits** : Mise en évidence de l'effet de convolution de la résolution en énergie qui brouille le signal théorique.
3. **Taux d'Événements** : Calcul du taux absolu d'événements attendus par jour sur une prise de données de 6 ans.
4. **Sensibilité Statistique** : Une courbe  illustrant mathématiquement la zone d'exclusion de la mauvaise hiérarchie selon les paramètres d'oscillation.

## 🔮 Évolutions prévues (Objectifs TIPP)

Conformément aux objectifs d'analyse, les prochaines étapes de développement incluent :

* L'intégration des bruits de fond majeurs (géoneutrinos, muons cosmiques, radioactivité locale).
* L'affinement de la fonction de résolution du détecteur (effets de non-linéarité de la réponse du scintillateur liquide).
* L'optimisation des routines de minimisation  pour contraindre davantage les paramètres de mélange (, ).

## 📚 Références

* F. An et al. [JUNO Collaboration], *Neutrino Physics with JUNO*, J. Phys. G 43 (2016) no.3, 030401 [arXiv:1507.05613].

```

***

Est-ce que tu aimerais que je te génère également un fichier `requirements.txt` contenant les versions exactes des bibliothèques (`numpy`, `scipy`, `matplotlib`) pour accompagner ton dépôt ?

```
