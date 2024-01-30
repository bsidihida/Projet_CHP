#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// Assurez-vous que le chemin vers poissonCos.cpp est correct
#include "../Validation Poisson/poissonCos.cpp"

// Structure pour stocker les résultats de la simulation
struct PoissonCosResult {
    double errInf;
    double errL1;
    // Vous pouvez également inclure d'autres données si nécessaire
};

// Supposons que cette fonction est définie dans poissonCos.cpp
PoissonCosResult runPoissonCosSimulation();

// Fonction pour comparer deux résultats de simulations
bool compareResults(const PoissonCosResult& resultA, const PoissonCosResult& resultB) {
    const double TOLERANCE = 1e-6;
    bool isErrInfClose = std::abs(resultA.errInf - resultB.errInf) < TOLERANCE;
    bool isErrL1Close = std::abs(resultA.errL1 - resultB.errL1) < TOLERANCE;
    return isErrInfClose && isErrL1Close;
}

int main() {
    const int NUM_TESTS = 10; // Nombre d'exécutions de simulation pour tester la reproductibilité
    std::vector<PoissonCosResult> results;

    // Exécuter les simulations et stocker les résultats
    for (int i = 0; i < NUM_TESTS; ++i) {
        results.push_back(runPoissonCosSimulation());
    }

    // Comparer les résultats de toutes les simulations
    bool isReproducible = true;
    for (int i = 1; i < NUM_TESTS; ++i) {
        if (!compareResults(results[0], results[i])) {
            isReproducible = false;
            break;
        }
    }

    // Afficher le résultat du test de reproductibilité
    if (isReproducible) {
        std::cout << "Test de reproductibilité réussi." << std::endl;
    } else {
        std::cout << "Test de reproductibilité échoué." << std::endl;
    }

    return 0;
}
