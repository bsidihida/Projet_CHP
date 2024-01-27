#include <iostream>
#include <cmath>
#include "Neos.hpp" // Assurez-vous d'inclure les fichiers nécessaires pour NPoint

bool test_taylorGreen_Ux() {
    NPoint pt = {M_PI / 2, M_PI / 2}; // Exemple de point
    double t = 1.0; // Exemple de temps
    double expected = sin(pt[0]) * cos(pt[1]) * exp(-2*0.01*t); // Valeur attendue
    double result = taylorGreen_Ux(pt, t);
    return std::abs(result - expected) < 1e-6; // Vérifiez la précision
}

bool test_taylorGreen_Uy() {
    NPoint pt = {M_PI / 2, M_PI / 2}; // Exemple de point
    double t = 1.0; // Exemple de temps
    double expected = -cos(pt[0]) * sin(pt[1]) * exp(-2*0.01*t); // Valeur attendue
    double result = taylorGreen_Uy(pt, t);
    return std::abs(result - expected) < 1e-6; // Vérifiez la précision
}

int main() {
    bool uxTestPassed = test_taylorGreen_Ux();
    bool uyTestPassed = test_taylorGreen_Uy();

    if (uxTestPassed) {
        std::cout << "test_taylorGreen_Ux PASSED." << std::endl;
    } else {
        std::cout << "test_taylorGreen_Ux FAILED." << std::endl;
    }

    if (uyTestPassed) {
        std::cout << "test_taylorGreen_Uy PASSED." << std::endl;
    } else {
        std::cout << "test_taylorGreen_Uy FAILED." << std::endl;
    }

    return 0;
}
