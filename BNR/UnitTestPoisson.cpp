#include <gtest/gtest.h>
#include "Neos.hpp" // Assurez-vous que ce chemin est correct pour inclure votre programme

// Test de la fonction u
TEST(NeosTest, FunctionUTest) {
    // Création d'un point test
    NPoint testPoint(2);
    testPoint[0] = 1.0; // x = 1.0
    testPoint[1] = 2.0; // y = 2.0

    // Appel de la fonction u avec le point test
    double result = u(testPoint);

    // Vérification si le résultat est correct
    // u = y^2 + 2*y, donc pour y = 2.0, u = 2*2 + 2*2 = 8
    EXPECT_DOUBLE_EQ(result, 8.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
