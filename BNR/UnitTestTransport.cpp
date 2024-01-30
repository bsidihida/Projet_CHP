#include <gtest/gtest.h>
#include "Neos.hpp"  // Assurez-vous que ce chemin est correct pour inclure votre programme

// Test de la fonction dist
TEST(NeosTest, FunctionDistTest) {
    double x = 0.0; // Exemple de coordonnée x
    double y = 1.0; // Exemple de coordonnée y

    double expectedDistance = 0.5; // La distance attendue (y - 0.5)

    // Appel de la fonction dist
    double result = dist(x, y);

    // Vérification si le résultat est égal à la distance attendue
    EXPECT_DOUBLE_EQ(result, expectedDistance);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
