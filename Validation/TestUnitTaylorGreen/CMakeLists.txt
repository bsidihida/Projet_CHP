cmake_minimum_required(VERSION 3.1)

project(TaylorGreen)

# Set project sources
set(SRCS
     TaylorGreen.cpp
     prediction-projection/Prediction.cpp
     prediction-projection/Projection.cpp
     UnitTestTaylor.cpp  # Ajout du fichier de test
)

# Trouver les paquets nécessaires (si nécessaire)
find_package(NEOS)

# Créer un exécutable pour le projet principal
add_executable(TaylorGreen ${SRCS})
target_link_libraries(TaylorGreen neos::neos)

# (Optionnel) Créer un exécutable séparé pour les tests si nécessaire
# add_executable(UnitTestTaylor UnitTestTaylor.cpp)
# target_link_libraries(UnitTestTaylor neos::neos)
