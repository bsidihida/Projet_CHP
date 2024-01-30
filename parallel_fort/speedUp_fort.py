import matplotlib.pyplot as plt

# Données fournies
procs = [1, 2, 3, 4]  # Nombre de processeurs
speedUps = [1, 0.367, 0.093, 0.044]  # SpeedUp correspondant pour chaque nombre de processeurs

# Création du graphique
plt.figure(figsize=(10, 6))
plt.plot(procs, speedUps, marker='o')

# Titre et étiquettes
plt.title("Courbe de SpeedUp en fonction du nombre de processeurs")
plt.xlabel("Nombre de processeurs")
plt.ylabel("SpeedUp")

# Ajout de grille
plt.grid(True)

# Affichage du graphique
plt.show()
