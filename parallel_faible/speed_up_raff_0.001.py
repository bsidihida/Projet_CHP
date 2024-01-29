import matplotlib.pyplot as plt

# Données : nombre de processeurs (np) et leur speed-up correspondant
np = [1, 2, 3, 4, 5, 6, 7, 8]
speed_up = [1.000, 1.536, 1.637, 1.817, 1.585, 1.475, 1.366, 1.318]

# Création de la figure et du graphique
plt.figure(figsize=(10, 6))
plt.plot(np, speed_up, marker='o', color='b')

# Ajout des titres et des étiquettes
plt.title('Courbe du Speed-Up en fonction du nombre de processeurs (Raffinement 0,001)')
plt.xlabel('Nombre de processeurs (np)')
plt.ylabel('Speed-Up')

# Ajout de la grille pour une meilleure lisibilité
plt.grid(True)

# Affichage du graphique
plt.show()
