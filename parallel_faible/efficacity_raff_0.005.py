import matplotlib.pyplot as plt

# Données : nombre de processeurs (np) et leur efficacité correspondante
np = [1, 2, 3, 4, 5, 6, 7, 8]
efficacite = [1.000, 0.73, 0.575, 0.467, 0.33, 0.217, 0.195, 0.16]

# Création de la figure et du graphique
plt.figure(figsize=(10, 6))
plt.plot(np, efficacite, marker='o', color='r')  # couleur rouge pour la distinction

# Ajout des titres et des étiquettes
plt.title('Courbe de l\'Efficacité en fonction du nombre de processeurs (Raffinement 0,005)')
plt.xlabel('Nombre de processeurs (np)')
plt.ylabel('Efficacité')

# Ajout de la grille pour une meilleure lisibilité
plt.grid(True)

# Affichage du graphique
plt.show()
