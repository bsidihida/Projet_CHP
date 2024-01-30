import matplotlib.pyplot as plt

# Données fournies
procs_efficacite = [1, 2, 3, 4]  # Nombre de processeurs
efficacites = [1, 0.183, 0.031, 0.011]  # Efficacité correspondante pour chaque nombre de processeurs

# Création du graphique pour l'efficacité
plt.figure(figsize=(10, 6))
plt.plot(procs_efficacite, efficacites, marker='o', color='green')

# Titre et étiquettes
plt.title("Courbe d'Efficacité en fonction du nombre de processeurs")
plt.xlabel("Nombre de processeurs")
plt.ylabel("Efficacité")

# Ajout de grille
plt.grid(True)

# Affichage du graphique
plt.show()
