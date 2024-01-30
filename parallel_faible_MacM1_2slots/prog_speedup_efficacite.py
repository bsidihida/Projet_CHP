import matplotlib.pyplot as plt

data = {1: 16, 2: 10, 3: 14, 4: 13}

procs = list(data.keys())
times = list(data.values())
T1 = times[0]  # Temps d'exécution avec 1 processeur
speedup = [T1 / time for time in times]
efficiency = [s / p for s, p in zip(speedup, procs)]

plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(procs, speedup, marker='o', linestyle='-', color='b')
plt.title('Speedup Curve')
plt.xlabel('Number of Processors')
plt.ylabel('Speedup')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(procs, efficiency, marker='o', linestyle='-', color='r')
plt.title('Efficiency Curve')
plt.xlabel('Number of Processors')
plt.ylabel('Efficiency')
plt.grid(True)

plt.tight_layout()
plt.show()

