import numpy as np
import matplotlib.pyplot as plt

# 1. Definizione dei vettori numpy
width = np.array([50, 47, 54, 52, 53, 51, 59, 45, 41, 33])
up    = np.array([32, 44, 33, 44, 39, 38, 36, 22, 12, 11])
down  = np.array([14, 18, 31, 24, 38, 21, 37, 14, 11, 7])

# 2. Ordinamento dei dati
sort_indices = np.argsort(width)
width_sorted = width[sort_indices]
up_sorted    = up[sort_indices]
down_sorted  = down[sort_indices]

# 3. Calcolo dell'Asimmetria (Y)
numerator = up_sorted - down_sorted
denominator = up_sorted + down_sorted
Y = numerator / denominator

# 4. Propagazione degli errori
sigma_Y = 2 * np.sqrt( (up_sorted * down_sorted) / (denominator**3) )

# 5. Plotting con Errorbar
plt.figure(figsize=(10, 6))

plt.errorbar(width_sorted, Y, yerr=sigma_Y, fmt='o--', color='blue', 
             ecolor='red', capsize=5, label='Dati con errore statistico')

plt.xlabel('Larghezza PMT4-discriminato [ns]')
plt.ylabel(r'Asimmetria $\frac{UP - DOWN}{UP + DOWN}$')
plt.title('Asimmetria vs Larghezza PMT4 (con barre d\'errore)')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

# Imposta il limite dell'asse Y
plt.ylim(-1, 1)

plt.show()
