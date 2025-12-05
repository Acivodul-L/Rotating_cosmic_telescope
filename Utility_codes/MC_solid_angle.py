import numpy as np

def geometric_factor_montecarlo(L1, W1, L2, W2, Z, N=1000000):
    """
    Calcola il Fattore Geometrico G e l'Angolo Solido Medio Effettivo Omega_eff.
    """
    
    # 1. Calcolo Aree
    A1 = L1 * W1
    A2 = L2 * W2

    # 2. Generazione casuale dei punti (omessa per brevità, è la stessa)
    x1 = np.random.uniform(-L1 / 2, L1 / 2, N)
    y1 = np.random.uniform(-W1 / 2, W1 / 2, N)
    z1 = np.zeros(N)

    x2 = np.random.uniform(-L2 / 2, L2 / 2, N)
    y2 = np.random.uniform(-W2 / 2, W2 / 2, N)
    z2 = np.full(N, Z)
    
    # 3. Vettore Distanza e Modulo (r)
    dX = x2 - x1
    dY = y2 - y1
    dZ = z2 - z1
    r_sq = dX**2 + dY**2 + dZ**2
    r = np.sqrt(r_sq)

    # 4. Calcolo degli Angoli (theta)
    cos_theta1 = dZ / r 
    cos_theta2 = dZ / r 
    
    # 5. Calcolo del Fattore Geometrico (G)
    Integrand = (cos_theta1 * cos_theta2) / r_sq
    G = A1 * A2 * np.mean(Integrand)
    
    # 6. Calcolo dell'Angolo Solido Medio Effettivo (Omega_eff)
    # L'area di riferimento è A1
    Omega_eff = G / A1 

    # 7. Calcolo dell'errore (opzionale ma fondamentale)
    # La Varianza dell'integrando I = (cos(theta1) * cos(theta2)) / r^2
    var_integrand = np.var(Integrand)
    mean_integrand = np.mean(Integrand)
    
    # Errore assoluto su G (prop. all'errore sulla media)
    delta_G = G * np.sqrt(var_integrand / N) / mean_integrand
    
    # Errore assoluto su Omega_eff
    delta_Omega = delta_G / A1
    
    return G, delta_G, Omega_eff, delta_Omega

# --- Esempio d'Uso (Assumiamo Rivelatori 20x20 a 10cm di distanza) ---
L1 = 40.0; W1 = 20.0
L2 = 40.0; W2 = 20.0
Z = 136.0
N_ITERATIONS = 500

G, delta_G, Omega_eff, delta_Omega = geometric_factor_montecarlo(L1, W1, L2, W2, Z, N=N_ITERATIONS)

print(f"Rivelatore 1 Area (A1): {L1*W1:.2f} cm^2")
print(f"Eseguito Monte Carlo con N = {N_ITERATIONS:,} iterazioni.")
print("---")
print(f"Fattore Geometrico G = {G:.4f} +/- {delta_G:.4f} cm^2 * sr")
print(f"Angolo Solido Medio Effettivo Ω_eff = {Omega_eff:.4f} +/- {delta_Omega:.4f} sr")