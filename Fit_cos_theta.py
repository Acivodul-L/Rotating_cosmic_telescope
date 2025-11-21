import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- 1. Definizione del Modello ---
def custom_model(x, A, B, alpha, C):
    """
    Modello matematico da fittare: A * cos(x + B)^alpha + C
    Si utilizza np.abs() per stabilità numerica.
    """
    return A * (np.abs(np.cos(x + B)) ** alpha) + C

# --- 2. Dati e Costanti ---
# Vettore x
angles = np.array([0, np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi])

# Vettori y (Conteggi)
triples = np.array([312, 287, 224, 132, 74, 16, 264])
pairs = np.array([363, 345, 270, 168, 97, 37, 307])
singles_1 = np.array([15214, 14794, 13799, 11703, 10643, 17634, 13607])
singles_4 = np.array([48655, 48802, 20104, 43267, 41415, 38915, 41242])

DT = 1000.0   # s, Tempo di integrazione
eff = 0.95    # Efficienza (senza errore associato come richiesto)
S = 0.2 * 0.4 # m^2, Area
Omega = 0.043 # str, Angolo solido

# Calcolo Accidentali
# Nota: La formula delle coincidenze casuali (Rate = R1*R2*DeltaT) dipende dalla larghezza delle finestre
# Assumo che i numeri (51e-9, etc.) siano le larghezze temporali corrette in secondi.
width_coinc = (51e-9 + 78e-9 - 2*2e-9) 
R_pairs_fake = (singles_1/DT) * (singles_4/DT) * width_coinc
R_triples_fake = 0 # Trascurabile o non definito

# --- 3. Calcolo Flussi e Incertezze ---

# Incertezze statistiche sui conteggi (Poisson: sqrt(N))
# NOTA: Usiamo sqrt(N) perché stiamo calcolando un flusso assoluto, 
# non un rapporto di efficienza binomiale.
d_pairs_counts = np.sqrt(pairs)
d_triples_counts = np.sqrt(triples) 

# Fattore di normalizzazione comune
norm_factor = Omega * S * DT * eff

# Calcolo Flussi
flux_pairs = (pairs - (DT * R_pairs_fake)) / norm_factor
flux_triples = (triples - (DT * R_triples_fake)) / norm_factor

# Propagazione errori (Semplificata senza d_eff)
# Se F = N / K, allora sigma_F = sigma_N / K
d_flux_pairs = d_pairs_counts / norm_factor
d_flux_triples = d_triples_counts / norm_factor

# Parametri iniziali (A, B, alpha, C)
p0_initial = [100, 0.0, 2.0, 0.0] 

# --- 4. Esecuzione dei Fit ---

results = {
    'pairs': {'data': flux_pairs, 'unc': d_flux_pairs, 'popt': None},
    'triples': {'data': flux_triples, 'unc': d_flux_triples, 'popt': None}
}

for name in results:
    y_data = results[name]['data']
    y_unc = results[name]['unc']
    
    print(f"\n--- Fit {name.capitalize()} ---")
    
    try:
        # absolute_sigma=True è CRUCIALE quando si forniscono errori fisici reali (y_unc)
        popt, pcov = curve_fit(
            f=custom_model, 
            xdata=angles, 
            ydata=y_data, 
            p0=p0_initial, 
            sigma=y_unc, 
            absolute_sigma=True 
        )
        
        results[name]['popt'] = popt
        perr = np.sqrt(np.diag(pcov)) 
        
        # Calcolo Chi-Quadro
        residuals = y_data - custom_model(angles, *popt)
        chi_squared = np.sum((residuals / y_unc)**2)
        dof = len(angles) - len(popt)
        chi_red = chi_squared / dof if dof > 0 else np.nan

        print(f"A:     {popt[0]:.4f} ± {perr[0]:.4f}")
        print(f"B:     {popt[1]:.4f} ± {perr[1]:.4f}")
        print(f"alpha: {popt[2]:.4f} ± {perr[2]:.4f}")
        print(f"C:     {popt[3]:.4f} ± {perr[3]:.4f}")
        print(f"Chi^2 ridotto: {chi_red:.3f}")

    except Exception as e:
        print(f"Fit fallito: {e}")

# --- 5. Plotting ---

if results['pairs']['popt'] is not None and results['triples']['popt'] is not None:
    plt.figure(figsize=(10, 6))
    
    # Griglia densa per le curve
    x_fit = np.linspace(np.min(angles), np.max(angles), 500)

    # Plot Pairs
    plt.errorbar(angles, flux_pairs, yerr=d_flux_pairs, fmt='o', color='blue', 
                 capsize=5, label='Dati Pairs')
    plt.plot(x_fit, custom_model(x_fit, *results['pairs']['popt']), 'b-', alpha=0.7,
             label=f"Fit Pairs ($\\alpha={results['pairs']['popt'][2]:.2f}$)")

    # Plot Triples
    plt.errorbar(angles, flux_triples, yerr=d_flux_triples, fmt='s', color='red', 
                 capsize=5, label='Dati Triples')
    plt.plot(x_fit, custom_model(x_fit, *results['triples']['popt']), 'r--', alpha=0.7,
             label=f"Fit Triples ($\\alpha={results['triples']['popt'][2]:.2f}$)")

    plt.title(r'Fit Flusso di Muoni: $I(\theta) = I_0 \cos^\alpha(\theta + \phi) + C$')
    plt.xlabel('Angolo Zenitale [rad]')
    plt.ylabel('Flusso [particelle / ($s \cdot m^2 \cdot sr$)]')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plt.show()