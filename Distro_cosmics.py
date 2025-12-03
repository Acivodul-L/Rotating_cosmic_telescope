import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Dati forniti
pairs = np.array([3464, 6025, 10708, 3272, 11643, 5266, 3278])
triples = np.array([1513, 4519, 8404, 2785, 9321, 3819, 1442])
times = np.array([1680, 855, 870, 144, 1020, 840, 1620]) * 60 # s
angles = np.array([-np.pi/2, -np.pi/3, -np.pi/4, 0, np.pi/4, np.pi/3, np.pi/2])
sigma = 0.043 # str, sottostimato
S = 0.2*0.4 # m^2

print(sigma*S)
rate_pairs = pairs/times
rate_triples = triples/times


flux_pairs = rate_pairs/(sigma*S)
flux_triples = rate_triples/(sigma*S)

# Assumendo che il dato 'sigma' sia la sezione d'urto [cm^2] e S sia l'area del detector [cm^2]
# Se sigma è in str (steradianti), e S in m^2, le unità di flusso sono (eventi/s)/(str*m^2)
d_flux_pairs = np.sqrt(pairs)/(times*sigma*S)
d_flux_triples = np.sqrt(triples*(1-triples/pairs))/(times*sigma*S)


print(f"Flusso Pairs (eventi/s/str/m^2):\n{flux_pairs}")
print("---")

# --- 1. Definizione della funzione modello ---
def fit_function(x, A, b, c, d):
    """Funzione da fittare: y = A * cos(x + b)^c + d"""
    # Uso np.power con np.abs(np.cos) per evitare problemi con esponenti non interi (c)
    # quando np.cos(x+b) è negativo.
    return A * np.power(np.abs(np.cos(x + b)), c) + d


# Parametri iniziali (p0): [A, b, c, d]
p0_guess = [200, 0.0, 2.5, 0.0] # Aggiustato 'c' a 2.0 per una stima più tipica

# --- Esecuzione dei Fit e Stampa Risultati (Codice invariato) ---

# Inizializza i parametri di fit a None per poterli usare nel plot solo se il fit riesce
popt_pairs, perr_pairs, popt_triples, perr_triples = None, None, None, None

### Fit con flux_pairs (in funzione di angles) ###
print("### Fit con Flux Pairs (X=angles, Y=flux_pairs) ###")

try:
    popt_pairs, pcov_pairs = curve_fit(
        fit_function, 
        angles, 
        flux_pairs, 
        p0=p0_guess,
        sigma=d_flux_pairs, 
        absolute_sigma=False # 'sigma' è interpretato come errore relativo
    )
    perr_pairs = np.sqrt(np.diag(pcov_pairs))

    print("Parametri ottimali [A, b, c, d]:")
    print(f"A = {popt_pairs[0]:.4f} +/- {perr_pairs[0]:.4f}")
    print(f"b = {popt_pairs[1]:.4f} +/- {perr_pairs[1]:.4f}")
    print(f"c = {popt_pairs[2]:.4f} +/- {perr_pairs[2]:.4f}")
    print(f"d = {popt_pairs[3]:.4f} +/- {perr_pairs[3]:.4f}")

except RuntimeError as e:
    print(f"Errore: Il fit con flux_pairs non è riuscito. Errore: {e}")


print("\n" + "---" + "\n")

### Fit con flux_triples (in funzione di angles) ###
print("### Fit con Flux Triples (X=angles, Y=flux_triples) ###")

try:
    popt_triples, pcov_triples = curve_fit(
        fit_function, 
        angles, 
        flux_triples, 
        p0=p0_guess,
        sigma=d_flux_triples, 
        absolute_sigma=False
    )
    perr_triples = np.sqrt(np.diag(pcov_triples))

    print("Parametri ottimali [A, b, c, d]:")
    print(f"A = {popt_triples[0]:.4f} +/- {perr_triples[0]:.4f}")
    print(f"b = {popt_triples[1]:.4f} +/- {perr_triples[1]:.4f}")
    print(f"c = {popt_triples[2]:.4f} +/- {perr_triples[2]:.4f}")
    print(f"d = {popt_triples[3]:.4f} +/- {perr_triples[3]:.4f}")

except RuntimeError as e:
    print(f"Errore: Il fit con flux_triples non è riuscito. Errore: {e}")


print("\n" + "---" + "\n")

# --- 5. Generazione del Plot ---
## Visualizzazione Grafica dei Risultati

# 1. Crea i punti per disegnare la curva fittata in modo liscio
x_fit = np.linspace(angles.min(), angles.max(), 100)

plt.figure(figsize=(10, 6))

# Aggiunge i dati sperimentali con barre d'errore (Error Bars)
plt.errorbar(angles, flux_pairs, yerr=d_flux_pairs, fmt='b+', capsize=5, label='Dati Pairs')
plt.errorbar(angles, flux_triples, yerr=d_flux_triples, fmt='r+', capsize=5, label='Dati Triples')

# Aggiunge le curve fittate
if popt_pairs is not None:
    y_fit_pairs = fit_function(x_fit, *popt_pairs)
    plt.plot(x_fit, y_fit_pairs, 'b--', label=f'Fit Pairs: $y=A\\cos(x+b)^c+d$')

if popt_triples is not None:
    y_fit_triples = fit_function(x_fit, *popt_triples)
    plt.plot(x_fit, y_fit_triples, 'r--', label=f'Fit Triples: $y=A\\cos(x+b)^c+d$')

# Configurazione del Grafico
plt.title('Fit del Flusso in funzione dell\'Angolo')
plt.xlabel('Angolo [rad]')
plt.ylabel('Flusso [eventi/s/str/m$^2$]')
plt.grid(True, linestyle='--')
plt.legend()
plt.tight_layout() # Adatta automaticamente i margini
plt.show()
