import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- 1. Definizione del Modello ---
def custom_model(x, A, B, alpha, C):
    """
    Modello matematico da fittare: A * cos(x + B)^alpha + C
    
    Si utilizza np.abs() per la stabilità numerica (evitare numeri complessi)
    nel caso in cui np.cos(x + B) sia negativo e alpha non sia un intero.
    """
    return A * (np.abs(np.cos(x + B)) ** alpha) + C

# --- 2. Dati Esempio con Incertezze ---
# Vettore x
x_data = np.array([0, np.pi/12, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi])

# Vettori y (Dati)
triples = np.array([312, 287, 224, 132, 74, 16, 264])
pairs = np.array([363, 345, 270, 168, 97, 37, 307])

# Vettori delle incertezze su y (sigma_y)
pairs_unc = np.sqrt(pairs)
# Questa formula di incertezza (sqrt(N * (1 - N/M))) richiede che pairs (M) sia sempre >= triples (N)
# Se triples > pairs, il fit fallirà. Assumiamo che i dati rispettino pairs >= triples.
triples_unc = np.sqrt(triples * (1 - triples / pairs) ) 

# Parametri iniziali (stima iniziale comune per i 4 parametri: A, B, alpha, C)
p0_initial = [500.0, 0.0, 2.5, 0.0] 

# --- 3. Esecuzione dei Fit ---

# Dizionari per memorizzare i risultati
results = {
    'pairs': {'data': pairs, 'unc': pairs_unc, 'popt': None, 'pcov': None},
    'triples': {'data': triples, 'unc': triples_unc, 'popt': None, 'pcov': None}
}

for name in results:
    y_data = results[name]['data']
    y_unc = results[name]['unc']
    
    print(f"\n Esecuzione del Fit ({name.capitalize()})...")
    
    try:
        popt, pcov = curve_fit(
            f=custom_model, 
            xdata=x_data, 
            ydata=y_data, 
            p0=p0_initial, 
            sigma=y_unc, 
            absolute_sigma=False
        )
        
        # Salvataggio dei risultati
        results[name]['popt'] = popt
        results[name]['pcov'] = pcov
        
        # Calcolo degli errori standard e del chi-quadro ridotto
        perr = np.sqrt(np.diag(pcov)) 
        
        residuals = y_data - custom_model(x_data, *popt)
        chi_squared = np.sum((residuals / y_unc)**2)
        dof = len(x_data) - len(popt)
        chi_squared_reduced = chi_squared / dof if dof > 0 else np.nan

        # Stampa dei risultati
        A_fit, B_fit, alpha_fit, C_fit = popt
        A_err, B_err, alpha_err, C_err = perr
        
        print(f"--- Risultati {name.capitalize()} ---")
        print(f"Parametri trovati (A, B, alpha, C):")
        print(f"A: {A_fit:.4f} ± {A_err:.4f}")
        print(f"B: {B_fit:.4f} ± {B_err:.4f}")
        print(f"alpha: {alpha_fit:.4f} ± {alpha_err:.4f}")
        print(f"C: {C_fit:.4f} ± {C_err:.4f}")
        print(f"Chi-quadro ridotto: {chi_squared_reduced:.3f}")

    except RuntimeError:
        print(f"Fit {name.capitalize()} fallito: Impossibile trovare i parametri ottimali.")
    except ZeroDivisionError:
        print(f"Fit {name.capitalize()} fallito: Errore di divisione per zero (controlla incertezze).")


# --- 4. Visualizzazione dei Risultati ---

# Controlla se entrambi i fit sono riusciti prima di plottare
if results['pairs']['popt'] is not None and results['triples']['popt'] is not None:
    plt.figure(figsize=(12, 7))

    # Generazione dei punti per le curve di fit
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 500)

    # Plot del Fit PAIRS
    popt_pairs = results['pairs']['popt']
    y_fit_pairs = custom_model(x_fit, *popt_pairs)
    
    plt.errorbar(x_data, pairs, yerr=pairs_unc, fmt='o', color='blue', 
                 ecolor='lightblue', capsize=4, label='Dati Pairs con incertezze')
    
    plt.plot(x_fit, y_fit_pairs, 'b-', 
             label=f'Fit Pairs: A={popt_pairs[0]:.2f}, \u03B1={popt_pairs[2]:.2f}')

    # Plot del Fit TRIPLES
    popt_triples = results['triples']['popt']
    y_fit_triples = custom_model(x_fit, *popt_triples)

    plt.errorbar(x_data, triples, yerr=triples_unc, fmt='s', color='red', 
                 ecolor='salmon', capsize=4, label='Dati Triples con incertezze')
                 
    plt.plot(x_fit, y_fit_triples, 'r--', 
             label=f'Fit Triples: A={popt_triples[0]:.2f}, \u03B1={popt_triples[2]:.2f}')


    plt.title(r'Fit Non Lineare $y = A \cdot |\cos(x+B)|^\alpha + C$ (Pairs vs Triples)')
    plt.xlabel('Angolo x (radianti)')
    plt.ylabel('Conteggi y')
    plt.legend()
    plt.grid(True, linestyle=':')
    plt.show()

else:
    print("\n Plot non generato perché almeno uno dei fit è fallito.")