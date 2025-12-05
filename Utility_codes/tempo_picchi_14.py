import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# array
dati = [3.2, 3, 8, .2, 2.4, 4, 7, 4, 3.6, 4.8, 4, 5.2, 7.2, 5.6, 2.4, 3.2, 4.8, 2.8, 4.4, 3.6, 6.4, 4.4, 3.6, 6.4, 6.4, 4.4, 3.6, 3.6, 7.6, 4.4, 5.6, 3.6, 4.4, 9.2, 4.8]

# punti
counts, edges = np.histogram(dati, bins=10)
centers = 0.5 * (edges[:-1] + edges[1:])

# gaussiana
def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# stime iniziali
p0 = [max(counts), np.mean(dati), np.std(dati)]

# fit
pars, cov = curve_fit(gauss, centers, counts, p0=p0)
A_fit, mu_fit, sigma_fit = pars
print(f"gaussiana con A = {A_fit:.2f}, mu = {mu_fit:.2f}, sigma(-> risoluzione temporale) = {sigma_fit:.2f}")


# grafico
x_fit = np.linspace(min(dati), max(dati), 300)
plt.hist(dati, bins=12, edgecolor='black', fill=False)
plt.plot(x_fit, gauss(x_fit, *pars), color = "blue", label='Fit gaussiano')

plt.xlabel('Misure [ns]')
plt.ylabel('Frequenza')
plt.title("Distribuzione intervallo di tempo tra i picchi di PMT1 e PMT4")
plt.legend()
plt.show()