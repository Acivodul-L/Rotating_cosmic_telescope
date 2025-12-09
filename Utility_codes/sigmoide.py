import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy.optimize import curve_fit

# array
t = np.array([33, 35, 37, 39, 41, 43, 45, 47, 49, 51]) # durata del segnale discriminato impostata
r_up = np.array([-6.5, -4, -2, 0, 2, 3, 4, 6, 8, 10])
r_up = r_up + 2.5 # aggiungiamo il ritardo intrinseco per avere il ritardo totale introdotto
r_down = np.array([-16.5, -12.5, -10.5, -8.5, -6.5, -4.5, -2, 0, 2])
r_down = r_down + 2.5 # aggiungiamo il ritardo intrinseco per avere il ritardo totale introdotto
'''
-4.5: senza scatola, 4 ns di cavo 1 + 2 ns di cavo 4
-10.5: senza scatola, 2 ns di cavo 1 + 6 ns di cavo 4
-12.5: con scatola per il 4, 2 ns di cavo 1 + (2 + 2) ns di cavo 4 + 2.5 di scatola + 1.5 di ritardo impostato per il 4
-16.5: con scatola per il 4, come prima ma altri 4 ns impostati alla scatola
'''

c = np.array([161, 132, 140, 158, 193, 292, 338, 397, 370, 399]) # doppie misurate
c_up = np.array([340, 376, 341, 294, 213, 145, 40, 4, 12,5]) # doppie misurate, variando il ritardo
c_down = np.array([296, 316, 328, 256, 129, 54, 10, 9, 5]) # doppie misurate, variando il ritardo


delta_c = np.sqrt(c) # errore statistico
delta_c_up = np.sqrt(c_up) # errore statistico
delta_c_down = np.sqrt(c_down) # errore statistico

# funzione di fit (sigmoide)
def sigmoid(x, a, b, c):
    return a / (1 + np.exp(-(b-x))) + c




# stime iniziali
p0_up = [344, 5, 0]
p0_down = [330, -5, 0]

# fit
pars_up, cov_up = curve_fit(sigmoid, r_up, c_up, p0=p0_up, sigma = delta_c_up)
a_fit_up, b_fit_up, c_fit_up = pars_up

pars_down, cov_down = curve_fit(sigmoid, r_down, c_down, p0=p0_down, sigma = delta_c_down)
a_fit_down, b_fit_down, c_fit_down = pars_down
print(f"sigmoide UP con a = {a_fit_up:.4f}, b = {b_fit_up:.4f}, c = {c_fit_up:.4f}")
print(f"sigmoide DOWN con a = {a_fit_down:.4f}, b = {b_fit_down:.4f}, c = {c_fit_down:.4f}")


# grafico
x_fit = np.linspace(min(r_down), max(r_up), 300)
plt.errorbar(r_up, c_up, delta_c_up, fmt = 'o', color = 'blue')
plt.errorbar(r_down, c_down, delta_c_down, fmt = 'o', color = 'red')

plt.plot(x_fit, sigmoid(x_fit, *pars_up), color = "blue", label='fit UP')
plt.plot(x_fit, sigmoid(x_fit, *pars_down), color = "red", label='fit DOWN')


plt.xlabel('Ritardo  [ns]')
plt.ylabel('Conteggi (doppie)')
plt.title("Conteggi (doppie) in funzione del ritardo introdotto")
plt.legend()
plt.show()

## asimmetria boh
plt.plot(x_fit, (sigmoid(x_fit, *pars_up)-sigmoid(x_fit, *pars_down))/(sigmoid(x_fit, *pars_up)+sigmoid(x_fit, *pars_down)), color = "orange", label='asimmetria')
plt.show()



