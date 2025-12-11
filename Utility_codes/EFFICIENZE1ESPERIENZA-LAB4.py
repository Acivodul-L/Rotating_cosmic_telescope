import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


#Dati
V2 = np.array([724, 750, 775, 788, 800,813, 825, 838, 850, 875])
triple2 = np.array([1, 8, 13, 21, 40, 34 ,35,24 ,37, 32])
doppie2 = np.array([31, 33, 28, 28, 44, 36, 38,28 ,40, 36])
epsilon2 = (triple2/doppie2)
sigmaepsilon2 = np.sqrt((epsilon2*(1-epsilon2))/(doppie2))


V1 = np.array([750, 775, 788, 800, 825, 850, 875])
triple1= np.array([8, 15, 26, 28, 33, 42, 37])
doppie1 = np.array([56, 55, 63, 57, 59, 63, 61])
varepsilon1 = (triple1/doppie1)
sigmavarepsilon1= np.sqrt((varepsilon1*(1-varepsilon1))/(doppie1))   #ERRORE DATO DALLE MISURE
epsilon1 = varepsilon1/0.565
sigmavarepsilon1_1 = 0.0035 #ERRORE RELATIVO DATO DALL'ACCETTANZA
sigmaepsilon1 = (np.sqrt((sigmavarepsilon1/varepsilon1)**2+(sigmavarepsilon1_1)**2)) * epsilon1

"""
triple3 = np.array([25 ,25 ,33 ,33, 24, 31,33])
doppie3 = np.array([257,281 ,300 ,259, 259, 275, 274])
varepsilon3 = (triple3/doppie3)
sigmavarepsilon3 = np.sqrt((varepsilon3*(1-varepsilon3))/(doppie3))

epsilon3 = varepsilon3/0.1037
#sigmaepsilon3 = np.full(len(epsilon3), 6)
"""


#RICALCOLO CON PMT2, 3, 4
V4 = np.array([ 650, 700, 725, 750, 775, 788, 800, 825, 850, 875])
triple4 = np.array([ 3, 39, 29, 52 ,51, 59, 53, 58, 61 ,57])
doppie4 = np.array([ 390, 432, 390, 370 ,369, 395, 368, 393, 374 ,373])
varepsilon4 = (triple4/doppie4)
sigmavarepsilon4 = np.sqrt((varepsilon4*(1-varepsilon4))/(doppie4))  #ERRORE DATO DALLE MISURE DI VAREPSILON
epsilon4 = varepsilon4/0.13
sigmavarepsilon4_1 = 0.0048     #ERRORE DATO DALL'ACCETTANZA (RELATIVO)
sigmaepsilon4 = (np.sqrt((sigmavarepsilon4/varepsilon4)**2+(sigmavarepsilon4_1)**2)) * epsilon4   #SOMMA IN QUADRATURA DEGLI ERR RELATIVI PER OTTENERE L'ERR RELATIVO TOT, MOLTIPLICATO POI PER LA MISURA STESSA


# funzione di fit (arctan)
def arctan(x, a, b, c, d):
    return (a/np.pi)*np.arctan(c*(x-b)) + d
    
def cost(x, k):
    return k
    
# Figura
plt.figure("Efficienza", figsize=(9, 8))
plt.grid(which="both", ls="dashed", color="gray")

x_fit = np.linspace(500, 900, 200)

# FIT PMT1
p0_1 = [0.7, 800, 0.1, 0.5]
popt1, pcov1 = curve_fit(
            f=arctan, 
            xdata=V1, 
            ydata=epsilon1, 
            p0=p0_1, 
            sigma=sigmavarepsilon1, 
            absolute_sigma=True 
        )
perr1 = np.sqrt(np.diag(pcov1))

print("==========================")
print("Stampa dei parametri PMT 1")
for i in range(len(popt1)):
    if (i== 0) :
        print(f'eff_1 = {popt1[i]:.3f} +/- {perr1[i]:.3f}') 
    else:
        print(f'{i}-param = {popt1[i]:.3f} +/- {perr1[i]:.3f}')

plt.plot(x_fit, arctan(x_fit, *popt1), color = "orange", label='Efficiency curve 1')
    
# FIT PMT2
p0_2 = [0.7, 800, 0.1, 0.5]
popt2, pcov2 = curve_fit(
            f=arctan, 
            xdata=V2, 
            ydata=epsilon2, 
            p0=p0_2, 
            sigma=sigmaepsilon2, 
            absolute_sigma=True 
        )
perr2 = np.sqrt(np.diag(pcov2))
        
plt.plot(x_fit, arctan(x_fit, *popt2), color = "green", label='Efficiency curve 2')

print("==========================")
print("Stampa dei parametri PMT 2")
for i in range(len(popt2)):
    if (i== 0) :
        print(f'eff_2 = {popt2[i]:.3f} +/- {perr2[i]:.3f}') 
    else:
        print(f'{i}-param = {popt2[i]:.3f} +/- {perr2[i]:.3f}')

# FIT PMT4
p0_4 = [1.7, 700, 0.13, 0.0]
popt4, pcov4 = curve_fit( #Fit con tutti i punti eccetto il primo (outlier) 
            f=arctan, 
            xdata=V4[1:], 
            ydata=epsilon4[1:], 
            p0=p0_4, 
            sigma=sigmavarepsilon4[1:], 
            absolute_sigma=True,
            maxfev = 5000
        )
perr4 = np.sqrt(np.diag(pcov4))

print("==========================")
print("Stampa dei parametri PMT 4")
for i in range(len(popt4)):
    if (i== 0) :
        print(f'eff_4 = {popt4[i]:.3f} +/- {perr4[i]:.3f}') 
    else:
        print(f'{i}-param = {popt4[i]:.3f} +/- {perr4[i]:.3f}')

plt.plot(x_fit, arctan(x_fit, *popt4), color = "purple", label='Efficiency curve 4')

plt.errorbar(V1, epsilon1, yerr=sigmaepsilon1,  fmt='o', color='orange' , label=r'$\epsilon_1$', capsize=3)

plt.errorbar(V2, epsilon2, yerr= sigmaepsilon2,  fmt='o',color= 'green', label=r'$\epsilon_2$', capsize=3)


plt.errorbar(V4, epsilon4, yerr=sigmaepsilon4, fmt='o', color='purple' , label=r'$\epsilon_4$', capsize=3)

# Etichette e legenda
plt.xlabel("Tensione V [V]")
plt.ylabel("Efficienza ")
plt.ylim(0, 1.3)
plt.title("Efficienza PMT1, PMT2, PMT4")
plt.legend()
plt.show()
