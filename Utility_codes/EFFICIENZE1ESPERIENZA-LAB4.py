import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


#Dati
V2 = np.array([724, 750, 775, 788, 800,813, 825, 838, 850, 875])

triple2 = np.array([1, 8, 13, 21, 40, 34 ,35,24 ,37, 32])
doppie2 = np.array([31, 33, 28, 28, 44, 36, 38,28 ,40, 36])
epsilon2 = (triple2/doppie2)
sigmaepsilon2 = np.sqrt((epsilon2*(1-epsilon2))/(doppie2))


V = np.array([750, 775, 788, 800, 825, 850, 875])
triple1= np.array([8, 15, 26, 28, 33, 42, 37])
doppie1 = np.array([56, 55, 63, 57, 59, 63, 61])
varepsilon1 = (triple1/doppie1)
sigmavarepsilon1= np.sqrt((varepsilon1*(1-varepsilon1))/(doppie1))   #ERRORE DATO DALLE MISURE
epsilon1 = varepsilon1/0.565
sigmavarepsilon1_1 = 0.0035 #ERRORE RELATIVO DATO DALL'ACCETTANZA

sigmaepsilon1 = (np.sqrt((sigmavarepsilon1/varepsilon1)**2+(sigmavarepsilon1_1)**2)) * epsilon1


triple3 = np.array([25 ,25 ,33 ,33, 24, 31,33])
doppie3 = np.array([257,281 ,300 ,259, 259, 275, 274])
varepsilon3 = (triple3/doppie3)
sigmavarepsilon3 = np.sqrt((varepsilon3*(1-varepsilon3))/(doppie3))

epsilon3 = varepsilon3/0.1037
#sigmaepsilon3 = np.full(len(epsilon3), 6)



#RICALCOLO CON PMT2, 3, 4
V4 = np.array([725, 750, 775, 788, 800, 825, 850, 875])
triple4 = np.array([29, 52 ,51, 59, 53, 58, 61 ,57])
doppie4 = np.array([390, 370 ,369, 395, 368, 393, 374 ,373])
varepsilon4 = (triple4/doppie4)
sigmavarepsilon4 = np.sqrt((varepsilon4*(1-varepsilon4))/(doppie4))  #ERRORE DATO DALLE MISURE DI VAREPSILON
epsilon4 = varepsilon4/0.13
sigmavarepsilon4_1 = 0.0048     #ERRORE DATO DALL'ACCETTANZA (RELATIVO)


sigmaepsilon4 = (np.sqrt((sigmavarepsilon4/varepsilon4)**2+(sigmavarepsilon4_1)**2)) * epsilon4   #SOMMA IN QUADRATURA DEGLI ERR RELATIVI PER OTTENERE L'ERR RELATIVO TOT, MOLTIPLICATO POI PER LA MISURA STESSA


#0.467
#Converti errori percentuali in errori assoluti
'''
sigmaepsilon1 = sigmaepsilon1/100 * epsilon1
sigmavarepsilon1 = sigmavarepsilon1/100 * varepsilon1
'''
#sigmaepsilon2 = sigmaepsilon2/100 * epsilon2
'''
sigmaepsilon3 = sigmaepsilon3/100 * epsilon3
sigmavarepsilon3 = sigmavarepsilon3/100 * varepsilon3
'''

# Figura
plt.figure("Efficienza", figsize=(9, 8))
plt.grid(which="both", ls="dashed", color="gray")

# Plot con barre d’errore

plt.errorbar(V, epsilon1, yerr=sigmaepsilon1,  fmt='o', color='orange' , label=r'$\epsilon_1$', capsize=3)
#plt.errorbar(V, varepsilon1, yerr=sigmavarepsilon1, fmt='.', color= 'red', label=r'$\bar{\epsilon}_1$', capsize=3)

plt.errorbar(V2, epsilon2, yerr= sigmaepsilon2,  fmt='.',color= 'green', label=r'$\epsilon_2$', capsize=3)

#plt.errorbar(V, varepsilon3, yerr=sigmavarepsilon3, fmt='.', color= 'dodgerblue', label=r'$\bar{\epsilon}_3$', capsize=3)
#plt.errorbar(V, epsilon3,  fmt='o',color='blue', label=r'$\epsilon_4$', capsize=3)


#CON IL RICALCOLO
plt.errorbar(V4, epsilon4, yerr=sigmaepsilon4, fmt='o', color='purple' , label=r'$\epsilon_4$', capsize=3)
#plt.errorbar(V, varepsilon4, fmt='.', color= 'red', label=r'$\bar{\epsilon}_3_1$', capsize=3)


# Etichette e legenda
plt.xlabel("Tensione V [V]")
plt.ylabel("Efficienza ")
plt.title("Efficienza PMT1, PMT2, PMT4")
plt.legend()
plt.show()