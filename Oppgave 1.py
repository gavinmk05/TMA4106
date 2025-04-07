import numpy as np
import matplotlib.pyplot as plt


h = 0.1

x = (3/2)

def f_derivert(x):
    return np.exp(x)

def f_derivert_approksimert(x,h):
    liste_h = []
    approksimert = []
    differanse = []

    for i in range(0,17):
        derivert = (np.exp(x+h)-np.exp(x))/(h)
        approksimert.append(derivert)
        differanse.append(abs(derivert-f_derivert(x)))
        liste_h.append(h)
        h = h/10
    return liste_h, approksimert, differanse

liste_h, approksimert, differanse = f_derivert_approksimert(x,h)

vanlig_h = [float(x) for x in liste_h] 
vanlig_approksimert = [float(x) for x in approksimert]
vanlig_differanse = [float(x) for x in differanse]

print("Liste med h: ", vanlig_h)
print("Liste med approksimert verdi: ", vanlig_approksimert)
print("Liste med differanse: ", vanlig_differanse)

# Første vindu: Approksimert derivert
plt.figure(1, figsize=(6, 4))
plt.plot(liste_h, approksimert, 'o', linestyle='', label='Approksimert')
plt.title('Approksimert derivert')
plt.xlabel('h')
plt.ylabel('Approksimert verdi')
plt.xscale('log')
plt.legend(loc ='lower right')
plt.grid(True)
plt.savefig('approksimert_derivert.pdf')

# Andre vindu: Differanse
plt.figure(2, figsize=(6, 4))
plt.plot(liste_h, differanse, 'o', linestyle='', label='Differanse')
plt.title('Differanse fra eksakt verdi')
plt.xlabel('h')
plt.ylabel('Absolutt differanse')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc ='lower right')
plt.grid(True)
plt.savefig('differanse.pdf')

# Vis begge vinduene
plt.show()
    

    



