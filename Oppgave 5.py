import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Definere konstanter

a = 127  # termisk diffusivitet til gull (mm^2/s)
lengde = 30 # lengde på stangen (mm)
tid = 1 # tid (s)
punkter = 31 # antall punkter i stangen
mål_tid = 0.1  # Tidspunkt vi vil plotte for (s)


dx = lengde / (punkter-1)       # avstand mellom punktene i stangen (mm)
dt = 0.25*dx**2 / a        # tidssteg (s) for stabilitet
t_punkter = int(tid / dt)   # antall tidssteg
mål_steg = int(mål_tid / dt)  # Tidssteg som tilsvarer mål_tid



# Initialbetingelser
x = np.linspace(0, lengde, punkter)     # x-koordinater for punktene i stangen
u = np.abs(np.sin(np.pi*x/lengde))             # u(x,0) = sin(pi*x)

# Bygg systemmatrisen for implisitt Euler
diagonals = [-a * dt / dx**2 * np.ones(punkter-1), 
             (1 + 2 * a * dt / dx**2) * np.ones(punkter), 
             -a * dt / dx**2 * np.ones(punkter+1)]
A = diags(diagonals, [-1, 0, 1], format='csc')

fig, ax = plt.subplots()
pcm = ax.pcolormesh([u], cmap=plt.cm.jet, vmin=0, vmax=1)
plt.colorbar(pcm, ax=ax, label='Temperatur (°C)')
ax.set_xlabel('Punkter')




# simulasjon, eulers implisitt metode

for n in range(t_punkter):    
    u[0], u[-1] = 0, 0    # randbetingelser
    u = spsolve(A, u)    # løse systemet Ax = b for x = u_ny
   

    print(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnittlig temperatur: {np.mean(u[1:-1]):.3f} celcius")

    if n % 5 == 0:
        pcm.set_array([u])
        ax.set_title(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnittlig temperatur: {np.mean(u[1:-1]):.3f} celcius")
        plt.pause(0.1) # pause for å oppdatere grafen

    # Lagre resultatet ved mål_tid
    if n == mål_steg:
        u_på_mål_tid = u.copy()
        print(f"Data lagret for tid: {(n+1)*dt:.4f} s")

plt.figure(figsize=(10, 6))
plt.plot(x, u_på_mål_tid, 'o-', label=f'Tid = {mål_tid:.3f} s')
plt.xlabel("Posisjon (mm)")
plt.ylabel("Temperatur (°C)")
plt.title(f"Temperaturfordeling ved t = {mål_tid:.3f} s")
plt.grid(True)
plt.legend()
plt.show()

print(t_punkter)