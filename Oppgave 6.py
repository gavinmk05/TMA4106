import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Konstanter
a = 127  # Termisk diffusivitet til gull (mm^2/s)
lengde = 30  # Lengde på stangen (mm)
tid = 1  # Simuleringstid (s)
punkter = 31  # Antall punkter i stangen
mål_tid = 0.1  # Tidspunkt vi vil plotte for (s)



# Diskretisering
dx = lengde /(punkter - 1)  # Avstand mellom punkter (mm)
dt = 0.25*dx**2 / a  # Tidssteg
t_punkter = int(tid / dt)  # Antall tidssteg
mål_steg = int(mål_tid / dt)  # Tidssteg som tilsvarer mål_tid

# Initialbetingelser
x = np.linspace(0, lengde, punkter)
u = np.abs(np.sin(np.pi * x/lengde))  # u(x,0) = |sin(πx/L)|

# Crank-Nicolson koeffisienter
alpha = a * dt / (2 * dx**2)

# Bygg matriser
# Venstre side (LHS): I + θ·α·A
diagonals_LHS = [
    - alpha * np.ones(punkter - 1),
    (1 + 2 * alpha) * np.ones(punkter),
    -alpha * np.ones(punkter - 1)
]
LHS = diags(diagonals_LHS, [-1, 0, 1], format='csc')

# Høyre side (RHS): I - (1-θ)·α·A
diagonals_RHS = [
    alpha * np.ones(punkter - 1),
    (1 - 2 * alpha) * np.ones(punkter),
    alpha * np.ones(punkter - 1)
]
RHS = diags(diagonals_RHS, [-1, 0, 1], format='csc')

# Visualisering
fig, ax = plt.subplots()
pcm = ax.pcolormesh([u], cmap=plt.cm.jet, vmin=0, vmax=1)
plt.colorbar(pcm, ax=ax, label='Temperatur (°C)')
ax.set_xlabel('Punkter')

# Simulering (Crank-Nicolson)
for n in range(t_punkter):
    
    # Randbetingelser
    u[0], u[-1] = 0, 0
    
    # Beregn høyre side: RHS @ u_n
    b = RHS @ u
    
    # Løs LHS @ u_{n+1} = b
    u = spsolve(LHS, b)
    
    print(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnittlig temperatur: {np.mean(u[1:-1]):.3f} celcius")


    # Oppdater plot
    if n % 5 == 0:  # Oppdater plot hvert 10. steg
        pcm.set_array([u])
        ax.set_title(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnitt: {np.mean(u[1:-1]):.3f} °C")
        plt.pause(0.1)

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