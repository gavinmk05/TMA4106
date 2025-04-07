import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Definere konstanter

a = 127  # termisk diffusivitet til gull (mm^2/s)
lengde = 30 # lengde av stang i noe enhet
tid = 1 # tid (s)
punkter = 31 # antall punkter på stangen
mål_tid = 0.1  # Tidspunkt vi vil plotte for (s)


dx = lengde/ (punkter-1) # avstand mellom punktene i stangen (mm)
dt = 0.25*dx**2 / a # tidssteg (s) for stabilitet
t_punkter = int(tid / dt) # antall tidssteg
mål_steg = int(mål_tid / dt)  # Tidssteg som tilsvarer mål_tid


# Initialbetingelser
x = np.linspace(0, lengde, punkter)     # x-koordinater for punktene i stangen
u = np.abs(np.sin(np.pi*x/lengde))             # u(x,0) = sin(pi*x)



fig,ax = plt.subplots()

pcm = ax.pcolormesh([u], cmap=plt.cm.jet, vmin=0, vmax=1)
plt.colorbar(pcm, ax=ax, label='Temperatur (°C)')
ax.set_xlabel('punkter')



# simulasjon, eulers eksplisitt metode
for n in range(t_punkter):  # Tidssteg
    u_ny = u.copy()  # Lag en kopi for å lagre nye verdier
    
    u_ny[0], u_ny[-1] = 0, 0  

    for i in range(1, punkter-1):  # Unngå randpunktene (0 og -1)
        u_ny[i] = u[i] + (a * dt / dx**2) * (u[i+1] - 2*u[i] + u[i-1])  
    

    
    u = u_ny  # Oppdater hele arrayet for neste tidssteg

    print(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnittlig temperatur: {np.mean(u[1:-1]):.3f} celcius")

    if n % 5 == 0:  # Oppdater plott hver 5. iterasjon
        pcm.set_array([u])
        ax.set_title(f"Tid: {(n+1)*dt:.4f} s, Gjennomsnittlig temperatur: {np.mean(u[1:-1]):.3f} °C")
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