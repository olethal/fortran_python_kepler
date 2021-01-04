import matplotlib.pyplot as plt
import numpy as np
import euler
import euler_cromer
import runge_kutta_2
import velocity_verlet
import time as t


h = 10000


#Nmax = 100000000
N = 10000


###############################################
#Bestimmung der Startwerte für die Erde im Schwerefeld der Sonne
#Astronomische Einheit AU in Metern
AU = 149597870700
#Gravitationskonstante
G = 6.67430E-11
#Sonnenmasse
M_sun = 1.98847E30



#Erde

M_earth = 2.99E-6*M_sun
extr = 0.017
a = 1*AU

X1_0 = 0.0
X2_0 = a*(1+extr)
V1_0 = np.sqrt(G*M_sun)*((1-extr)/(a*(1+extr))*(1+M_earth/M_sun))**0.5
V2_0 = 0.0


#Merkur
"""
M_mercury = 1.65E-7*M_sun
extr = 0.207
a = 0.387*AU

X1_0 = 0.0
X2_0 = a*(1+extr)
V1_0 = np.sqrt(G*M_sun)*((1-extr)/(a*(1+extr))*(1+M_mercury/M_sun))**0.5
V2_0 = 0.0
"""







x11, x12 = X1_0, X2_0
x21, x22 = V1_0, V2_0



###############################################

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, constrained_layout=True, figsize=(12,3.5))



ax1.axis("equal")
ax2.axis("equal")
ax3.axis("equal")
ax4.axis("equal")

print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((4, N))
tra= euler.euler(tra, x11, x12, x21, x22, h, N)
print(tra)
ax1.plot(tra[1], tra[2], lw= 0.4)
ax1.set_title("Euler")
np.savetxt("euler.txt", np.transpose(tra))
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %s"%np.std(tra[0]))

print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((4, N))
tra = euler_cromer.euler(tra, x11, x12, x21, x22, h, N)
print(tra)
ax2.plot(tra[1], tra[2], lw= 0.4)
ax2.set_title("Euler Cromer")
np.savetxt("euler_cromer.txt", np.transpose(tra))
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %s"%np.std(tra[0]))




print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((4, N))
tra = runge_kutta_2.rk2(tra, x11, x12, x21, x22, h, N)
print(tra)
ax3.plot(tra[1], tra[2], lw= 0.4)
ax3.set_title("RK2")
np.savetxt("runge_kutta_2.txt", np.transpose(tra))
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %s"%np.std(tra[0]))

print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((4, N))
tra = velocity_verlet.velocity_verlet(tra, x11, x12, x21, x22, h, N)
print(tra)
ax4.plot(tra[1], tra[2], lw= 0.4)
ax4.set_title("Velocity Verlet")
np.savetxt("velocity_verlet.txt", np.transpose(tra))
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %s"%np.std(tra[0]))


fig.suptitle("h=%ss, N=%s, T=%.3ga"%(h, N, h*N/(3600*24*365.25)))



plt.show()


