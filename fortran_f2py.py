import matplotlib.pyplot as plt
import numpy as np
import euler
import euler_cromer
import runge_kutta_2
import velocity_verlet
import time as t



start_time_fortran = t.time()


h = 1000
#Nmax = 100000000
N = 1000000
###############################################
#Bestimmung der Startwerte für die Erde im Schwerefeld der Sonne
#Astronomische Einheit AU in Metern
AU = 149597870700
AU = 0.467*AU
extr = 0.206
M_sun = 1.98847E30

#M_earth = 2.99E-6*M_sun
M_earth = 1.65E-7*M_sun

G = 6.67430E-11
X1_0 = 0.0
X2_0 = AU
V1_0 = np.sqrt(M_sun*G/AU)-np.sqrt(2*M_sun*G*(1/(AU*np.sqrt(1-extr**2))-1/AU))
V2_0 = 0.0









x11, x12 = X1_0, X2_0
x21, x22 = V1_0, V2_0



###############################################

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)#, sharex=True, sharey=True)


print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((3, N))

#x11, x12, x21, x22 = 0.0, 150E9, 3E4, 0.0
tra= euler.euler(tra, x11, x12, x21, x22, h, N)
print(tra)
ax1.plot(tra[1], tra[2], lw= 0.4)
ax1.set_title("Euler")


print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((3, N))
tra = euler_cromer.euler(tra, x11, x12, x21, x22, h, N)
print(tra)
ax2.plot(tra[1], tra[2], lw= 0.4)
ax2.set_title("Euler Cromer")


print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((3, N))
tra = runge_kutta_2.rk2(tra, x11, x12, x21, x22, h, N)
print(tra)
ax3.plot(tra[1], tra[2], lw= 0.4)
ax3.set_title("RK2")


print("%s Umläufe in Zeitintervall"%(h*N/31536000))
tra = np.zeros((3, N))
tra = velocity_verlet.velocity_verlet(tra, x11, x12, x21, x22, h, N)
print(tra)
ax4.plot(tra[1], tra[2], lw= 0.4)
ax4.set_title("Velocity Verlet")


fig.suptitle("test")

end_time_fortran = t.time()
plt.show()




start_time_python = t.time()
#euler-cromer zum vergleich in reinem python geschrieben
v0 = np.array((V1_0, V2_0))
x0 = np.array((X1_0, X2_0))
tra = np.zeros((3, N))
x = x0
v = v0
for i in range(0, N):
    x += h*v
    v -= h*M_sun*G*x/(np.dot(x,x)**(3.0/2.0))
    tra[0][i] = h*i
    tra[1][i] = x[0]
    tra[2][i] = x[1]
end_time_python = t.time()




plt.plot(tra[1], tra[2], lw = 0.5)
plt.show()

print("it took fortran modules a total of %s seconds to perform %s iterations of various numerical integration routines" % ((end_time_fortran-start_time_fortran), N*4))

print("it took python a total of %s seconds to perform %s iterations of one numerical integration routines" % ((end_time_python-start_time_python), N))


speed = 4*(end_time_python-start_time_python)/(end_time_fortran-start_time_fortran)
print("speed factor is about: %s" % speed)

