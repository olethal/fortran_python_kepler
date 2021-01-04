python3 -m numpy.f2py -c euler.f90 -m euler --debug-capi
python3 -m numpy.f2py -c euler_cromer.f90 -m euler_cromer --debug-capi
python3 -m numpy.f2py -c runge_kutta_2.f90 -m rk2 --debug-capi
python3 -m numpy.f2py -c velocity_verlet.f90 -m velocity_verlet --debug-capi
