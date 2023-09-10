import numpy as np
import math as m

# values from table
R = 31
B = 3
rho = 1.225
V0 = 8
omega = 2.61
theta_p = -3 * m.pi / 180
beta = 2 * m.pi / 180
c = 0.5
Cl = 0.5
Cd = 0.01

r = 24.5

# tolerance
eps = 1e-10

# initialize a and a_prime
a = 0
a_prime = 0
a_old = 0
a_prime_old = 0

# loop
while a == 0 or (abs(a - a_old) > eps and abs(a_prime - a_prime_old) > eps):
    phi = np.arctan((1-a_old) * V0 / ((1 + a_prime_old) * omega * r))
    theta = theta_p + beta
    alpha = phi - theta
    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
    F = 2/m.pi * np.arccos(np.exp(- B*(R-r)/(2 * r * np.sin(abs(phi)))))
    sigma = c * B / (2 * m.pi * r)

    CT = (1 - a_old)**2 * Cn * sigma / np.sin(phi)**2

    a_old = a
    a_prime_old = a_prime

    if a_old <= 0.33:
        a = CT / (4 * F * (1-a_old))
    else:
        a = CT / (4 * F * (1 - 1/4 * (5 - 3 * a_old) * a_old))
    a_prime = sigma * Ct * (1 + a_prime_old) / (4 * F * np.sin(phi) * np.cos(phi))

print(f"a ={a:.3f}, a_prime = {a_prime:.3f}, F = {F:.3f}")
