
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# parameters dereived from the BEM calculations
V_r = 11.177
rho_air = 1.225
D = 89.17*2
A_rotor = np.pi*(D/2)**2
C_P = 0.498
P_0 = 10.64e6
A = 9
k = 1.9

# =============================================================================
# Define the power function P(V)
def P(V):
    if V < V_r:
        return 1 / 2 * rho_air * C_P * A_rotor * V ** 3
    else:
        return P_0


# Define the Weibull PDF f(v)
def f(v):
    return k/A * (v/A)**(k-1) * np.exp(-(v/A)**k)

# Compute the AEP using numerical integration
def P_f(v):
    return P(v) * f(v)
# =============================================================================


# =============================================================================
# Wind turbine stopped at V0 = 25

V_cut_off_25 = 25

def P_f_new(v):
    if v <= V_cut_off_25:
        return P(v) * f(v)
    else:
        return 0.0

AEP_25_mps, _ = quad(P_f_new, 4, 25)
AEP_25_mps *= 8760/1000000  # Convert to MWh

print(f'Annual Energy Production at {V_cut_off_25} m/s: {AEP_25_mps:.2f} MWh')
# =============================================================================



# =============================================================================
# Wind turbine stopped at V0 = 20

V_cut_off_20 = 20

def P_f_new(v):
    if v <= V_cut_off_20:
        return P(v) * f(v)
    else:
        return 0.0

AEP_20_mps, _ = quad(P_f_new, 4, 20)
AEP_20_mps *= 8760/1000000  # Convert to MWh

print(f'Annual Energy Production at {V_cut_off_20} m/s: {AEP_20_mps:.2f} MWh')
# =============================================================================




# =============================================================================
# Energy loss

loss_MWh = AEP_25_mps - AEP_20_mps
loss_percentage = (loss_MWh / AEP_25_mps) * 100

print(f'Energy lost by stopping at {V_cut_off_20} m/s: {loss_MWh:.2f} MWh')
print(f'Energy loss as a percentage: {loss_percentage:.2f}%')
# =============================================================================


# =============================================================================
# plot

# Plot the power curve
V = np.linspace(4, 25, 1000)
P_V = np.array([P(v) for v in V])
plt.plot(V, P_V, label='Power curve')

AEP, _ = quad(P_f, 4, 25)
AEP *= 8760/1000000 # convert to Wh
plt.title(f'Annual Energy Production = {AEP:.2f} Wh')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Power (W)')
plt.legend()
plt.show()


# Compute the Weibull PDFs for both cutoffs
V_cut_off_25 = 25
f_25 = [f(v) if v <= V_cut_off_25 else 0.0 for v in V]

V_cut_off_20 = 20
f_20 = [f(v) if v <= V_cut_off_20 else 0.0 for v in V]

# Plot the Weibull PDFs
plt.plot(V, f_25, label=f'Weibull PDF (V_cut_off={V_cut_off_25} m/s)')
plt.plot(V, f_20, label=f'Weibull PDF (V_cut_off={V_cut_off_20} m/s)')

plt.title('Weibull PDFs at Different Wind Speed Cutoffs')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Probability density')
plt.legend()
plt.show()













