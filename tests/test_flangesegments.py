
from pyflange.flangesegments import PolynomialLFlangeSegment
from pyflange.bolts import MetricBolt

from math import pi
import numpy as np
import matplotlib.pyplot as plt

m = 1
mm = 0.001*m

kg = 1
t = 1000*kg

s = 1

N = kg*m/s**2
kN = 1000*N

Pa = 1
GPa = 1e9*Pa


fseg = PolynomialLFlangeSegment(

    a = 190*mm,           # distance between inner face of the flange and center of the bolt hole
    b = 140*mm,           # distance between center of the bolt hole and center-line of the shell
    s =  80*mm,           # shell thickness
    t = 170*mm,           # flange thickness
    c = 2*pi/120 * 6.5*m, # shell arc length
    R = 6.5*m,            # shell curvature radius

    Zg = -(700*t + 800*t)*9.81*m/s**2 / 120, # load applied to the flange segment shell at rest
                                             # (normally dead weight of tower + RNA, divided by the number of bolts)

    bolt = MetricBolt("M64", "8.8"),
    Fv = 1200*kN,                            # applied bolt preload

    Do = 66*mm,     # bolt hole diameter
    Dw = 115*mm,    # washer diameter

    gap_height = 1.4*mm,   # maximum longitudinal gap height
    gap_length = 2000*mm)  # longitudinal gap length



Z1 = fseg.shell_force_at_rest
Fs1 = fseg.bolt_force_at_rest

Z2 = fseg.shell_force_at_tensile_ULS
Z2A = fseg._ideal_shell_force_at_tensile_ULS
Z2B = fseg._cantilever_shell_force_at_tensile_ULS
Fs2 = fseg.bolt_force_at_tensile_ULS

Z3 = fseg.shell_force_at_small_displacement
Fs3 = fseg.bolt_force_at_small_displacement

poli = fseg._tensile_force_polynomial
cpoli = fseg._compressive_force_polynomial


print("")

print(f"Z1 = {Z1/kN:.2f} kN")
print(f"Fs1 = {Fs1/kN:.2f} kN")
print("")

print(f"Z0 = {Z2A/kN:.2f} kN")
print(f"Z2_ = {Z2B/kN:.2f} kN")
print(f"Z2 = {Z2/kN:.2f} kN")
print(f"Fs2 = {Fs2/kN:.2f} kN")
print("")

print(f"Z3 = {Z3/kN:.2f} kN")
print(f"Fs3 = {Fs3/kN:.2f} kN")
print("")

print("Polynomial in tensile domain:")
print(poli)

Z = np.linspace(Z1, Z2, 1000)
F = fseg.bolt_axial_force(Z)
# F = np.array(poli(Zi) for Zi in Z])

plt.xlabel("Z [kN]")
plt.ylabel("Fs [kN]")
plt.grid(visible=True)
plt.plot(Z/kN, F/kN)
# plt.plot(np.array([Z1, Z3, Z2])/kN, np.array([Fs1, Fs3, Fs2])/kN)
plt.show()

print("")
print("Polynomial in compressive domain:")
print(cpoli)

Z = np.linspace(cpoli.domain[0], cpoli.domain[1], 1000)
F = cpoli(Z)
# F = np.array(poli(Zi) for Zi in Z])

plt.xlabel("Z [kN]")
plt.ylabel("Fs [kN]")
plt.grid(visible=True)
plt.plot(Z/kN, F/kN)
# plt.plot(np.array([Z1, Z3, Z2])/kN, np.array([Fs1, Fs3, Fs2])/kN)
plt.show()


print("")
print("Total polynomial:")

Z = np.linspace(1.2*cpoli.domain[0], 1.2*poli.domain[1], 2000)
F = fseg.bolt_axial_force(Z)
# F = np.array(poli(Zi) for Zi in Z])

plt.xlabel("Z [kN]")
plt.ylabel("Fs [kN]")
plt.grid(visible=True)
plt.plot(Z/kN, F/kN)
# plt.plot(np.array([Z1, Z3, Z2])/kN, np.array([Fs1, Fs3, Fs2])/kN)
plt.show()
