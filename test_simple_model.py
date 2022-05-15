# -*- coding: utf-8 -*-

"""
@author: pmenaq
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def LobDyn(state, t):
    """
    Modelo presentado por garcia en:

    García, Oscar, Harold E. Burkhart, and Ralph L. Amateis. 2011.
    “A Biologically-Consistent Stand Growth Model for Loblolly
    Pine in the Piedmont Physiographic Region, USA.”
    Forest Ecology and Management 262(11): 2035–41.
    """
    global q, p
    H, N, O, W = state
    dHdt = q*(48.76/(H**0.07860) - 0.9271*H)
    dNdH = -1.754*10**-11 * (H**3.642) * (N**2.518) / (p**3.037)
    dOdH = 0.1344*p*H*(1.0 - O)
    dWdH = 0.2474*p*O*H*(N**0.4) + 0.4*(W/N)*dNdH
    return dHdt, dNdH, dOdH, dWdH


t = np.arange(1, 31, 1)
q = 0.0224452084897823
c2 = 15582
p = 1
H0 = 1.3
N0 = 2600
G0 = 0
O0 = 1-(1 - min([N0/c2, 1]))**2.4

state0 = [H0, N0, O0, G0]
statest = odeint(LobDyn, state0, t, full_output=True)[0]

fig ,ax = plt.subplots(2,2,figsize=(8,8))
ax[0][0].set_ylabel("H (m)")
ax[0][1].set_ylabel("N (tress/ha)")
ax[1][0].set_ylabel("O ")
ax[1][1].set_ylabel("G (m2/ha)")
ax[0][0].set_xlabel("H (m)")
ax[0][1].set_xlabel("N (tress/ha)")
ax[1][0].set_xlabel("O ")
ax[1][1].set_xlabel("G (m2/ha)")

for q in np.arange(0.01,0.035,0.001):
    states = odeint(LobDyn, [0.6, 1250, 0.04, 0], t, full_output=True)[0]
    #el area basal se obtiene desde la relacion W = G*H
    #se requiere revisar salida ya que no coincide con
    #implementacion excel.
    ax[0][0].plot(t,states[:,0])
    ax[0][1].plot(t,states[:,1])
    ax[1][0].plot(t,states[:,2])
    ax[1][1].plot(t,states[:,3]/states[:,0])

plt.tight_layout()
