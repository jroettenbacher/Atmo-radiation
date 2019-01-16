# -*- coding: utf-8 -*-
"""
Created on Mon Dec  25 10:46:25 2018

@author: Johannes Roettenbacher

Script for a Monte-Carlo simulation for calculating radiative transfer through a single cloud with variable number 
of layers.
Goal in Task A: Calculate and plot F upward TOA, F downward BOA in dependence of optical thickness tau (C_ext),
single scattering albedo omega and scattering phase function represented by the asymmetry parameter g.
Goal in Task B: Differentiate FBOA into direct transmitted and scattered diffuse radiation
"""
# %%
import math
import random as rn
import numpy as np
import matplotlib.pyplot as plt
# %%
# Starting variables
# cext = 2*10**-10
surface_albedo = 0.05
r = 5*10**-6  # droplet radius
# %%


def layer(tau_layer, angle, fun_layer, fun_omega, fun_g):
    """
    Computes whether a photon is extinct or not and whether it is scattered or absorbed.
    If it is scattered the scatter direction is computed.

    Parameters
    ----
    tau_layer: float
        thickness of layer
    angle: float
        incoming angle of photon in radiant
    fun_layer: integer
        layer photon is currently in. 0 = ground
    fun_omega: float
        defines single scattering albedo
    fun_g: float
        defines asymmetry parameter g
    Returns
    ----
    angle: float
           the angle at which the photon continues its path through the cloud
    fun_layer: integer
               the layer the photon is in after the processes
    """
    tau1 = -math.cos(angle)*math.log(rn.random())
    if tau1 < tau_layer:
        if rn.random() < fun_omega:
            my_scatter = 0.5/fun_g*((1+fun_g**2)-((1-fun_g**2)/(1+fun_g-2*fun_g*rn.random()))**2)
            angle = angle + math.acos(my_scatter)
            if angle > (2*math.pi):
                angle = angle - 2*math.pi
        else:
            fun_layer = None
    if angle < (math.pi/2) or angle > (3*math.pi/4):
        direction = -1
    else:
        direction = 1
    if fun_layer is not None:
        fun_layer = fun_layer + direction
    return angle, fun_layer


# # %%
#
# # test what number of layers is appropriate
# omega = 1.  # single scattering albedo
# g = 0.9  # asymmetry parameter
# tau = 10  # cloud optical thickness
# n_layer = np.array(range(10, 101))
# n_photon = 5000
# delta_tau = (tau / n_layer)  # optical thickness of each layer
# output = np.empty((len(n_layer), 4))
# for i in range(len(n_layer)):
#     print(i)
#     ground_count = 0  # counter for photons reaching the ground
#     toa_up_count = 0  # counter for photons reflected into space
#     in_cloud_count = 0  # counter for photons absorbed in the cloud
#     for n_p in range(1, n_photon+1):
#         theta = 20/180*math.pi
#         photon_in_cloud = True
#         in_layer = n_layer[i]
#         # print("***************************")
#         while photon_in_cloud == True:
#             # print(in_layer)
#             # print(theta/math.pi*180)
#             if in_layer is None:
#                 photon_in_cloud = False
#                 in_cloud_count += 1
#             if in_layer == 0:
#                 if rn.random() > surface_albedo:
#                     photon_in_cloud = False
#                     ground_count += 1
#                 else:
#                     in_layer += 1
#             theta, in_layer = layer(tau_layer=delta_tau[i], angle=theta, fun_layer=in_layer)
#             if in_layer == n_layer[i] + 1:
#                 photon_in_cloud = False
#                 toa_up_count += 1
#     print("F_TOA: %8.6f \nF_BOA: %8.6f \nAbsorbed: %8.6f" % (toa_up_count/n_photon, ground_count/n_photon,
#           in_cloud_count/n_photon))
#     output[i] = (n_layer[i], toa_up_count/n_photon, ground_count/n_photon, in_cloud_count/n_photon)
#     print(output[i])
#
# plt.xlabel("Number of layers")
# # plt.axes(output[:, 0])
# plt.plot(output[:, 1])
# plt.plot(output[:, 2])
# plt.show()  # 50 layers should be enough

# %%


# vary cloud optical thickness, single scattering albedo and asymmetry parameter
tau = np.array(range(1, 21))  # cloud optical thickness
omega = np.arange(0.5, 1.01, 0.1)  # single scattering albedo
g = np.arange(-0.9, 1, 0.1)  # asymmetry parameter
n_layer = 50
n_photon = 5000
delta_tau = (tau / n_layer)  # optical thickness of each layer
plot_data = np.empty((len(tau)*len(omega)*len(g), 6))
h = 0
for i in range(len(tau)):
    for j in range(len(omega)):
        for k in range(len(g)):
            print("tau = %2d, omega = %2.1f, g = %2.1f" % (tau[i], omega[j], g[k]))
            ground_count = 0  # counter for photons reaching the ground
            toa_up_count = 0  # counter for photons reflected into space
            in_cloud_count = 0  # counter for photons absorbed in the cloud
            for n_p in range(1, n_photon+1):
                theta = 20/180*math.pi
                photon_in_cloud = True
                in_layer = n_layer
                # print("***************************")
                while photon_in_cloud == True:
                    # print(in_layer)
                    # print(theta/math.pi*180)
                    if in_layer is None:
                        photon_in_cloud = False
                        in_cloud_count += 1
                    if in_layer == 0:
                        if rn.random() > surface_albedo:
                            photon_in_cloud = False
                            ground_count += 1
                        else:
                            in_layer += 1
                    theta, in_layer = layer(tau_layer=delta_tau[i], angle=theta, fun_layer=in_layer, fun_omega=omega[j],
                                            fun_g=g[k])
                    if in_layer == n_layer + 1:
                        photon_in_cloud = False
                        toa_up_count += 1
            # print("F_TOA: %8.6f \nF_BOA: %8.6f \nAbsorbed: %8.6f" % (toa_up_count/n_photon, ground_count/n_photon,
            #                                                          in_cloud_count/n_photon))
            plot_data[h] = (tau[i], omega[j], g[k], toa_up_count/n_photon, ground_count/n_photon, in_cloud_count/n_photon)
            h += 1
