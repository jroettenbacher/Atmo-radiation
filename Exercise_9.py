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
# %%
# Starting variables
# cext = 2*10**-10
surface_albedo = 0.05
omega = 1.  # single scattering albedo
g = 0.9  # asymmetry parameter
r = 5*10**-6  # droplet radius
tau = 1  # cloud optical thickness
n_layer = 20
n_photon = 5
ground_count = 0  # counter for photons reaching the ground
toa_up_count = 0  # counter for photons reflected into space
in_cloud_count = 0  # counter for photons absorbed in the cloud
delta_tau = tau / n_layer  # optical thickness of each layer
# %%


def layer(tau_layer, angle, fun_layer):
    """
    Computes whether a photon is extinct or not and if it is scattered or absorbed.
    If it is scattered the scatter direction is computed.

    Parameters
    ----
    tau_layer: float
               thickness of layer
    angle: float
           incoming angle of photon in radiant
    fun_layer: integer
               layer photon is currently in. 0 = ground

    Returns
    ----
    angle: float
           the angle at which the photon continues its path through the cloud
    fun_layer: integer
               the layer the photon is in after the processes
    """
    tau1 = -math.cos(angle)*math.log(rn.random())
    if tau1 < tau_layer:
        if rn.random() < omega:
            my_scatter = 0.5/g*((1+g**2)-((1-g**2)/(1+g-2*g*rn.random()))**2)
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


# %%


for n_p in range(1, n_photon+1):
    theta = 20/180*math.pi
    photon_in_cloud = True
    in_layer = n_layer
    print("***************************")
    while photon_in_cloud == True:
        print(in_layer)
        print(theta/math.pi*180)
        if in_layer is None:
            photon_in_cloud = False
            in_cloud_count += 1
        if in_layer == 0:
            if rn.random() > surface_albedo:
                photon_in_cloud = False
                ground_count += 1
            else:
                in_layer += 1
        theta, in_layer = layer(tau_layer=delta_tau, angle=theta, fun_layer=in_layer)
        if in_layer == n_layer + 1:
            photon_in_cloud = False
            toa_up_count += 1
print("F_TOA: %8.6f \nF_BOA: %8.6f \nAbsorbed: %8.6f" % (toa_up_count/n_photon, ground_count/n_photon, in_cloud_count/n_photon))
