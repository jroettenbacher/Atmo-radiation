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
import numpy as np
import pandas


# %% layer function


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
    fun_layer: float
        layer photon is currently in. 0 = ground
    fun_omega: float
        defines single scattering albedo
    fun_g: float
        defines asymmetry parameter g
    Return
    ----
    angle: float
           the angle at which the photon continues its path through the cloud
    fun_layer: integer
               the layer the photon is in after the processes
    scattered: bool
               variable which tells whether photon has been scattered or not
    """
    tau1 = abs(np.cos(angle) * np.log(np.random.random()))
    if tau1 < tau_layer:
        scattered = True
        if np.random.random() < fun_omega:
            my_scatter = (0.5 / fun_g) * ((1 + fun_g ** 2) - ((1 - fun_g ** 2) / (1 + fun_g - 2 * fun_g * np.random.random())) ** 2)
            angle = angle + np.arccos(my_scatter)
            if angle > (2 * np.pi):
                angle = angle - 2 * np.pi
        else:
            fun_layer = None
    else:
        scattered = False
    if angle < (np.pi / 2) or angle > (3 * np.pi / 2):
        direction = -1
    else:
        direction = 1
    if fun_layer is not None:
        fun_layer = fun_layer + direction
    return angle, fun_layer, scattered


# %% test what layer thickness is appropriate

surface_albedo = 0.05
omega = 1  # single scattering albedo
g = 0.9  # asymmetry parameter
tau = 10  # cloud optical thickness
n_layer = np.array(range(10, 101))
n_photon = 5000
delta_tau = (tau / n_layer)  # optical thickness of each layer
output = np.empty((len(n_layer), 3))
for i in range(len(n_layer)):
    print(i / len(n_layer))
    ground_count = 0  # counter for photons reaching the ground
    toa_up_count = 0  # counter for photons reflected into space
    in_cloud_count = 0  # counter for photons absorbed in the cloud
    for n_p in range(1, n_photon + 1):
        theta = 20 / 180 * np.pi
        photon_in_cloud = True
        in_layer = n_layer[i]
        while photon_in_cloud is True:
            if in_layer == 0:
                if np.random.random() > surface_albedo:
                    photon_in_cloud = False
                    ground_count += 1
                else:
                    in_layer += 1
            if in_layer is None:
                photon_in_cloud = False
                in_cloud_count += 1
            theta, in_layer, scattered = layer(tau_layer=delta_tau[i], angle=theta, fun_layer=in_layer,
                                               fun_omega=1, fun_g=0.9)
            if in_layer == n_layer[i] + 1:
                photon_in_cloud = False
                toa_up_count += 1
    output[i] = (n_layer[i], toa_up_count / n_photon, ground_count / n_photon)

df_output = pandas.DataFrame({"n_layer": output[:, 0], "FTOA": output[:, 1], "FBOA": output[:, 2]})
df_output.to_csv(
    "C:/Users/Johannes/Nextcloud/Studium/Atmosphaerische_Strahlung/Atmo_Strahlung/delta_tau_layer_depend.csv",
    index=False)


# %% monte_carlo method


def monte_carlo(tau, omega, g, n_photon=5000, theta_in=20, delta_tau=0.2):
    """
    Runs a MonteCarlo simulation with the given parameters.
    Parameters
    ----
    tau: integer
        optical thickness of cloud
    omega: float
        single scattering albedo
    g: float
        asymmetry parameter of the phase function
    n_photon: integer
        number of photons with which the model should run
    theta_in: float
        Incoming solar zenith angle for photon hitting cloud top in degree
    delta_tau: float
        thickness of each individual cloud layer, default = 0.2
    Returns
    ----
    Input parameters tau, omega, g and theta_in. Additionally the fraction of FTOA, FBOA,
    FBOA_diffuse1/2, FBOA_direct1/2 and in cloud absorbed photons is returned.
    """
    ground_count = 0  # counter for photons reaching the ground
    toa_up_count = 0  # counter for photons reflected into space
    in_cloud_count = 0  # counter for photons absorbed in the cloud
    boa_direct1 = 0  # counter for photons directly transmitted to surface, method 1
    boa_diffuse1 = 0  # counter for photons diffusely transmitted to surface, method 1
    boa_direct2 = 0  # counter for photons directly transmitted to surface, method 2
    boa_diffuse2 = 0  # counter for photons diffusely transmitted to surface, method 2
    n_layer = tau / delta_tau  # number of layers in cloud
    surface_albedo = 0.05
    for n_p in range(1, n_photon + 1):
        theta = theta_in / 180 * np.pi
        photon_in_cloud = True
        in_layer = n_layer
        diffuse = False
        while photon_in_cloud is True:
            if in_layer == 0:
                if np.random.random() > surface_albedo:
                    photon_in_cloud = False
                    ground_count += 1
                    if diffuse is True:
                        boa_diffuse1 += 1
                    else:
                        boa_direct1 += 1
                    if (theta_in - (1 / 180 * np.pi)) < theta < (theta_in + (1 / 180 * np.pi)):
                        boa_direct2 += 1
                    else:
                        boa_diffuse2 += 1
                else:
                    theta = theta + np.pi
                    if theta > (2 * np.pi):
                        theta = theta - (2 * np.pi)
                    in_layer += 1
            if in_layer is None:
                photon_in_cloud = False
                in_cloud_count += 1
            theta, in_layer, scattered = layer(tau_layer=delta_tau, angle=theta, fun_layer=in_layer, fun_omega=omega,
                                               fun_g=g)
            if scattered is True:
                diffuse = True
            if in_layer == n_layer + 1:
                photon_in_cloud = False
                toa_up_count += 1
    return tau, omega, g, toa_up_count / n_photon, ground_count / n_photon, in_cloud_count / n_photon, \
           boa_diffuse1 / (ground_count + 1), boa_direct1 / (ground_count + 1), \
           boa_diffuse2 / (ground_count + 1), boa_direct2 / (ground_count + 1), theta_in


# %% vary number of photons
n_photon = np.arange(1, 20002, 5000)
plot_data1 = np.empty((len(n_photon), 11))
h = 0
for x in range(len(n_photon)):
    plot_data1[h] = monte_carlo(tau=10, omega=0.999, g=0.9, n_photon=n_photon[x])
    h += 1
    print(h / (len(n_photon)))
df = pandas.DataFrame({"tau": plot_data1[:, 0], "omega": plot_data1[:, 1], "g": plot_data1[:, 2],
                       "FTOA": plot_data1[:, 3], "FBOA": plot_data1[:, 4], "Absorbed": plot_data1[:, 5]})
df.to_csv("C:/Users/Johannes/Nextcloud/Studium/Atmosphaerische_Strahlung/Atmo_Strahlung/ex09_plot_data1.csv",
          index=False)

# %% vary cloud optical thickness, single scattering albedo and asymmetry parameter
tau = np.array(range(1, 11))  # cloud optical thickness
omega = np.array([0.9, 0.99, 1])  # single scattering albedo
g = np.array([-0.9, -0.5, -0.01, 0.5, 0.9])  # asymmetry parameter
plot_data2 = np.empty((len(tau) * len(omega) * len(g), 11))
h = 0
for i in range(len(tau)):
    for j in range(len(omega)):
        for k in range(len(g)):
            plot_data2[h] = monte_carlo(tau=tau[i], omega=omega[j], g=g[k], n_photon=5000)
            h += 1
            print(h/(len(tau) * len(omega) * len(g)))

df1 = pandas.DataFrame({"tau": plot_data2[:, 0], "omega": plot_data2[:, 1], "g": plot_data2[:, 2],
                        "FTOA": plot_data2[:, 3], "FBOA": plot_data2[:, 4], "Absorbed": plot_data2[:, 5],
                        "BOA_diffuse1": plot_data2[:, 6], "BOA_direct1": plot_data2[:, 7],
                        "BOA_diffuse2": plot_data2[:, 8], "BOA_direct2": plot_data2[:, 9]})
df1.to_csv(
    "C:/Users/Johannes/Nextcloud/Studium/Atmosphaerische_Strahlung/Atmo_Strahlung/ex09_plot_data2_B.csv",
    index=False)

# %% vary cloud optical thickness, solar zenith angle and asymmetry parameter
tau = np.array(range(1, 11))  # cloud optical thickness
g = np.array([-0.9, -0.5, -0.01, 0.5, 0.9])  # asymmetry parameter
theta = np.arange(1, 82, 20)
plot_data3 = np.empty((len(tau) * len(theta) * len(g), 11))
h = 0
for i in range(len(tau)):
    for j in range(len(g)):
        for k in range(len(theta)):
            plot_data3[h] = monte_carlo(tau=tau[i], omega=1, g=g[j], n_photon=5000,
                                        theta_in=theta[k])
            h += 1
            print(h/(len(tau) * len(theta) * len(g)))

df2 = pandas.DataFrame({"tau": plot_data3[:, 0], "omega": plot_data3[:, 1], "g": plot_data3[:, 2],
                        "FTOA": plot_data3[:, 3], "FBOA": plot_data3[:, 4], "Absorbed": plot_data3[:, 5],
                        "BOA_diffuse1": plot_data3[:, 6], "BOA_direct1": plot_data3[:, 7],
                        "BOA_diffuse2": plot_data3[:, 8], "BOA_direct2": plot_data3[:, 9],
                        "theta": plot_data3[:, 10]})
df2.to_csv(
    "C:/Users/Johannes/Nextcloud/Studium/Atmosphaerische_Strahlung/Atmo_Strahlung/ex09_plot_data3_B.csv",
    index=False)
