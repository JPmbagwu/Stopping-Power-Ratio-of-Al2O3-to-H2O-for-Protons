#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:14:47 2023

@author: johnpaulmbagwu
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
elementary_charge = 1.60217663e-19  # Coulombs
proton_mass = 1.6726219e-27  # kg
speed_of_light = 299792458  # m/s
avogadro_number = 6.02214076e23  # mol^(-1)

# Material-specific parameters (Let's replace these with accurate values)
z_al2o3 = 13  # Aluminum oxide atomic number
a_al2o3 = 101.96  # Aluminum oxide atomic mass (g/mol)
i_al2o3 = 166.0  # Aluminum oxide mean excitation potential (eV)

z_h2o = 1  # Water atomic number
a_h2o = 18.01528  # Water atomic mass (g/mol)
i_h2o = 75.0  # Water mean excitation potential (eV)

# Proton energies (MeV)
proton_energies = np.linspace(1, 100, 100)

# Let's Calculate stopping power ratios
stopping_power_ratios = []
for energy in proton_energies:
    beta = (energy * 1e6 * elementary_charge) / (proton_mass * speed_of_light ** 2)
    gamma = 1 / np.sqrt(1 - beta**2)
    
    dEdx_al2o3 = (4 * np.pi * elementary_charge**4 * z_al2o3**2) / (proton_mass * speed_of_light**2 * a_al2o3) * (z_al2o3 / i_al2o3) * (avogadro_number / z_al2o3) * (1 / (beta**2) * (0.5 * np.log((2 * proton_mass * speed_of_light**2 * beta**2 * gamma**2) / (i_al2o3**2)) - beta**2))
    dEdx_h2o = (4 * np.pi * elementary_charge**4 * z_h2o**2) / (proton_mass * speed_of_light**2 * a_h2o) * (z_h2o / i_h2o) * (avogadro_number / z_h2o) * (1 / (beta**2) * (0.5 * np.log((2 * proton_mass * speed_of_light**2 * beta**2 * gamma**2) / (i_h2o**2)) - beta**2))
    
    stopping_power_ratio = dEdx_al2o3 / dEdx_h2o
    stopping_power_ratios.append(stopping_power_ratio)

# Let's Plot the stopping power ratio vs proton energy
plt.figure()
plt.plot(proton_energies, stopping_power_ratios, label='Al2O3/H2O Stopping Power Ratio')
plt.xlabel('Proton Energy (MeV)')
plt.ylabel('Stopping Power Ratio')
plt.title('Stopping Power Ratio of Al2O3 to H2O for Protons')
plt.legend()
plt.grid()
plt.show()
