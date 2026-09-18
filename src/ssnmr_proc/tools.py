# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 15:27:29 2023
Usefull NMR tools.

@author: Marcos de Oliveira Jr.
"""
"""  """
import numpy as np
import re
import csv
from scipy.constants import pi, hbar, mu_0, physical_constants

# Nuclear magneton (J/T)
mu_N = physical_constants["nuclear magneton"][0]

# ---------------------------------------------------------------------
# Nuclear parameters
# spin: I
# g: nuclear g factor
# abundance: natural abundance (fraction)
# ---------------------------------------------------------------------
NUCLEAR_DATA = {
    "1H":   {"spin": 0.5, "g": 5.5856946893, "abundance": 0.999885},
    "2H":   {"spin": 1.0, "g": 0.8574382338, "abundance": 0.000115},
    "7Li":  {"spin": 1.5, "g": 2.170951,     "abundance": 0.9241},
    "11B":  {"spin": 1.5, "g": 1.7924326,    "abundance": 0.801},
    "13C":  {"spin": 0.5, "g": 1.404825,     "abundance": 0.0107},
    "15N":  {"spin": 0.5, "g": -0.566378,    "abundance": 0.00364},
    "17O":  {"spin": 2.5, "g": -0.757516,    "abundance": 0.00038},
    "19F":  {"spin": 0.5, "g": 5.257736,     "abundance": 1.0000},
    "23Na": {"spin": 1.5, "g": 1.478348,     "abundance": 1.0000},
    "25Mg": {"spin": 2.5, "g": -0.342181,    "abundance": 0.1000},
    "27Al": {"spin": 2.5, "g": 1.456602,     "abundance": 1.0000},
    "29Si": {"spin": 0.5, "g": -1.11058,     "abundance": 0.04685},
    "31P":  {"spin": 0.5, "g": 2.26320,      "abundance": 1.0000},
    "45Sc": {"spin": 3.5, "g": 1.3590,       "abundance": 1.0000},
    "51V":  {"spin": 3.5, "g": 1.47106,      "abundance": 0.9975},
    "67Zn": {"spin": 2.5, "g": 0.350192,     "abundance": 0.0410},
    "69Ga": {"spin": 1.5, "g": 1.34439,      "abundance": 0.6011},
    "71Ga": {"spin": 1.5, "g": 1.70818,      "abundance": 0.3989},
    "111Cd":{"spin": 0.5, "g": -1.18977,     "abundance": 0.1280},
    "113Cd":{"spin": 0.5, "g": -1.24460,     "abundance": 0.1222},
}

#%% Calculate Diso and SOQE from MQMAS position
def mqmascalc(X1, X2, Spin, freq):
    """
    Calculate delta_cs and SOQE values based on input parameters.
    
    Parameters:
    X1, X2 (float): Input values
    Nucleus (str): Nucleus name
    freq (float): Frequency value
    
    Returns:
    tuple: (delta_cs, SOQE) or (None, None) if inconsistent input
    """
    # Get nuclear spin (this requires implementation of nucspin function)
    I = Spin
    
    # Calculate delta_cs
    delta_cs = (X2 * 20.0 / 34.0 + X1) * 34.0 / 54.0
    
    # Check for inconsistent input
    if (delta_cs - X2) < 0.0:
        print('Inconsistent input!')
        return None, None
    
    # Calculate SOQE using numpy operations
    SOQE = (delta_cs - X2) * np.power(freq, 2) * 1e6 * 10.0 * np.power((2.0 * I * (2.0 * I - 1.0)), 2) / (I * (I + 1.0) - 0.75) / 3.0
    SOQE = np.sqrt(SOQE) / 1e6

    return delta_cs, SOQE

def center_of_mass(x,y):
    cg=0
    for i in range(y.size):
        cg = cg + x[i]*y[i]
    return cg/y.sum()

def read_distances(filename):

    for encoding in ("utf-8", "cp1252", "latin1"):

        try:

            with open(filename, newline="", encoding=encoding) as f:

                reader = csv.reader(f, delimiter=";")

                next(reader)

                data = {}

                current_site = None

                for row in reader:

                    if not row:
                        continue

                    atom1 = row[0].strip()
                    atom2 = row[1].strip()
                    mult  = row[2].strip()
                    dist  = row[3].strip()

                    if atom1:

                        current_site = atom1

                        data[current_site] = {
                            "neighbors": [],
                            "distances": [],
                            "multiplicities": []
                        }

                    data[current_site]["neighbors"].append(atom2)
                    data[current_site]["multiplicities"].append(
                        int(mult.rstrip("x"))
                    )
                    data[current_site]["distances"].append(
                        float(re.sub(r"\(.*\)", "", dist))
                    )

                break

        except UnicodeDecodeError:
            continue

    else:
        raise UnicodeDecodeError(
            "Could not determine the file encoding."
        )

    for site in data:

        data[site]["distances"] = np.asarray(
            data[site]["distances"], dtype=float
        )

        data[site]["multiplicities"] = np.asarray(
            data[site]["multiplicities"], dtype=int
        )

    return data

def M2(nucleus1, nucleus2, r, n=None):
    """
    Calculate the heteronuclear dipolar second moment.

    Parameters
    ----------
    nucleus1 : str
        Observed nucleus (e.g. '113Cd').
    nucleus2 : str
        Coupled nucleus (e.g. '1H').
    r : array_like
        Distances (Å).
    n : array_like, optional
        Multiplicity of each distance.

    Returns
    -------
    m2 : float
        Second moment (rad²/s²)
    summatory : float
        Σ n_i / r_i^6
    """

    try:
        iso1 = NUCLEAR_DATA[nucleus1]
        iso2 = NUCLEAR_DATA[nucleus2]
    except KeyError as exc:
        raise ValueError(f"Unknown isotope: {exc.args[0]}")

    gamma1 = iso1["g"] * mu_N / hbar
    gamma2 = iso2["g"] * mu_N / hbar

    I = iso2["spin"]
    abundance = iso2["abundance"]

    r = np.asarray(r, dtype=float) * 1e-10

    if n is None:
        n = np.ones_like(r)
    else:
        n = np.asarray(n, dtype=float)

    summatory = np.sum(n / r**6)

    if nucleus1 == nucleus2:
        prefactor = (
            abundance
            * (3/5)
            * (mu_0/(4*pi))**2
            * I*(I+1)
            * gamma1**4
            * hbar**2
        )
    else:
        prefactor = (
            abundance
            * (4/15)
            * (mu_0/(4*pi))**2
            * I*(I+1)
            * gamma1**2
            * gamma2**2
            * hbar**2
        )

    m2 = prefactor * summatory
    
    # Effective heteronuclear dipolar coupling
    D_eff_rad = np.sqrt(3 * m2 / (I * (I + 1)))   # rad/s
    D_eff_hz = D_eff_rad / (2 * pi)               # Hz

    return m2, D_eff_hz

    