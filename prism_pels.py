r"""

Reference: 

Luigi Cannavacciuolo, Jan Skov Pedersen and Peter Schurtenberger, Langmuir 2002, 18, 7, 2922

----------

Authorship and Verification
----------------------------

* **Author: Andi X**
* **Last Modified by: Nov 2024**
* **Last Reviewed by:** **Date:**
"""

import numpy as np
from numpy import inf

name = "prism_pel"
title = """
	Solution of semi-flexible polymer or polyelectrolyte polymer chains, PRISM,
        """
description = """
	This model should work well for semi-flexible polymer or polyelectrolytes with
	self-avodiance.
	This model is based on built-in flexible_cylinder model. 
              """

category = "plugin"
single = False  

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["length",      "Ang",       1000.0, [0, inf],    "volume", "Length of the flexible cylinder"],
    ["kuhn_length", "Ang",        100.0, [0, inf],    "volume", "Kuhn length of the flexible cylinder"],
    ["radius",      "Ang",         20.0, [0, inf],    "volume", "Radius of the flexible cylinder"],
    ["sld",         "1e-6/Ang^2", 10.72, [-inf, inf], "sld",    "Cylinder scattering length density"],
    ["sld_solvent", "1e-6/Ang^2",  9.42, [-inf, inf], "sld",    "Solvent scattering length density"],
    ["parav",       "",            0.01, [-inf, inf], "",       "Parameter v in PRISM model"  ],
    ["sigma",       "",             140, [0, inf],    "",       "sigma fit parameter from PRISM cq"],
    ["rc",          "",              10, [0, inf],    "",       "rc fit parameter from PRISM cq"],
    ["conc",        "100*w/w",      1.0, [0, 100],    "",       "Concentration of polymer"],
    ["wamw",        "g mol^-1",   40000, [0, inf],    "",       "Weight-averaged molecular weight" ]
    ]
# scale should be number density with unit 1e-7 A^-3 
# pylint: enable=bad-whitespace, line-too-long
source = ["lib/polevl.c", "lib/sas_J1.c", "lib/wrc_cyl.c", "prism_pels.c"]

def random():
    """Return a random parameter set for the model."""
    length = 10**np.random.uniform(2, 6)
    radius = 10**np.random.uniform(1, 3)
    kuhn_length = 10**np.random.uniform(-2, 0)*length
    parav = 10**np.random.uniform(0, 1)
    sigma = 10**np.random.uniform(1, 2)
    rc = 10**np.random.uniform(1, 3)
    pars = dict(
        length=length,
        radius=radius,
        kuhn_length=kuhn_length,
        parav=parav,
        sigma=sigma,
        rc=rc,
    )
    return pars


# Test values are not verified in SASView but have been test externally
# Besides are wrong number, ignore them
tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    [{'length':     1000.0,  # test T1
      'kuhn_length': 100.0,
      'radius':       20.0,
      'sld':           1.0,
      'sld_solvent':   6.3,
      'parav':        0.01,
      'conc':          1.0,
      'wamw':        40000,
      'background':    0.0001,
     }, 0.001, 3509.2187],

    # Additional tests with larger range of parameters
    [{'length':    1000.0,  # test T2
      'kuhn_length': 100.0,
      'radius':       20.0,
      'sld':           1.0,
      'sld_solvent':   6.3,
      'parav':        0.01,
      'conc':          1.0,
      'wamw':        40000,
      'background':    0.0001,
     }, 1.0, 0.000595345],
    [{'length':        10.0,  # test T3
      'kuhn_length': 800.0,
      'radius':        2.0,
      'sld':           6.0,
      'sld_solvent':  12.3,
      'parav':        0.01,
      'conc':          1.0,
      'wamw':        40000,
      'background':    0.001,
     }, 0.1, 1.55228],
    [{'length':        100.0,  # test T4
      'kuhn_length': 800.0,
      'radius':       50.0,
      'sld':           0.1,
      'sld_solvent':   5.1,
      'parav':        0.01,
      'conc':          1.0,
      'wamw':        40000,
      'background':    0.0,
     }, 1.0, 0.000938456]
    ]

