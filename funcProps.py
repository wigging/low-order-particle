"""
Functions to determine thermal properties of woody biomass and char.

heatcap - heat capacity, Cp, kJ/(kg*K)
thermalcond - thermal conductivity, k, W/(m*K)
"""

# Heat Capacity from Wood Handbook 2010
# -----------------------------------------------------------------------------

def heatcap(x, T):
    """
    Calculate heat capacity of wood at temperature and moisture content.

    Example:
        cp = heatcap(12, 300)
    Inputs:
        x = moisture content, %
        T = temperature, K
    Output:
        cp_wet = heat capacity wet wood, kJ/(kg*K)

    Reference:
        Glass and Zelinka, 2010. Wood Handbook, Ch. 4, pp. 1-19.
    """

    cpw = 4.18  # heat capacity of water, kJ/(kg*K)

    # coefficients for adjustment factor Ac
    b1 = -0.06191
    b2 = 2.36e-4
    b3 = -1.33e-4

    # adjustment factor for additional energy in wood-water bond, Eq. 4-18
    Ac = x*(b1 + b2*T + b3*x)

    # heat capacity of dry wood, Eq. 4-16a, kJ/(kg*K)
    cp_dry = 0.1031 + 0.003867*T

    # heat capacity of wood that contains water, Eq. 4-17, kJ/(kg*K)
    cp_wet = (cp_dry + cpw*x/100) / (1 + x/100) + Ac

    return cp_wet


# Thermal Conductivity from Wood Handbook 2010
# -----------------------------------------------------------------------------

def thermalcond(x, So, Gb):
    """
    Calculate thermal conductivity of wood at moisture content, volumetric
    shrinkage, and basic specific gravity.

    Example:
        k = thermalcond(12, 12.3, 0.54)
    Inputs:
        x = moisture content, %
        So = volumetric shrinkage, Table 4-3, %
        Gb = basic specific gravity, Table 4-7 or Table 5-3
    Outputs:
        k = thermal conductivity, W/(m*k)

    Reference:
        Glass and Zelinka, 2010. Wood Handbook, Ch. 4, pp. 1-19.
    """

    mcfs = 30   # fiber staturation point estimate, %

    # shrinkage from green to final moisture content, Eq. 4-7, %
    Sx = So*(1 - x/mcfs)

    # specific gravity based on volume at given moisture content, Eq. 4-9
    Gx = Gb / (1 - Sx/100)

    # thermal conductivity, Eq. 4-15, W/(m*K)
    A = 0.01864
    B = 0.1941
    C = 0.004064
    k = Gx*(B + C*x) + A

    return k
