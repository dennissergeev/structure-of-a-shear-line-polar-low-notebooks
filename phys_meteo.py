# -*- coding: utf-8 -*-
"""
@author: D. E. Sergeev
"""
import numpy as np

#
# Physical parameters and constants for Earth's atmosphere
#
p0       = 1e5            # Reference pressure                      [Pa]
T0       = 273.           # Melting point,                          [K]
T27      = 300.           # Standard temperature                    [K]
day_s    = 86400.         # Day length                              [s]
omg      = 2.*np.pi/day_s # Earth angular velocity,                 [s^-1]
g        = 9.8            # Gravity acceleration,                   [m s^-2]
cp       = 1004.5         # Air heat capacity at constant pressure, [J kg^-1 K^-1]
L        = 2.501e6        # latent heat of vaporization             [J kg^-1]
rho_s    = 1.225          # Air density at the surface layer,       [kg m^-3]
Rd       = 287.04         # Specific gas constant for dry air,      [J kg^-1 K^-1]
kappa    = Rd/cp          # 
Rv       = 461.5          # Specific gas constant for water vapour
eps      = Rd/Rv          #
r_earth  = 6.371e6        # Earth radius [m]
r_sun    = 7.05e8         # Sun radius [m]
S0       = 1360.          # Solar constant                          [W m^-2]
stefan_boltzmann = 5.673e-8 # Stefan-Boltzmann constant             [J m^-2 s^-1 K^-4]


def lon2dist(lon,lat,R=r_earth):
    return R*np.cos(np.radians(lat))*np.radians(lon)

def lat2dist(lat,R=r_earth):
    return R*np.radians(lat)

def uv2wdir(u, v):
    return 180.+180./np.pi*np.arctan2(u, v)

def uv2wspd(u, v):
    return np.sqrt(u**2 + v**2)

def wswd2uv(ws, wd):
    u = ws*np.cos(np.deg2rad(90-wd))
    v = ws*np.sin(np.deg2rad(90-wd))
    return u, v

def deg2name(num_deg, nbins=16):
    assert type(num_deg) != np.str, 'Input cannot be of string type'
    assert nbins==16 or nbins==8 or nbins==4, 'Number of bins must be 4, 8 or 16'
    db = 16//nbins
    deg_lev = np.linspace(0,360,nbins+1)
    deg_bound = deg_lev[1::] - deg_lev[1]/2.
    compass = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'][::db]

    if not hasattr(num_deg,'__iter__'):
        for j in range(len(deg_bound)):
            if deg_bound[j-1] < num_deg <= deg_bound[j]:
                out_deg = compass[j]
            if deg_bound[-1] < num_deg or num_deg <= deg_bound[0]:
                out_deg = compass[0]
                
    elif type(num_deg) is list:
        out_deg = []
        for i in num_deg:
            for j in range(len(deg_bound)-1):
                if deg_bound[j] < i <= deg_bound[j+1]:
                    out_deg.append(compass[j+1])
            if deg_bound[-1] < i or i <= deg_bound[0]:
                out_deg.append(compass[0])
            
    else:
        raise TypeError('Handling input of type '+type(num_deg).__name__+' is not implemented yet')

    return out_deg

def calc_theta(t,p):
    """
    Calculate potential temperature.

    **Inputs/Outputs**

    Variable   I/O   Description       Units
    --------   ---   -----------       -----
    t           I    Air pressure      Pa
    p           I    Air temperature   K
    theta       O    Potential temp.   K

    """
    return t*(p0/p)**(Rd/cp)

def calc_theta2t(theta,p):
    return theta*(p/p0)**(Rd/cp)

# constants used to calculate moist adiabatic lapse rate
# See formula 3.16 in Rogers&Yau
a = 2./7.
b = eps*L*L/(Rd*cp)
c = a*L/Rd

def calc_fcor(deg_lat):
    """Coriolis parameter for a given latitude in degrees"""
    return 2.*omg*np.sin(np.radians(deg_lat))

def calc_gamma_s(T,p):
    """Calculates moist adiabatic lapse rate for T (Celsius) and p (Pa)
    Note: We calculate dT/dp, not dT/dz
    See formula 3.16 in Rogers&Yau for dT/dz, but this must be combined with
    the dry adiabatic lapse rate (gamma = g/cp) and the 
    inverse of the hydrostatic equation (dz/dp = -RT/pg)"""
    esat = calc_es(T)
    wsat = eps*esat/(p-esat) # Rogers&Yau 2.18
    numer = a*(T+T0) + c*wsat
    denom = p * (1 + b*wsat/((T+T0)**2)) 
    return numer/denom # Rogers&Yau 3.16

def calc_es(T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*np.exp(17.67*T/(T+243.5))

def calc_e(r,p):
    """Returns vapor pressure (Pa) at mixing ratio r (kg/kg) and pressure p (Pa)
    Formula 2.18 in Rogers&Yau"""
    return r*p/(r+eps)

def calc_td(e):
    """Returns dew point temperature (C) at vapor pressure e (Pa)
    Insert Td in 2.17 in Rogers&Yau and solve for Td"""
    return 243.5 * np.log(e/611.2)/(17.67-np.log(e/611.2))

def calc_theta_e(t,p,r):
    """
    Calculate pseudo-equivalent potential temperature.
    Reference: [Bolton, 1980]

    **Inputs/Outputs**

    Variable   I/O   Description       Units
    --------   ---   -----------       -----
    t           I    Air temperature   K
    p           I    Air pressure      Pa
    r           I    Mixing ratio      kg/kg
    thetae      O    Eq. pot. temp.    K

    **Dependencies**
    calc_e, calc_theta
    """
    # Firstly calculate Temp at LCL, note input e must be in Pa
    t_lcl = 2840./(3.5*np.log(t) - np.log(calc_e(r,p)/100.) - 4.805) + 55 # eqn 21, Bolton
    thetae = calc_theta(t,p)*np.exp(((3.376/t_lcl) - 0.00254)*r*1e3*(1+0.81*r)) # eqn 38, Bolton
    return thetae

def calc_theta_w(t,p,r):
    """
    Calculate wet-bulb temperature.
    Reference: [Stull, 2011]

    **Inputs/Outputs**

    Variable   I/O   Description       Units
    --------   ---   -----------       -----
    t           I    Air temperature   K
    p           I    Air pressure      Pa
    r           I    Mixing ratio      kg/kg
    thetaw      O    Wet bulb temp.    K

    **Dependencies**
    calc_e, calc_es
    """
    rh = calc_e(r,p) / calc_es(t)
    
    thetaw = t*np.atan(0.151977*(rh+8.313659)**0.5) + np.atan(t+rh) - np.atan(t-1.676331)
    + 0.00391838*(rh)**1.5 * np.atan(0.023101*rh) - 4.686035
    return thetaw

def spechum2mixr(q):
    """
    Calculate mixing ratio from specific humidity.

    **Inputs/Outputs**

    Variable   I/O   Description       Units
    --------   ---   -----------       -----
    q           I    Specific humidity kg/kg
    r           O    Mixing ratio      kg/kg
    """
    return q/(1-q)

def mixr2spechum(r):
    """
    Calculate specific humidity from mixing ratio.

    **Inputs/Outputs**

    Variable   I/O   Description       Units
    --------   ---   -----------       -----
    r           I    Mixing ratio      kg/kg
    q           O    Specific humidity kg/kg
    """
    return r/(r+1)

def angle_uv(u1, v1, u2, v2):
    """
    Calculate angle between two wind velocity vectors.
    Input is given as zonal and meridional components of wind vector.
    Input can be as numbers or arrays.

    **Inputs/Outputs**

    Variable   I/O   Description               Units
    --------   ---   -----------               -----
    u1          I    1st vector, x-component    m/s
    v1          I    1st vector, y-component    m/s
    u1          I    2nd vector, x-component    m/s
    v1          I    2nd vector, y-component    m/s
    angle       O    Angle between the vectors  rad

    Generic methods:  http://stackoverflow.com/a/13849249
    """
    angle = np.arccos((u1*u2+v1*v2)/np.sqrt((u1**2+v1**2)*(u2**2+v2**2)))
    return angle
