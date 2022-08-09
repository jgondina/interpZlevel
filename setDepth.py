import numba
from numba import njit, jit

from misc import msgError
from stretching import *

class VTransform:
    Original = 1   # original transformation
    New      = 2   # new transformation

    strDesc = ('original ROMS', 'ROMS-UCLA')
    valid = [1, 2]  # Valid choices for Vtransform

class IGrid:  # Staggered grid C-type (integer):
    density    = 1  # density points
    streamFunc = 2  # streamfunction points
    uVel       = 3  # u-velocity points
    vVel       = 4  # v-velocity points
    wVel       = 5  # w-velocity points

    strHPointType = ('RHO', 'PSI', 'U', 'V', 'RHO')
    strIGridPointType = ('Density', 'Streamfunction', 'uVel', 'vVel', 'wVel')
    valid = [1, 2, 3, 4, 5]  # Valid choices for igrid


# Transform functions
# @njit([(numba.float64, numba.float64, numba.float64, numba.float64[:,:])])
def oldZ0(s, C, hc, h):
    return (s - C)*hc + C*h

# @njit([(numba.float64[:,:], numba.float64[:,:], numba.float64[:,:])])
def oldZ(z0, zeta, h):
    return z0 + zeta*(1.0 + z0/h)

# @njit([(numba.float64, numba.float64, numba.float64, numba.float64[:,:])])
def newZ0(s, C, hc, h):
    return (hc*s + C*h)/(hc + h)

# @njit([(numba.float64[:,:], numba.float64[:,:], numba.float64[:,:])])
def newZ(z0, zeta, h):
    return zeta + (zeta + h)*z0

def setDepth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta = None, report = False):
# Compute ROMS grid depth from vertical stretched variables
#
# Given a bathymetry (h), free-surface (zeta) and terrain-following
# parameters, this function computes the 3D depths for the requested
# C-grid location. If the free-surface is not provided, a zero value
# is assumed resulting in unperturbed depths.  This function can be
# used when generating initial conditions or climatology data for
# an application. Check the following link for details:
#
#    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
#
# On Input:
#
#    Vtransform    Vertical transformation equation: [Original, New] original transformation: oldZ(), oldZ0(), new transformation:      newZ(), newZ0()
#    Vstretching   Vertical stretching function: [vsOriginal, vaShchepetkin2005, vaGeyer, vaShchepetkin2010]
#    theta_s       S-coordinate surface control parameter (scalar)
#    theta_b       S-coordinate bottom  control parameter (scalar)
#    hc            Width (m) of surface or bottom boundary layer in which
#                    higher vertical resolution is required during
#                    stretching (scalar)
#    N             Number of vertical levels (scalar)
#    igrid         Staggered grid C-type (integer): [density, streamFunc, uVel, vVel, wVel]
#    h             Bottom depth, 2D array at RHO-points (m, positive), h(1:Lp+1,1:Mp+1)
#    zeta          Free-surface, 2D array at RHO-points (m, OPTIONAL), zeta(1:Lp+1,1:Mp+1)
#    report        Flag to report detailed information (OPTIONAL):
#
# On Output:
#
#    z             Depths (m, negative), 3D array
#
# Adapted from MATLAB code by:
#=========================================================================#
#  Copyright (c) 2002-2012 The ROMS/TOMS Group                            #
#    Licensed under a MIT/X style license                                 #
#    See License_ROMS.txt                           Hernan G. Arango      #
#=========================================================================%


    hmin = h.min()

    #--------------------------------------------------------------------------
    #  Check preconditions
    #--------------------------------------------------------------------------

    if Vtransform not in VTransform.valid:
        msgError('Error:  SET_DEPTH - Illegal parameter Vtransform = i' % Vtransform)

    if Vstretching not in VStretch.valid:
        msgError('Error:  SET_DEPTH - Illegal parameter Vstretching = i' % Vstretching)

    if igrid not in IGrid.valid:
        msgError('Error:  SET_DEPTH - Illegal parameter igrid = i' % igrid)

    if hc > hmin and Vtransform == VTransform.Original:

        msgError('Error:  SET_DEPTH - critical depth exceeds minimum\n'
                 ' bathymetry value.\n'
                 '            Vtransform = %i\n'
                 '            hc         = %f\n'
                 '            hmin       = %f\n' % (Vtransform, hc, hmin))

    print('>>>>>', zeta.shape, h.shape)
    if zeta is None:
        zeta = np.zeros(h.shape)


    #--------------------------------------------------------------------------
    # Compute vertical stretching function, C(k):
    #--------------------------------------------------------------------------

    if report:
        print('')
        print('Vtransform   = %i (%s)' % (Vtransform, VTransform.strDesc[Vtransform-1]))
        print('igrid        = %i (at horizontal %s-points)' % (igrid, IGrid.strHPointType[igrid - 1]))


    if igrid == IGrid.wVel:
        kgrid = KGrid.W
    else:
        kgrid = KGrid.Rho


    [s, C, numLevels] = stretching(Vstretching, theta_s, theta_b, N, kgrid, report)

    #--------------------------------------------------------------------------
    #  Bathymetry and free-surface interpolated at requested C-grid node type.
    #--------------------------------------------------------------------------

    if igrid in [IGrid.density, IGrid.wVel]:
            pass

    elif igrid == IGrid.streamFunc:
        h    = 0.25*(h   [:-1,:-1] + h   [1:,:-1] + h   [:-1,1:] + h   [1:,1:])
        zeta = 0.25*(zeta[:-1,:-1] + zeta[1:,:-1] + zeta[:-1,1:] + zeta[1:,1:])

    elif igrid == IGrid.uVel:
        h    = 0.5*(h   [:-1,:] + h   [1:,:])
        zeta = 0.5*(zeta[:-1,:] + zeta[1:,:])

    elif igrid == IGrid.vVel:
        h    = 0.5*(h   [:,:-1] + h   [:,1:])
        zeta = 0.5*(zeta[:,:-1] + zeta[:,1:])


    #--------------------------------------------------------------------------
    # Compute depths (m) at requested C-grid location.
    #--------------------------------------------------------------------------

    z = np.zeros((numLevels,) + h.shape)

    if Vtransform == VTransform.Original:
        Z0 = oldZ0
        Z  = oldZ
    else:  # VTransform.New
        Z0 = newZ0
        Z  = newZ

    for k in range(N):
        if igrid == IGrid.wVel:
            z0 = Z0(s[k], C[k], hc, h)
            if k == 0:
              z[0, :,:] = -h
            else:
              z[k, :,:] = Z(z0, zeta, h)

        else:
            z0 = Z0(s[k], C[k], hc, h)
            aaaa = Z(z0, zeta, h)
            print(z0.shape)
            print(aaaa.shape)
            z[k, :,:] = Z(z0, zeta, h)

    return z
