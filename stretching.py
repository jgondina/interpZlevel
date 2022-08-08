import numpy as np

class KGrid:
    Rho             = 0  # function at vertical RHO-points
    W               = 1  # function at vertical W-points
    strDesc         = ['at vertical RHO-points', 'at vertical W-points']

class VStretch:
    Original        = 1  # original (Song and Haidvogel, 1994)
    Shchepetkin2005 = 2  # A. Shchepetkin (UCLA-ROMS, 2005)
    Geyer           = 3  # R. Geyer BBL refinement
    Shchepetkin2010 = 4  # A. Shchepetkin (UCLA-ROMS, 2010)

    valid           = [1, 2, 3, 4]  # Valid choices for Vstreching
    strDesc         = ['original (Song and Haidvogel, 1994)', 'A. Shchepetkin (UCLA-ROMS, 2005)',
                       'R. Geyer BBL refinement', 'A. Shchepetkin (UCLA-ROMS, 2010)']


def stretching(Vstretching, theta_s, theta_b, N, kgrid, report = False):
# Compute ROMS vertical coordinate stretching function
#
# Given vertical terrain-following vertical stretching parameters, this
# routine computes the vertical stretching function used in ROMS vertical
# coordinate transformation. Check the following link for details:
#
#    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
#
# On Input:
#
#    Vstretching   Vertical stretching function: [Original, Shchepetkin2005, Geyer, Shchepetkin2010]
#    theta_s       S-coordinate surface control parameter (scalar)
#    theta_b       S-coordinate bottom  control parameter (scalar)
#    N             Number of vertical levels (scalar)
#    kgrid         Depth grid type: Rho, W]
#    report        Flag to report detailed information (OPTIONAL):
#
# On Output:
#
#    s             S-coordinate independent variable, [-1 <= s <= 0] at
#                    vertical RHO- or W-points (vector)
#    C             Nondimensional, monotonic, vertical stretching function,
#                    C(s), 1D array, [-1 <= C(s) <= 0]
#    numLevels     Number of vertical levels.
#
# Adapted from MATLAB code by:
#=========================================================================#
#  Copyright (c) 2002-2012 The ROMS/TOMS Group                            #
#    Licensed under a MIT/X style license                                 #
#    See License_ROMS.txt                           Hernan G. Arango      #
#=========================================================================#


    #--------------------------------------------------------------------------
    # Compute ROMS S-coordinates vertical stretching function
    #--------------------------------------------------------------------------

    if   kgrid == KGrid.W:
        numLevels = N + 1
        levels = np.arange(N+1)
    elif kgrid == KGrid.Rho:
        numLevels = N
        print('HHHHJSDSDSDSJKDSJ', N)
        levels = np.arange(N) + 0.5

    s = (1.0/N)*(levels - N)

    if Vstretching == VStretch.Original:

        if theta_s > 0:
            Ptheta = np.sinh(theta_s*s)/np.sinh(theta_s)
            Rtheta = np.tanh(theta_s*(s + 0.5))/(2.0*np.tanh(0.5*theta_s)) - 0.5
            C = (1.0-theta_b)*Ptheta + theta_b*Rtheta
        else:
            C = s


    elif Vstretching == VStretch.Shchepetkin2005:

        alfa = 1.0
        beta = 1.0

        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s*s))/(np.cosh(theta_s) - 1.0)
            if theta_b > 0:
                Cbot = -1.0 + np.sinh(theta_b*(s + 1.0))/np.sinh(theta_b)
                weigth = ((s + 1.0)**alfa)*(1.0 + (alfa/beta)*(1.0 - (s + 1.0)**beta))
                C = weigth*Csur+(1.0 - weigth)*Cbot
            else:
                C = Csur

        else:
            C = s


    elif Vstretching == VStretch.Geyer:

        if theta_s > 0:
            alpha = 3            # scale factor for all hyperbolic functions
            Cbot =  np.log(np.cosh(alpha*(s + 1)**theta_b))/np.log(np.cosh(alpha)) - 1
            Csur = -np.log(np.cosh(alpha*abs(s) **theta_s))/np.log(np.cosh(alpha))
            weight = (1 - np.tanh(alpha*(s + 0.5)))/2
            C = weight*Cbot + (1 - weight)*Csur
        else:
            C = s


    elif Vstretching == VStretch.Shchepetkin2010:

        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s*s))/(np.cosh(theta_s) - 1.0)
        else:
            Csur = -s**2

        if theta_b > 0:
            Cbot = (np.exp(theta_b*Csur) - 1.0)/(1.0 - np.exp(-theta_b))
            C = Cbot
        else:
            C = Csur


    # Report S-coordinate parameters.
    if report:
        print('')
        print('Vstretching = %i: %s' % (Vstretching, VStretch.strDesc[Vstretching - 1]))
        print('kgrid       = %i: %s' % (kgrid, KGrid.strDesc[kgrid - 1]))
        print('theta_s     = %.3f' % theta_s)
        print('theta_b     = %.3f' % theta_b)
        print('S-coordinate curves: k, s(k), C(k)')

        for k in range(len(s)-1, -1, -1):
            print('  %3g  %20.12e  %20.12e' % (k, s[k], C[k]))


    return s, C, numLevels
