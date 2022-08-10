# encoding: utf-8
import sys

import numpy as np

import pyroms
import pyroms._interp
import pyroms._remapping
import matplotlib.pyplot as plt

import multiprocessing
queue = multiprocessing.JoinableQueue()



def z22roms(varz, grdz, grd, Cpos='rho', irange=None, jrange=None, \
           spval=1e37, flood=True, dmax=0, cdepth=0, kk=0, \
           mode='linear'):
    """
    var = z2roms(var, grdz, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'	     specify the C-grid position where
				     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e37                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
      - cdepth=0                     critical depth for flooding
                                     if depth<cdepth => no flooding
      - kk
      - mode='linear' or 'spline'    specify the type of interpolation

    Interpolate the variable from z vertical grid grdz to ROMS grid grd
    """

    varz = varz.copy()

    assert len(varz.shape) == 3, 'var must be 3D'

    if mode == 'linear':
        imode = 0
    elif mode == 'spline':
        imode = 1
    else:
        raise  'Warning, %s not supported, defaulting to linear' % mode

    print('2')
    if Cpos == 'rho':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_r[0, :]
        mask = grd.hgrid.mask_rho
    elif Cpos == 'u':
        z = 0.5 * (grdz.vgrid.z[:, :, :-1] + grdz.vgrid.z[:, :, 1:])
        depth = 0.5 * (grd.vgrid.z_r[0, :, :, :-1] + grd.vgrid.z_r[0, :, :, 1:])
        mask = grd.hgrid.mask_u
    elif Cpos == 'v':
        z = 0.5 * (grdz.vgrid.z[:, :-1, :] + grdz.vgrid.z[:, 1:, :])
        depth = 0.5 * (grd.vgrid.z_r[0, :, :-1, :] + grd.vgrid.z_r[0, :, 1:, :])
        mask = grd.hgrid.mask_v
    elif Cpos == 'w':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_w[0, :]
        mask = grd.hgrid.mask_rho
    else:
        raise 'Warning,%s bad position. Use depth at Arakawa-C rho points instead.' % Cpos



    nlev, Mm, Lm = varz.shape
    Nm = depth.shape[0]


    if irange is None:
        irange = (0, Lm)
    else:
        assert varz.shape[2] == irange[1] - irange[0], 'var shape and irange must agreed'

    if jrange is None:
        jrange = (0, Mm)
    else:
        assert varz.shape[1] == jrange[1] - jrange[0], 'var shape and jrange must agreed'

    # flood varz if requested
    if flood is True:
        varz = pyroms.remapping.flood(varz, grdz, Cpos=Cpos,
                                      irange=irange, jrange=jrange, spval=spval,
                                      dmax=dmax, cdepth=cdepth, kk=kk)


    varz = np.concatenate((varz[0:1, :, :], varz, varz[-1:, :, :]), 0)

    a = -9999 * np.ones((1, z.shape[1], z.shape[2]))

    b = z.copy()

    c = 100 * np.ones((1, z.shape[1], z.shape[2]))


    z = np.concatenate((a,b,c), 0)

    var = np.ma.zeros((Nm, Mm, Lm))

    def worker(k, queue, varz, z, depth, mask, imode, spval, irange, jrange):
        """thread worker function"""
        print('    Process %i, started' % k)
        # var[k, :, :] =\
        aaa = pyroms._interp.xhslice(varz,                                               z[:, jrange[0]:jrange[1], irange[0]:irange[1]],
                                              depth[k, jrange[0]:jrange[1], irange[0]:irange[1]],
                                              1 + 0*mask[jrange[0]:jrange[1], irange[0]:irange[1]],
                                              imode, spval)
        print('    Process %i, finished' % k)
        # if k<2:
        #     plt.imshow(aaa)
        #     plt.show()

        queue.put(aaa)
        # queue.put({'idx': k, 'data': aaa})
        print(k)

    print('Creating processes for vertical interpolation')
    jobs = []
    for k in range(Nm):
        p = multiprocessing.Process(target=worker, args=(k, queue, varz, z, depth, mask, imode, spval, irange, jrange))
        jobs.append(p)
        p.start()


    print('  Waiting for processes to finish')


    idx = 0
    while len(jobs)>0:

        # if idx<2:
        jobs[0].join()
        aaa = queue.get()
        # print(jobs[0].__dict__)
        # print('>>>>>>   ', idx, aaa['idx'], aaa['data'])
        print('>>>>>>   ', idx, aaa)
            # plt.imshow(aaa)
            # plt.show()
        # var[k,:,:] = jobs[0].join()
        jobs = jobs[1:]
        idx += 1

    print('there')

    # for k in range(Nm):
    #
    #     var[k, :, :] = pyroms._interp.xhslice(varz,
    #                                           z[:, jrange[0]:jrange[1], irange[0]:irange[1]],
    #                                           depth[k, jrange[0]:jrange[1], irange[0]:irange[1]],
    #                                           mask[jrange[0]:jrange[1], irange[0]:irange[1]],
    #                                           imode, spval)
    #     # mask
    #     var = np.ma.masked_values(var, spval, rtol=1e-5)

    var = np.ma.masked_values(var, spval, rtol=1e-5)

    return var
