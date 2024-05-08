# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:21:39 2024

@author: pahi9557/
"""

import matplotlib
matplotlib.rcParams['text.usetex'] = True

matplotlib.use("QtAgg") # these is needed for offline graphics saving
# use Agg if you just want to save the plots. 
# use Qt5Agg if you want to see them pop up interactively 
from matplotlib import pyplot as plt, patches
#matplotlib.rcParams.update(matplotlib.rcParamsDefault)
import sys
import os
sys.path.append('/Users/veronicaestrada/Downloads/UnderGrad/Spring_2024/Independent_studies/IS_code/')
from analysator import pytools as pt
import numpy as np
plt.switch_backend("agg")
# same thing here, use agg if you want to save plots, use qt5agg if you want plots to show up in the notebook 
import pandas as pd
#sys.path.append('/Users/veronicaestrada/anaconda3/')



plt.close('all')
plt.rc('font', size=22) #controls default text size
plt.rc('axes', titlesize=26) #fontsize of the title
plt.rc('axes', labelsize=22) #fontsize of the x and y labels
plt.rc('xtick', labelsize=22) #fontsize of the x tick labels
plt.rc('ytick', labelsize=22) #fontsize of the y tick labels
plt.rc('legend', fontsize=22) #fontsize of the legend



#state_0 = "/Users/veronicaestrada/Downloads/UnderGrad/Spring_2024/Independent_studies/IS_code/state00090000.vlsv/"
state_0 = "/Users/veronicaestrada/Downloads/UnderGrad/Spring_2024/Independent_studies/IS_code/state00010000.vlsv/"

# this says "where do we want to look in the grid?"
Rp = 6371e3
xPlane = 0.0*Rp # which yz (x = const.) plane to plot
yPlane = 0.0*Rp # which xz (y = const.) plane to plot
zPlane = 0.0*Rp # which xy (z = const.) plane to plot
colormap_name = "inferno"
#colormap_name = "seismic"
#colormap_name = "inferno"

# this just sets up spatial dimensions 
sm_0 = pt.vlsvfile.VlsvReader(state_0)
[xmin,ymin,zmin,xmax,ymax,zmax] = sm_0.get_spatial_mesh_extent() # simulation box dimensions 
[mx,my,mz] = sm_0.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = sm_0.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
dy = (ymax-ymin)/ny
dx = (zmax-zmin)/nz

# normalization things 
x = np.linspace(xmin, xmax, nx)/Rp
y = np.linspace(ymin, ymax, ny)/Rp
z = np.linspace(zmin, zmax, nz)/Rp

axisLims = np.divide([xmin,xmax,ymin,ymax,zmin,zmax],Rp)

cellids_sorted = sm_0.read_variable("CellID").argsort() # CellID variable should always be saved in VLSV files

# find the number of column for given coordinate value a between amin and amax with na cells (a can be x, y or z)
def chooseNcol(a,amin,amax,na):
    if (a > amax) or (a < amin):
          error("coordinate out of domain")
    da = (amax-amin)/na
    Ncol = int(np.floor((a - amin)/da))
    return Ncol

#these change depending on where you want to look in the cube 
XZ_Ncol = chooseNcol(yPlane,ymin,ymax,ny)
XY_Ncol = chooseNcol(zPlane,zmin,zmax,nz)
YZ_Ncol = chooseNcol(xPlane,xmin,xmax,nx)
# YZ_Ncol = 120
######################matplotlib.use("QtAgg") 
###plotting



# load B Field 
D_i = sm_0.read_variable_info("cellBAverage")
D = np.sqrt( D_i.data[:,0]**2 + D_i.data[:,1]**2 + D_i.data[:,2]**2 )

Dx = D_i.data[:,0][cellids_sorted].reshape(nz, ny, nx)
Dy = D_i.data[:,1][cellids_sorted].reshape(nz, ny, nx)
Dz = D_i.data[:,2][cellids_sorted].reshape(nz, ny, nx)


# reshape variable into a proper 3D data cube
D = D[cellids_sorted].reshape(nz,ny,nx)

# get 2d slices from 3D data cube
D_xz = D[:,XZ_Ncol,:]
D_xy = D[XY_Ncol,:,:]
D_yz = D[:,:,YZ_Ncol]




plt.figure(figsize=(20,10))
#plt.suptitle(fr'B = 0 nT')

       
plt.subplot(1,3,1)

# this is just how i get the quivers to be unit-length. I have log scaled ones WAY below if you want 
U2 = Dx[:,XZ_Ncol,:]/1e-9
W2 = Dz[:,XZ_Ncol,:]/1e-9
U3 = U2 / np.sqrt(U2**2 + W2**2);
W3 = W2 / np.sqrt(U2**2 + W2**2);

# a = plt.imshow(D_xz/1e-9,vmin=0,vmax=10,cmap='viridis',extent=[axisLims[0],axisLims[1],axisLims[4],axisLims[5]],aspect="equal",origin="lower")

a = plt.imshow(D_xz/1e-9,cmap='viridis',vmin = 0,
                       vmax = 10, extent=[axisLims[0],axisLims[1],axisLims[4],axisLims[5]],aspect="equal",origin="lower")

plt.title("xz")
#plt.xlabel(r"x $R_p$", labelpad = 0.002)
plt.xlim(-12,4)
plt.ylim(-8,8)
#plt.ylabel(r"z $R_p$", labelpad = 0.00002)

plt.colorbar(label = 'Magnetic Field Strength, nT', orientation='horizontal')
plt.gca().add_patch(plt.Circle((0, 0), 1, color='blue'))
plt.quiver(x[::10], z[::10], U3[::10,::10], W3[::10,::10], color='white', pivot = 'middle')
# plt.axis("scaled")

plt.subplot(1,3,2)

U4 = Dx[XY_Ncol,:,:]/1e-9
W4 = Dy[XY_Ncol,:,:]/1e-9
U5 = U4 / np.sqrt(U4**2 + W4**2);
W5 = W4 / np.sqrt(U4**2 + W4**2);

a = plt.imshow(D_xy/1e-9,vmin=0,vmax=10,cmap='viridis',extent=[axisLims[0],axisLims[1],axisLims[2],axisLims[3]],aspect="equal",origin="lower")
#Q = plt.streamplot(x, y, U5[:,:].T, W5[:,:].T, color='white')
plt.quiver(x[::10], y[::10], U5[::10,::10], W5[::10,::10], color='white', pivot = 'middle')
plt.title("xy")
plt.xlim(-12,4)
plt.ylim(-8,8)
#plt.xlabel(r"x $R_p$")
#plt.ylabel(r"y $R_p$",labelpad = 0.00002)
plt.colorbar(label = 'Magnetic Field Strength, nT', orientation='horizontal')
# plt.plot(a,b, c = 'white', lw = 4)

plt.gca().add_patch(plt.Circle((0, 0), 1, color='blue'))
# plt.axis("scaled")

  
plt.subplot(1,3,3)

U6 = Dy[:,:,YZ_Ncol]/1e-9
W6 = Dz[:,:,YZ_Ncol]/1e-9
U7 = U6 / np.sqrt(U6**2 + W6**2);
W7 = W6 / np.sqrt(U6**2 + W6**2);

a = plt.imshow(D_yz/1e-9,vmin=0,vmax=10,cmap='viridis',extent=[axisLims[2],axisLims[3],axisLims[4],axisLims[5]],aspect="equal",origin="lower")
#Q = plt.streamplot(x, y, U7[:,:], W7[:,:], color='white')
plt.quiver(y[::10], z[::10], U7[::10,::10], W7[::10,::10], color='white', pivot = 'middle')
plt.title("yz")
plt.xlim(-8,8)
plt.ylim(-8,8)
#plt.xlabel(r"y $R_p$")
#plt.ylabel(r"z $R_p$", labelpad = 0.00002)
plt.gca().add_patch(plt.Circle((0, 0), 1, color='blue'))
plt.colorbar(label = 'Magnetic Field Strength, nT', orientation='horizontal')
# plt.axis("scaled")

plt.tight_layout()
plt.show()
