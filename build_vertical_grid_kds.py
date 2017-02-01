# build a vertical grid for ocean models
# this is scripted for MOM which uses the "super grid"
# the "super grid" is a netcdf variable called "zeta"
# which is a vector of the depths at the top and middle of each cell
# this vector is of length 2N-1 where N is the number of z levels in the model
# the functional form of the vertical grid is hyperbolic tangent

# written by Kial Stewart
# this script relates to the Ocean Modelling paper "Vertical resolution of baroclinic modes in global ocean models" by Stewart et al.

# choose some numbers for ocean depth and vertical resolution

# what is the maximum depth of your ocean (approximately)? in meters
H = 6000.0

# what is the maximum grid spacing (the grid spacing at the deepest point in the ocean)? in meters
dzd = 100.0

# what is the minimum grid spacing (the grid spacing at the ocean surface)? in meters
min_dz = 1.0

# how sharp/gentle would you like the hyperbolic tangent (<1 is sharp, 1 is neutral, >1 is gentle)? this is used to "tune" the number of levels to an acceptable amount
depfac = 1.03

# what filename would you like?

grid_filename = "ocean_vertical_grid_for_MOM5.nc"

################
# start the build
################

import netCDF4 as nc
import numpy as np

# define the functional form of the vertical grid
epsilon = 0.001 # this is a small number needed to begin the iteration
def f_all(kk):
	return np.tanh(np.pi*((kk)/(H*depfac)))*(dzd)+epsilon # the function is {tanh(pi*m/H)*dz_max + epsilon}, which is epsilon at the surface and dz_max at H

# make the first two entries of the initial grid; these will be 0 and epsilon for both z and dz
delta_z = [0,epsilon*1.0]
prop_z = [0,epsilon*1.0]

# this is where the magic happens: an iterative process that takes a step from the current end (deepest point) of the grid along the function to find the next point
while prop_z[-1]+delta_z[-1] < 1.2*H:
	aa = np.linspace(1.0,1.5,10000)
	bb = np.zeros([len(aa)])
	loopkill = 1.0
	ii = 0
	while loopkill > 0:
		bb[ii] = (f_all(prop_z[-1]+(delta_z[-1]*aa[ii])))-(delta_z[-1]*aa[ii])
		loopkill = bb[ii]
		ii += 1
	aa_bb = np.polyfit(aa[:ii-1],bb[:ii-1],1)
	dznew = (delta_z[-1]*(np.abs(aa_bb[1]/aa_bb[0])))
	delta_z = np.append(delta_z,dznew)
	prop_z = np.append(prop_z,(prop_z[-1]+delta_z[-1]))

# now that we have an initial grid that follows the desired functional form we need to relocate it vertically so that the grid spacing at the surface is min_dz
new_surf = np.max(np.where(delta_z<min_dz)) # find where the initial grid is min_dz (the surface resolution)
real_prop_z = prop_z[new_surf:]-prop_z[new_surf] # make a new grid that shifts the initial grid vertically
real_delta_z = delta_z[new_surf:] # make a new dz for this new grid
real_prop_z = real_prop_z[np.where(real_prop_z<H)] # cut the new grid off at desired depth, H
real_delta_z = real_delta_z[np.where(real_prop_z<H)] # and the new dz too

################
# for MOM we need the "super grid"
################

super_vgrid = np.zeros([(2*(len(real_prop_z)-1))+1])

for zz in np.arange(len(real_prop_z)):
	super_vgrid[zz*2] = real_prop_z[zz]

for zz in np.arange(len(real_prop_z)-1):
	super_vgrid[1+(zz*2)] = real_prop_z[zz+1] - ((real_prop_z[zz+1] - real_prop_z[zz])/2)

eddyfile = nc.Dataset('./'+grid_filename, 'w', format='NETCDF4')
eddyfile.createDimension('nzv', len(super_vgrid))
zeta = eddyfile.createVariable('zeta','f8',('nzv',))
zeta.units = 'meters'
zeta.standard_name = 'vertical_grid_vertex'
zeta.long_name = 'vgrid'
zeta.author = 'Kial Stewart'
eddyfile.variables['zeta'][:] = super_vgrid
eddyfile.close()

print "SUCCESS!! You now have a vertical grid of ", len(real_prop_z), " levels with grid spacing ranging from ", real_delta_z[0]," to ", real_delta_z[-1]," written to file ",grid_filename

