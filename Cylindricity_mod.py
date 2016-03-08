# Calculates Cylindricity including drift in the z coordinate

import numpy as np
import math

f = open("20_Layers_0.5ns.gro", "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

no_atoms = 137  # number of atoms in a single monomer
no_layers = 20  # number of layers in the membrane structure
mon_per_layer = 6  # number of monomers in each layer
no_pores = 4  # number of pores
# no_monomers =  # in case there are random distributions of monomers in each layer
no_ion = no_layers*mon_per_layer*no_pores
no_monomers = no_ion  # number of monomers (= number of ions since there is one monomer per ion in the case of sodium)

x = []  # list to hold x positions of ions
y = []  # list to hold y positions of ions
z = []  # list to hold z positions of ions. Z axis runs parallel to pore
start = no_atoms*no_layers*mon_per_layer*no_pores + 2  # This is the line number where the sodium ions will be listed
for i in range(start, start + no_ion):  # adds all sodium coordinates to respective x, y and z lists
    x.append(float(a[i][20:28]))
    y.append(float(a[i][28:36]))
    z.append(float(a[i][36:44]))

# find the values of z where the plane of each layer should be located (approximately)
f = open("20_Layers_0.5ns.gro", "r")  # .gro file whose positions of Na ions will be read
b = []
for line in f:
    b.append(line)

carb_z = []  # keep track of the z coordinates of the carbonyl carbon. Use this as a basis since it moves less than the
# ion
for i in range(0, no_monomers):
    carb_z.append(float(b[11 + i*no_atoms][36:44]))  # the 11th line contains information about the carbonyl carbon in
    # addition to every 137 lines following it

min = min(carb_z)  # minimum z value of carbonyl carbons
max = max(carb_z)  # maximum z value of carbonyl carbons
avg_dist_between_layers = (max - min)/(no_layers - 1)
z_pts = []  # list to hold values of z - positions which will be fixed along center axes of pores for measurement of
# deviation
for i in range(0, no_layers):
    z_pts.append(i*avg_dist_between_layers)

# find the central axis in each pore:
x_axis = []
y_axis = []
for k in range(0, no_pores):  # calculates average x and y values in each pore. Taken to be the center of the pore
    sum_x = 0
    sum_y = 0
    for i in range(k * no_ion/4, (k + 1) * no_ion/4):
        sum_x += x[i]
        sum_y += y[i]
    x_avg = sum_x / (no_layers*mon_per_layer)
    y_avg = sum_y / (no_layers*mon_per_layer)
    x_axis.append(x_avg)  # average x coordinates [x_coord_pore1, x_coord_pore2, x_coord_pore3, x_coord_pore4]
    y_axis.append(y_avg)  # average y coordinates [y_coord_pore1, y_coord_pore2, y_coord_pore3, y_coord_pore4]

# Find the deviation in distance from central axis
deviations = []  # list to hold the standard deviation in distance from center of pore
avg = []  # hold the average distance of ions from the central axis in each pore
for k in range(0, no_pores):
    for j in range(0, no_layers):
        dist_list = []  # list to hold distances of ions from center of pore. Reset on each loop
        for i in range(0, mon_per_layer):  # loops through each pore individually
            dist = math.sqrt((x[i + mon_per_layer*j + mon_per_layer*no_layers*k] - x_axis[k])**2 +
                             (y[i + mon_per_layer*j + mon_per_layer*no_layers*k] - y_axis[k])**2 +
                             (z[i + mon_per_layer*j + mon_per_layer*no_layers*k] - z_pts[j])**2)
            # magnitude of vector pointing from center to point where ion is stationed
            dist_list.append(dist)  # add the distance from the center to a list
    avg_dist = np.mean(dist_list)
    avg.append(avg_dist)
    deviations.append(np.std(dist_list))  # add standard deviation for pore to list [pore1, pore2, pore3, pore4]

Average_Deviation = np.mean(deviations)
print 'The standard deviation in distance of ions from the center of the pore are as follows:'
print 'Pore 1: %s angstroms' %deviations[0]
print 'Pore 2: %s angstroms' %deviations[1]
print 'Pore 3: %s angstroms' %deviations[2]
print 'Pore 4: %s angstroms' %deviations[3]
print ''
print 'The average deviation from the pore center is %s angstroms' %Average_Deviation
