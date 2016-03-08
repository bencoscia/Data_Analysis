import numpy as np
import math

# Pore Numbering : Each corner is the location of the center of the pore
#      Pore 1 ------------- Pore 2'
#           /            /
#          /            /
#         /            /
#        /            /
# Pore 4 ------------- Pore 3'

f = open("wiggle_traj.gro", "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

no_atoms = 137  # number of atoms in a single monomer
no_layers = 5  # number of layers in the membrane structure
mon_per_layer = 6  # number of monomers in each layer
no_pores = 4  # number of pores
# no_monomers =  # in case there are random distributions of monomers in each layer
no_ion = no_layers*mon_per_layer*no_pores
no_monomers = no_ion  # number of monomers (= number of ions since there is one monomer per ion in the case of sodium)
traj_points = int(len(a) / ((no_atoms + 1)*(no_monomers)))  # rounds down to the nearest integer. Won't be a problem
# unless there are thousands of trajectory points in which case it'll start rounding up. In that case you can just
# look at the associated .mdp file for the value of traj_points. This is here to be a time saver for the most part

x = []  # list to hold x positions of ions
y = []  # list to hold y positions of ions
z = []  # list to hold z positions of ions. Z axis runs parallel to pore
start = no_atoms*no_layers*mon_per_layer*no_pores + 2  # This is the line number where the sodium ions will be listed

# adds all sodium coordinates to respective x, y and z lists
for k in range(0, traj_points):
    for i in range(k*(start + 1 + no_ion) + start, k*(start + 1 + no_ion) + start + no_ion):
        x.append(float(a[i][20:28]))
        y.append(float(a[i][28:36]))
        z.append(float(a[i][36:44]))

# find the central axis in each pore:
sum_x_traj = np.zeros((no_pores, 1))
sum_y_traj = np.zeros((no_pores, 1))
traj_start = 25  # which trajectory to begin analysis at (may need to wait for system to equilibrate)
for j in range(traj_start, traj_points):
    for k in range(0, no_pores):  # calculates average x and y values in each pore. Taken to be the center of the pore
        sum_x = 0
        sum_y = 0
        for i in range(j*no_ion + (k * no_ion/4), j*no_ion + ((k + 1) * no_ion/4)):
            sum_x += x[i]
            sum_y += y[i]
        sum_x_traj[k] += sum_x
        sum_y_traj[k] += sum_y

# average positions to find the location of the central axis:
x_axis = sum_x_traj * (1/float(((traj_points - traj_start)*no_layers*mon_per_layer)))
y_axis = sum_y_traj * (1/float(((traj_points - traj_start)*no_layers*mon_per_layer)))

# find distance between pores
pore12 = math.sqrt((x_axis[0] - x_axis[1])**2 + (y_axis[0] - y_axis[1])**2)
pore13 = math.sqrt((x_axis[0] - x_axis[2])**2 + (y_axis[0] - y_axis[2])**2)
pore34 = math.sqrt((x_axis[2] - x_axis[3])**2 + (y_axis[2] - y_axis[3])**2)
pore42 = math.sqrt((x_axis[1] - x_axis[3])**2 + (y_axis[1] - y_axis[3])**2)
pore14 = math.sqrt((x_axis[0] - x_axis[3])**2 + (y_axis[0] - y_axis[3])**2)
pore23 = math.sqrt((x_axis[1] - x_axis[2])**2 + (y_axis[1] - y_axis[2])**2)
avg_dist_bw_pores = 0.25*(pore12 + pore23 + pore34 + pore14)

# Find the deviation in distance from central axis
dist_sum = np.zeros((no_pores, 1))  # list to hold sum of distances of ions from center of pore
dist_list = [[], [], [], []]  # need to save each value of distance in order to calculate deviation. 4 lists in a list

for l in range(0, traj_points):
    for k in range(0, no_pores):
        for j in range(0, no_layers):
            for i in range(0, mon_per_layer):  # loops through each pore individually
                dist = math.sqrt((x[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - x_axis[k])**2 +
                                 (y[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - y_axis[k])**2)
                # magnitude of vector pointing from center to point where ion is stationed
                dist_list[k].append(dist)

avg = np.zeros((no_pores, 1))
for i in range(0, no_pores):
    avg[i] = np.mean(dist_list[i])

Average_Distance = np.mean(avg)

std = np.zeros((no_pores, 1))
for i in range(0, no_pores):
    std[i] = np.std(dist_list[i])

Average_Deviation = np.mean(std)
"{0:.2f}".format(13.949999999999999)

print 'The average distance of ions from the pore center are as follows:'
print 'Pore 1: {:2.2f} nm' .format(float(avg[0]))
print 'Pore 2: {:2.2f} nm' .format(float(avg[1]))
print 'Pore 3: {:2.2f} nm' .format(float(avg[2]))
print 'Pore 4: {:2.2f} nm' .format(float(avg[3]))
print ''
print 'The standard deviation in distance of ions from the center of the pore are as follows:'
print 'Pore 1: {:2.2f} nm' .format(float(std[0]))
print 'Pore 2: {:2.2f} nm' .format(float(std[1]))
print 'Pore 3: {:2.2f} nm' .format(float(std[2]))
print 'Pore 4: {:2.2f} nm' .format(float(std[3]))
print ''
print 'The average distance of ions from the pore center is %s nm' %Average_Distance
print 'The average deviation from the pore center is %s nm' %Average_Deviation
print 'So on average, the pore is %s nm wide' %(2*Average_Distance)
print ''
print 'The distance between each pores is as follows:'
print 'Distance from pore 1 to pore 2: {:2.2f} nm' .format(pore12)
print 'Distance from pore 1 to pore 3: {:2.2f} nm' .format(pore13)
print 'Distance from pore 3 to pore 4: {:2.2f} nm' .format(pore34)
print 'Distance from pore 4 to pore 2: {:2.2f} nm' .format(pore42)
print 'Distance from pore 4 to pore 1: {:2.2f} nm' .format(pore14)
print 'Distance from pore 2 to pore 3: {:2.2f} nm' .format(pore23)
print ''
print 'The average pore to pore distance (not including diagonals) is: {:2.2f} nm' .format(avg_dist_bw_pores)
print 'Pore Orientation for reference: '
print '      Pore 1 ------------- Pore 2'
print '           /            /  '
print '          /            /  '
print '         /            /  '
print '        /            /  '
print ' Pore 4 ------------- Pore 3'