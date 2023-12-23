import numpy as np
import sys
import steps.geom as stetmesh
import steps.rng as srng
import steps.solver as ssolver
import ip3r_model_mesh
from random import randrange as randr

# Set simulation parameters
# sampling timing was similar to experimental data
T_END = 2.5
DT = 0.0036
POINTS = int(T_END / DT)

tpnts = np.arange(0.0, T_END, DT)
ntpnts = tpnts.shape[0]

# Create random number generator
seed = int(sys.argv[1])
r = srng.create('mt19937', 512)
r.initialize(seed)

# Import model
mdl = ip3r_model_mesh.getModel()

# Import geometry
name_mesh_file = str(sys.argv[2])
mesh, central_tets,  to_bleach_tets, cyt_tris, er_tris, cyto_tets = ip3r_model_mesh.gen_geom(name_mesh_file)

# Create solver object
sim = ssolver.Tetexact(mdl, mesh, r)

# set the central node and the region to be bleached as regions of interest
mesh.addROI('central', stetmesh.ELEM_TET, central_tets)
mesh.addROI('tobleach', stetmesh.ELEM_TET, to_bleach_tets)

# Run the simulation and record data
# tables to store nb zsgreen in the region to be bleached and in the central node at each dt
zs_table = []

# Reset the simulation object
sim.reset()

# Set initial concentration of ZSGreen
sim.setCompConc('cyto', 'ZSGreen', 25e-6)

fna = 'out.' + name_mesh_file + '.' + str(seed) 
f = open(fna, "w")

# run the simulation
for i in range(ntpnts):
    # bleaching at t=1000 ms:
    if i == 1000:
	# get basal and current number of ZSGreen
        prebleach_zs_median = np.median(zs_table)
        nb_zs = sim.getROICount('tobleach', 'ZSGreen')

        # bleach 60% of ZSGreen molecules in the region to be bleached
        nb_zs_to_bleach = int(nb_zs - 0.4 * prebleach_zs_median)
        zs_remaining = nb_zs - nb_zs_to_bleach

	# update the nb of ZSGreen and ZSGreen_bleached molecules accordingly
        sim.setROICount('tobleach', 'ZSGreen', zs_remaining)
        sim.setROICount('tobleach', 'ZSGreen_bleached', nb_zs_to_bleach)

    zs_table.append(sim.getROICount('tobleach', 'ZSGreen'))

    # record nb of ZSGreen molecules in central node and in the region to be bleached
    zs_count_node = sim.getROICount('central', 'ZSGreen')
    zs_count_ztb = sim.getROICount('tobleach', 'ZSGreen')

    sim.run(tpnts[i])

    f.write("%d %d\n" % (zs_count_ztb, zs_count_node))
    f.flush()

f.close()

