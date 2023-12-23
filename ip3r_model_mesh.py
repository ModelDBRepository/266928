# coding=utf-8
import steps.model as smod
import steps.geom as stetmesh
import steps.utilities.meshio as smeshio
import steps.rng as rng
import steps.utilities.meshio as meshio
import steps.utilities.meshctrl as meshctrl

# Model of bleaching experiments performed on diffusing ZSGreen molecules
def getModel():
    # Create model container
    mdl = smod.Model()

    # Create ZSGreen molecules
    ZSGreen = smod.Spec('ZSGreen', mdl)
    ZSGreen_bleached = smod.Spec('ZSGreen_bleached', mdl)

    # ER surface sys
    ssys = smod.Surfsys('ssys', mdl)

    # Plasma mb surface
    mb_surf = smod.Surfsys('mb_surf', mdl)

    # Create volume system
    # cyt vol sys
    vsys = smod.Volsys('vsys', mdl)

    # ER vol system
    er_vsys = smod.Volsys('er_vsys', mdl)

    # Diffusion constant of ZSGreen
    DZSGreen = 0.09e-10

    # Create diffusion rule 
    diff_ZSGreen = smod.Diff('diff_ZSGreen', vsys, ZSGreen, DZSGreen)
    diff_ZSGreen_bleached = smod.Diff('diff_ZSGreen_bleached', vsys, ZSGreen_bleached, DZSGreen)

    return mdl


########################################################################
def gen_geom(name_mesh_file):
    # import the tetrahedral mesh
    mesh, nodeproxy, tetproxy, triproxy = meshio.importAbaqus(name_mesh_file, 1e-9)
    # create a compartment comprising all mesh tetrahedrons
    ntets = mesh.countTets()
    tet_groups = tetproxy.blocksToGroups()

    # define the different compartments of the branchlet created with Trelis software https://www.csimsoft.com/trelis
    cyto_tets = tet_groups["EB36"] + tet_groups["EB37"] + tet_groups["EB38"] + tet_groups["EB39"] + tet_groups["EB40"]
    er_tets = tet_groups["EB25"]
    central_tets = tet_groups["EB36"]

    # the id of the bleached region varies depending on the geometry
    if name_mesh_file == '5nodes_ratio1.inp':
        to_bleach_tets = tet_groups["EB36"] + tet_groups["EB39"] + tet_groups["EB40"]

    elif name_mesh_file == '5nodes_ratio2.inp':
        to_bleach_tets = tet_groups["EB36"] + tet_groups["EB37"] + tet_groups["EB40"]
    else:
        to_bleach_tets = tet_groups["EB36"] + tet_groups["EB37"] + tet_groups["EB38"]

    # create cyto compartment
    cyto = stetmesh.TmComp('cyto', mesh, cyto_tets)
    # add volume system to cytosol
    cyto.addVolsys('vsys')

    # Define surfaces
    # ER surface triangles can be defined as the overlap between the cytosolic volume and the ER volume
    ER_TRIS = meshctrl.findOverlapTris(mesh, cyto_tets, er_tets)

    # create the patch for er membrane
    er_patch = stetmesh.TmPatch('er_patch', mesh, ER_TRIS, cyto)
    er_patch.addSurfsys('ssys')

    # Plasma membrane surface triangles
    ASTRO_TRIS = mesh.getSurfTris()

    # create the patch for plasma membrane
    cyto_patch = stetmesh.TmPatch('cyto_patch', mesh, ASTRO_TRIS, icomp=cyto)
    cyto_patch.addSurfsys('mb_surf')

    # return geometry container object
    return mesh, central_tets, to_bleach_tets, ASTRO_TRIS, ER_TRIS, cyto_tets
