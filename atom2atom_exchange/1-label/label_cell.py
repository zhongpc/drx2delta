from pymatgen.core import Structure, Species, PeriodicSite
from pymatgen.transformations.advanced_transformations import OrderDisorderedStructureTransformation
from smol.cofe import ClusterSubspace
from smol.cofe.space import get_allowed_species, Vacancy
from copy import deepcopy

import numpy as np
import os

########################################################
# 1. define a simple cluster expansion with smol
########################################################

tet_prim = Structure.from_file('../orderings/tet_prim.cif')
ro =  {'Li':{'Li+':1/5, 'Mn3+':1/5, 'Ti4+':1/5, 'Co3+': 1/5}, 'O':{'O2-': 1/2, 'F-': 1/2}, 'H':{'H+':1/3, 'He+': 1/3}}
tet_prim.replace_species(ro)


prim = Structure.from_file('../orderings/LiF.cif')
ro = {'Li+':{'Li+':1/5, 'Mn3+':1/5, 'Ti4+':1/5, 'Co3+': 1/5}, 'F-':{'O2-': 1/2, 'F-': 1/2}}
prim.replace_species(ro)


subspace = ClusterSubspace.from_cutoffs(prim,
                                      cutoffs={2: 4}, # will include orbits of 2 and 3 sites.
                                      basis='indicator',
                                      supercell_size='volume')

tet_subspace = ClusterSubspace.from_cutoffs(tet_prim,
                                      cutoffs={2: 4}, # will include orbits of 2 and 3 sites.
                                      basis='indicator',
                                      supercell_size='volume')


spinel_super = Structure.from_file('../orderings/LiCoO2_super_refine.cif')
supercell_matrix = subspace.scmatrix_from_structure(structure= spinel_super)

order = OrderDisorderedStructureTransformation(algo = 2)


########################################################
# 2. read in the DRX structures
########################################################

s0_structure_list = []

for file_name in os.listdir('./MC_structures'):
    if file_name.endswith('.cif'):
        s0_structure_list.append(Structure.from_file(f'./MC_structures/{file_name}'))


########################################################
# 3. label the fusion supercell with spinel and drx
########################################################

for s0 in s0_structure_list: #  in np.arange(int(0.2*total_Li), int(0.5*total_Li), 4):

    s_DRX = s0.copy()
    print(s_DRX.composition)

    # s_DRX = order.apply_transformation(s_DRX)

    num_Li = int(s_DRX.composition['Li'])
    
    supercell_matrix = subspace.scmatrix_from_structure(structure= spinel_super)
    tet_occu = tet_subspace.occupancy_from_structure(s_DRX, scmatrix = supercell_matrix)

    spinel_occu = tet_subspace.occupancy_from_structure(spinel_super)
    drx_occu = tet_subspace.occupancy_from_structure(s_DRX)

    tet_cluster_supercell = tet_subspace.expansion_structure.copy()
    tet_cluster_supercell.make_supercell(supercell_matrix)


    fusion_occu = deepcopy(drx_occu)

    Li_specie = Species('Li', oxidation_state= 1)
    Mn_specie = Species('Mn', oxidation_state= 3)
    Co_specie = Species('Co', oxidation_state= 3)
    Ti_specie = Species('Ti', oxidation_state= 4)
    Zr_specie = Species('Zr', oxidation_state= 4)
    H_specie = Species('H', oxidation_state= 1)
    Na_specie = Species('Na', oxidation_state= 1)

    # Li is the Li in the DRX structure, but also on 16c site in the spinel structure
    # Mn is the Mn in the DRX structure, but also on 16d site in the spinel structure

    # Co is the Li in the DRX structure, but on 16d site in the spinel structure
    # Zr is the Mn in the DRX structure, but on 16c site in the spinel structure

    # Ti is the Ti in the DRX structure


    flip_tet_occu = []


    allowed_species = get_allowed_species(tet_cluster_supercell)

    for ii, specie in enumerate(tet_occu):
        if H_specie in allowed_species[ii]:
            flip_tet_occu.append(H_specie)

        elif ('vac' in str(specie)) and (not (H_specie in allowed_species[ii])):
            flip_tet_occu.append(Na_specie)

        elif (spinel_occu[ii] == Li_specie) and (drx_occu[ii] == Li_specie):        
            flip_tet_occu.append(Li_specie)

        elif (spinel_occu[ii] == Co_specie) and (drx_occu[ii] == Mn_specie):        
            flip_tet_occu.append(Mn_specie)

        elif (spinel_occu[ii] == Co_specie) and (drx_occu[ii] == Li_specie):        
            flip_tet_occu.append(Co_specie)

        elif (drx_occu[ii] == Mn_specie) and (spinel_occu[ii] == Li_specie):        
            flip_tet_occu.append(Zr_specie)

        elif (drx_occu[ii] == Ti_specie):        
            flip_tet_occu.append(Ti_specie)

        else:
            flip_tet_occu.append(specie)

    sites = []



    for specie, site in zip(flip_tet_occu, tet_cluster_supercell):
        if not isinstance(specie, Vacancy):  # skip vacancies
            site = PeriodicSite(
                specie, site.frac_coords, tet_cluster_supercell.lattice
            )
            sites.append(site)

    s_fusion = Structure.from_sites(sites)

    s_fusion.to(filename= './drx_fusion_' + s_DRX.composition.reduced_formula + '_Li_' + str(num_Li) + '.cif')
    


