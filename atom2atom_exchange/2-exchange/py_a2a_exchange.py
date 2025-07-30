from pymatgen.core import Structure, Species, PeriodicSite
import random
import numpy as np
from copy import deepcopy
import argparse
import os

def mkdir(path: str):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
    else:
        print("Folder exists")
    return path


def transform():
    s_label = Structure.from_file(args.input)

    if not ('Co' in str(s_label.composition)):
        return 0

    tet_indices = []
    all_sites= s_label.sites

    for ii, site in enumerate(all_sites):
        if 'H+' in site:
            tet_indices.append(ii)


    oct_indices = []


    for trial in range(1000):

        tet_index = random.choice(tet_indices)
        tet_site = all_sites[tet_index]

        neighbor_list = s_label.get_neighbors(tet_site, r = 1.82)
    #     print(neighbor_list)
    #     print()

        num_Mn = 0
        num_Zr = 0
        num_Co = 0
        num_Li = 0
        num_Ti = 0

        neighbor_indices = []

        for neighbor in neighbor_list:
            if ('Zr' in str(neighbor)):
                num_Zr += 1
                Zr_index = neighbor.index
                neighbor_indices.append(neighbor.index)

            elif ('Co' in str(neighbor)):
                num_Co += 1
                Co_index = neighbor.index
                neighbor_indices.append(neighbor.index)

            elif ('Mn' in str(neighbor)):
                num_Mn += 1
                Mn_index = neighbor.index
                neighbor_indices.append(neighbor.index)

            elif ('Li' in str(neighbor)):
                num_Li += 1
                neighbor_indices.append(neighbor.index)

            elif ('Ti' in str(neighbor)):
                num_Ti += 1
                Ti_index = neighbor.index
                neighbor_indices.append(neighbor.index)

        # print("Mn = {}, Zr = {}, Co = {}, Li = {}".format(num_Mn, num_Zr, num_Co, num_Li))

        current_Mn = s_label.composition['Mn3+']
        current_Zr = s_label.composition['Zr4+']
        current_Ti = s_label.composition['Ti4+']
        prob_Mn = 0.25 * current_Mn / (current_Mn + current_Zr + current_Ti)
        prob_Ti = 1 - current_Ti / (current_Mn + current_Zr + current_Ti)

        flip_Mn = False
        flip_Ti = False
        flip_Zr = False

        if (np.random.rand(1) < prob_Mn)[0]:
            flip_Mn = True
        elif (np.random.rand(1) > prob_Ti)[0]:
            flip_Ti = True
        else:
            flip_Zr = True

        save_next = False
        if (num_Mn == 0) and (num_Zr == 1) and (num_Co >= 1) and (num_Ti == 0) and flip_Zr:
            print(tet_index)
            oct_indices = deepcopy([Zr_index, Co_index])
            valid_tet_index = deepcopy(tet_index)
            valid_neighbor_indices = deepcopy(neighbor_indices)
            save_next = True
            break

        elif (num_Mn == 0) and (num_Zr == 0) and (num_Co >= 1) and (num_Ti == 1) and flip_Ti:
            print(tet_index)
            oct_indices = deepcopy([Ti_index, Co_index])
            valid_tet_index = deepcopy(tet_index)
            valid_neighbor_indices = deepcopy(neighbor_indices)
            save_next = True
            break

        elif (num_Mn == 1) and (num_Zr == 0) and (num_Co >= 1) and (num_Ti == 0) and flip_Mn:
            print(tet_index)
            oct_indices = deepcopy([Mn_index, Co_index])
            valid_tet_index = deepcopy(tet_index)
            valid_neighbor_indices = deepcopy(neighbor_indices)
            save_next = True
            break

        else:
            save_next = False

    if save_next == False:
        return 0

    #########  generate the next protocol structure #########
    s_next = s_label.copy()
    if flip_Zr and save_next:
        next_sites = s_next.sites
        print("flip Zr")
        print(next_sites[oct_indices[0]], next_sites[oct_indices[1]])
        s_next.replace(oct_indices[0], 'Li+')
        s_next.replace(oct_indices[1], 'Mn3+')
        s_next.to(filename = args.input)

    elif flip_Mn and save_next:
        next_sites = s_next.sites
        print("flip Mn")
        print(next_sites[oct_indices[0]], next_sites[oct_indices[1]])
        s_next.replace(oct_indices[0], 'Co3+')
        s_next.replace(oct_indices[1], 'Mn3+')
        s_next.to(filename = args.input)

    elif flip_Ti and save_next:
        next_sites = s_next.sites
        print("flip Ti")
        print(next_sites[oct_indices[0]], next_sites[oct_indices[1]])
        s_next.replace(oct_indices[0], 'Li+')
        s_next.replace(oct_indices[1], 'Ti4+')
        s_next.to(filename = args.input)

    print(s_next.composition)
    print()



    ##### generate the start structure #####
    s_start = s_label.copy()
    start_site = s_start.sites[oct_indices[0]]

    s_start.remove_sites(valid_neighbor_indices) # remove neighbor cations
    sites = s_start.sites
    sites.append(start_site) # add start cation to the end of list

    s_start = Structure.from_sites(sites)

    s_start.replace_species({'Zr4+': {'Mn3+': 1}, 'Co3+': {'Li+':1}})
    s_start.remove_species({'H+'})

    mkdir('./TM_migration_drx/' + args.input +'step_' + args.step +  '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_start/')
    s_start.to(filename= './TM_migration_drx/' + args.input +'step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_start/s.cif')



    ##### generate the intermediate structure #####
    s_inter = s_label.copy()
    if flip_Mn or flip_Zr:
        s_inter.replace(valid_tet_index, 'Mn3+')
    elif flip_Ti:
        s_inter.replace(valid_tet_index, 'Ti4+')

    inter_site = s_inter.sites[valid_tet_index]

    s_inter.remove_sites(valid_neighbor_indices + [valid_tet_index])

    sites = s_inter.sites
    sites.append(inter_site) # add intermediate cation to the end of list

    s_inter = Structure.from_sites(sites)

    s_inter.replace_species({'Zr4+': {'Mn3+': 1}, 'Co3+': {'Li+':1}})
    s_inter.remove_species({'H+'})


    mkdir('./TM_migration_drx/' + args.input +'step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_inter/')
    s_inter.to(filename= './TM_migration_drx/' + args.input +'step_' + args.step +  '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_inter/s.cif' )



    ##### generate the end structure #####
    s_end = s_label.copy()

    if flip_Mn or flip_Zr:
        s_end.replace(oct_indices[1], 'Mn3+')
    elif flip_Ti:
        s_end.replace(oct_indices[1], 'Ti4+')


    end_site = s_end.sites[oct_indices[1]]

    s_end.remove_sites(valid_neighbor_indices + [valid_tet_index])


    sites = s_end.sites
    sites.append(end_site)

    s_end = Structure.from_sites(sites)



    s_end.replace_species({'Zr4+': {'Mn3+': 1}, 'Co3+': {'Li+':1}})
    s_end.remove_species({'H+'})

    mkdir('./TM_migration_drx/' + args.input +'step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_end/')
    s_end.to(filename= './TM_migration_drx/' + args.input +'step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '/INPUT_end/s.cif')


    s_list_1 = s_start.interpolate(end_structure= s_inter,
                                         nimages= 5, pbc= True,
                                         interpolate_lattices = True) # oct to tet path

    s_list_2 = s_inter.interpolate(end_structure= s_end,
                                 nimages= 5, pbc= True,
                                 interpolate_lattices = True) # tet to oct path

    print("#### composition check ####")
    print(s_start.composition)
    print(s_inter.composition)
    print(s_end.composition)

    # save the unrelaxed intermediate states
    # for jj, s in enumerate(s_list_1[1:-1]):
    #     mkdir('./TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_path1_' + str(jj))
    #     s.remove_species({'Na', 'Na+'})
    #     s.to(filename = './TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_path1_' + str(jj) + '/s.cif')
    #
    # for jj, s in enumerate(s_list_2[1:-1]):
    #     mkdir('./TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_path2_' + str(jj))
    #     s.remove_species({'Na', 'Na+'})
    #     s.to(filename = './TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_path2_' + str(jj) + '/s.cif')

    # mkdir('./TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_init')
    # s = s_list_1[0].copy()
    # s.remove_species({'Na', 'Na+'})
    # s.to(filename = './TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_init' + '/s.cif')
    #
    # mkdir('./TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_tet')
    # s = s_list_1[-1].copy()
    # s.remove_species({'Na', 'Na+'})
    # s.to(filename = './TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_tet' + '/s.cif')
    #
    # mkdir('./TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_final' )
    # s = s_list_2[-1].copy()
    # s.remove_species({'Na', 'Na+'})
    # s.to(filename = './TM_migration_drx/INPUT_step_' + args.step + '_16dMn_' + str(s_next.composition['Mn3+']) + '_final' + '/s.cif')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default= "0", help="input structure")

    parser.add_argument("-n", "--step", type=str, default= "0", help="number of step")

    args = parser.parse_args()

    transform()
