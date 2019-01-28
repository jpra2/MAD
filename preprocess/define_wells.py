import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import time
import os
import yaml

key_wells = 'Wells'
key_id = 'id'
key_prescription = 'Prescription'
key_value = 'Value'
key_type_region = 'type_region'
key_delimitation = 'delimitation'
key_gama = 'gama'
name_tag_dirichlet = 'DIRICHLET'
name_tag_Neumann = 'NEUMANN'
kkeys = [key_id, key_prescription, key_value, key_type_region, key_delimitation]
prescs = [name_tag_dirichlet, name_tag_Neumann]


class DefineWells:

    @staticmethod
    def wells_mono(mb, name_yml_file, all_centroids, all_volumes, gravity, Ltot):
        global kkeys
        global key_wells
        global prescs
        global key_gama
        DefineWells.get_bound_tags(mb)
        with open(name_yml_file, 'r') as stream:
            info = yaml.load(stream)

        for well in info[key_wells]:
            props = info[key_wells][well]
            value = props[kkeys[2]]
            if props[kkeys[3]] == 'box':
                box = props[kkeys[4]]
                box = np.array([np.array(box[0]), np.array(box[1])])
                volumes = DefineWells.get_wells_by_box(all_centroids, all_volumes, box)
            if props[kkeys[1]] == 'Dirichlet':
                DefineWells.set_dirichlet(mb, volumes, value, gravity, info[key_gama], Ltot[2])
            elif props[kkeys[1]] == 'Neumann':
                DefineWells.set_neumann(mb, volumes, value)

    @staticmethod
    def get_wells_by_box(all_centroids, all_volumes, box):
        box_volumes = box
        inds0 = np.where(all_centroids[:,0] > box_volumes[0,0])[0]
        inds1 = np.where(all_centroids[:,1] > box_volumes[0,1])[0]
        inds2 = np.where(all_centroids[:,2] > box_volumes[0,2])[0]
        c1 = set(inds0) & set(inds1) & set(inds2)
        inds0 = np.where(all_centroids[:,0] < box_volumes[1,0])[0]
        inds1 = np.where(all_centroids[:,1] < box_volumes[1,1])[0]
        inds2 = np.where(all_centroids[:,2] < box_volumes[1,2])[0]
        c2 = set(inds0) & set(inds1) & set(inds2)
        inds_vols = list(c1 & c2)
        volumes = rng.Range(np.array(all_volumes)[inds_vols])
        return volumes




        # import pdb; pdb.set_trace()

    @staticmethod
    def set_dirichlet(mb, volumes, value, gravity, gama=None, Lz=None):
        mtu = topo_util.MeshTopoUtil(mb)
        cents = np.array([mtu.get_average_position([v]) for v in volumes])
        if gravity == False:
            pressao = np.repeat(value, len(volumes))

        elif gravity == True:
            z_elems_d = -1*cents[:,2]
            delta_z = z_elems_d + Lz
            pressao = gama*(delta_z) + value
        ###############################################
        else:
            print("Defina se existe gravidade (True) ou nao (False)")

        DefineWells.mb.tag_set_data(DefineWells.dirichlet_tag, volumes, pressao)

    @staticmethod
    def set_neumann(mb, volumes, value):
        global name_tag_Neumann
        val = value/len(volumes)
        q = np.repeat(val, len(volumes))
        DefineWells.mb.tag_set_data(DefineWells.neumann_tag, volumes, q)

    @classmethod
    def get_bound_tags(cls, mb):
        global name_tag_Neumann
        global name_tag_dirichlet
        cls.neumann_tag = mb.tag_get_handle(name_tag_Neumann, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        cls.dirichlet_tag = mb.tag_get_handle(name_tag_dirichlet, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        cls.mb = mb

    @staticmethod
    def get_volumes_dirichlet():
        volumes = DefineWells.mb.get_entities_by_type_and_tag(
        DefineWells.mb.get_root_set(), types.MBHEX, np.array([DefineWells.dirichlet_tag]),
        np.array([None]))

        return volumes

    @staticmethod
    def get_volumes_neumann():
        volumes = DefineWells.mb.get_entities_by_type_and_tag(
        DefineWells.mb.get_root_set(), types.MBHEX, np.array([DefineWells.neumann_tag]),
        np.array([None]))

        return volumes

    @staticmethod
    def get_all_wells():
        vols1 = DefineWells.get_volumes_neumann()
        vols2 = DefineWells.get_volumes_dirichlet()
        return rng.unite(vols1, vols2)
