import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import time
import os
import yaml
import math
# import itertools
# import levels_generator.LevelsGenerator as LevelsGenerator
# from preprocess.levels_generator import LevelsGenerator
import sys
from utils.pymoab_utils import UtilsPymoab as utpy


name_dirichlet_tag = 'DIRICHLET'
name_neumann_tag = 'NEUMANN'
name_tag_lv = 'LEVEL'

def CalculateDistance(dist):
    return math.sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)


class CreateAdmMesh:

    def __init__(self, gen_lv=None, dfw=None, *args, **kwargs):
        global name_tag_lv
        name_tag_lv_0_id = 'l' + str(0) + '_ID'
        self.tags_lvs_id = [gen_lv.mb.tag_get_handle(name_tag_lv_0_id, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)]
        gid_tag = gen_lv.mb.tag_get_handle('GLOBAL_ID')
        gids = gen_lv.mb.tag_get_data(gid_tag, gen_lv.all_volumes, flat=True)
        gids -= gids.min()
        gen_lv.mb.tag_set_data(gid_tag, gen_lv.all_volumes, gids)
        gen_lv.mb.tag_set_data(self.tags_lvs_id[0], gen_lv.all_volumes, gids)
        ## informar quais volumes permanecem num determinado nivel
        self.level_tag = gen_lv.mb.tag_get_handle(name_tag_lv, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
        self.meshsets_with_levels = set()



        self.gen_lv = gen_lv
        self.dfw = dfw

        with open('adm_mesh.yml', 'r') as stream:
            self.info = yaml.load(stream)
        rs = []
        for key in self.info['Rs']:
            rs.append(float(self.info['Rs'][key]))

        self.rs = rs

    def GenerateAdmMesh(self):
        self.gen_lv.tags_primals.reverse()
        # print(self.gen_lv.tags_primals)
        # print(self.gen_lv.all_centroids_tag)
        # print(self.gen_lv.all_fine_to_primal_tags)
        # import pdb; pdb.set_trace()
        for i, tag_primal in enumerate(self.gen_lv.tags_primals):
            name_tag_lv_id = 'l' + str(i+1) + '_ID'
            tag_lv_id = self.gen_lv.mb.tag_get_handle(name_tag_lv_id, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
            self.tags_lvs_id.append(tag_lv_id)
            if i == 0:
                self.SetLv0()
                self.GenerateLv1(tag_primal, tag_lv_id, i+1, self.gen_lv.all_centroids_tag[i])
            else:
                self.GenerateOthersLevels(tag_primal, tag_lv_id, i+1, self.gen_lv.all_centroids_tag[i], self.gen_lv.all_fine_to_primal_tags[i])

    def GenerateLv1(self, tag_primal, tag_lv_id, level, centroid_tag):
        vols_lv_0 = self.gen_lv.mb.get_entities_by_type_and_tag(
        self.gen_lv.mb.get_root_set(), types.MBHEX, np.array([self.level_tag]),
        np.array([level-1]))

        n1 = len(vols_lv_0)

        ms0 = set()
        ms_current = set()

        meshsets_current_level = self.gen_lv.mb.get_entities_by_type_and_tag(
        self.gen_lv.mb.get_root_set(), types.MBENTITYSET, np.array([tag_primal]),
        np.array([None]))
        centroids = self.gen_lv.mb.tag_get_data(centroid_tag, meshsets_current_level)

        ## verificando quais volumes permanecem no nivel anterior
        for meshset in meshsets_current_level:
            elems_in_meshset = self.gen_lv.mb.get_entities_by_handle(meshset)
            inter = rng.intersect(elems_in_meshset, vols_lv_0)
            if len(inter) > 0:
                ms0.add(meshset)
                self.meshsets_with_levels.add(meshset)
                self.gen_lv.mb.tag_set_data(self.level_tag, elems_in_meshset, np.repeat(level-1, len(elems_in_meshset)))
                vols_lv_0 = rng.unite(vols_lv_0, elems_in_meshset)

        self.elems_with_level = rng.Range(vols_lv_0)

        ## verificando quais volumes sao vizinhos dos do nivel anterior
        ms2 = []
        for meshset in ms0:
            vizs = utpy.get_meshsets_viz_vertex(self.gen_lv.mb, meshset, meshsets_current_level)
            for m in vizs:
                if m in ms0:
                    continue

                elems = self.gen_lv.mb.get_entities_by_handle(m)
                self.meshsets_with_levels.add(m)
                self.gen_lv.mb.tag_set_data(self.level_tag, elems, np.repeat(level, len(elems)))
                ms2.append(m)
                self.elems_with_level = rng.unite(self.elems_with_level, elems)

        centroids2 = self.gen_lv.mb.tag_get_data(centroid_tag, ms2)

        ################################################################################
        ## verificando quais volumes que ficam no nivel atual dada uma certa distancia
        ## obs: descomentar
        # for meshset, centroid in zip(ms2, centroids2):
        #     dists = centroids - centroid
        #     dists = np.array(list(map(CalculateDistance, dists)))
        #     dists = dists < self.rs[level]
        #     meshsets = np.array(meshsets_current_level)[dists]
        #     if len(meshsets) > 0:
        #         for m in meshsets:
        #             if m in ms2 or m in ms0:
        #                 continue
        #             self.meshsets_with_levels.add(m)
        #             elems = self.gen_lv.mb.get_entities_by_handle(m)
        #             self.gen_lv.mb.tag_set_data(self.level_tag, elems, np.repeat(level, len(elems)))
        #             self.elems_with_level = rng.unite(self.elems_with_level, elems)
        ##################################################################################
        # n = len(vols_lv_0)
        # self.gen_lv.mb.tag_set_data(tag_lv_id, vols_lv_0, np.arange(0, n))
        ################################################################

        n = 0
        for m in meshsets_current_level:
            elems = self.gen_lv.mb.get_entities_by_handle(m)
            if m in ms0:
                n2 = len(elems)
                self.gen_lv.mb.tag_set_data(tag_lv_id, elems, np.arange(n, n+n2))
                n+=n2
            else:
                self.gen_lv.mb.tag_set_data(tag_lv_id, elems, np.repeat(n, len(elems)))
                n+=1

    def GenerateOthersLevels(self, tag_primal, tag_lv_id, level, centroid_tag, fine_to_primal_tag):
        meshsets_current_level = self.gen_lv.mb.get_entities_by_type_and_tag(
        self.gen_lv.mb.get_root_set(), types.MBENTITYSET, np.array([tag_primal]),
        np.array([None]))
        centroids = self.gen_lv.mb.tag_get_data(centroid_tag, meshsets_current_level)

        ms0 = set() # meshsets locais que permanecem no nivel local caso nao haja mais niveis
        ms02 = set()


        ## verificando quais volumes permanecem no nivel anterior
        for meshset in meshsets_current_level:
            # childs = self.gen_lv.mb.get_child_meshsets(meshset)
            elements_in_meshset = self.gen_lv.mb.get_entities_by_handle(meshset)
            inter1 = rng.intersect(elements_in_meshset, self.elems_with_level)
            if len(inter1) == 0:
                ms0.add(meshset)
                continue
            elif len(inter1) == len(elements_in_meshset):
                ms02.add(meshset)
                continue
            else:
                # len(inter1) > 0 and len(inter1) < len(elements_in_meshset):
                elements_for_set_level = rng.subtract(elements_in_meshset, inter1)
                self.gen_lv.mb.tag_set_data(self.level_tag, elements_for_set_level, np.repeat(level-1, len(elements_for_set_level)))
                self.elems_with_level = rng.unite(self.elems_with_level, elements_for_set_level)
                ms02.add(meshset)

        #################################################################################
        ## verificando quais volumes sao vizinhos dos do nivel anterior
        elements_level_ant = self.gen_lv.mb.get_entities_by_type_and_tag(
        self.gen_lv.mb.get_root_set(), types.MBHEX, np.array([self.level_tag]),
        np.array([level-1]))
        adjs = self.gen_lv.mtu.get_bridge_adjacencies(elements_level_ant, 2, 3)
        adjs = rng.subtract(adjs, self.elems_with_level)

        ms = set()
        ms2 = []
        for adj in adjs:
            if adj in ms:
                continue
            nc = self.gen_lv.mb.tag_get_data(fine_to_primal_tag, adj, flat=True)[0]
            elements = self.gen_lv.mb.get_entities_by_type_and_tag(
            self.gen_lv.mb.get_root_set(), types.MBHEX, np.array([fine_to_primal_tag]),
            np.array([nc]))
            meshset = self.gen_lv.mb.get_entities_by_type_and_tag(
            self.gen_lv.mb.get_root_set(), types.MBENTITYSET, np.array([fine_to_primal_tag]),
            np.array([nc]))
            ms0 = ms0 - set(meshset)
            ms = ms | set(elements)
            ms2.append(meshset)
            ms02.add(meshset)
            self.gen_lv.mb.tag_set_data(self.level_tag, elements, np.repeat(level, len(elements)))
        ############################################################################################

        # centroids2 = self.gen_lv.mb.tag_get_data(centroid_tag, ms2)
        # ###############################################################################
        # # verificando quais volumes que ficam no nivel atual dada uma certa distancia
        # # obs: descomentar
        # for meshset, centroid in zip(ms2, centroids2):
        #     dists = centroids - centroid
        #     dists = np.array(list(map(CalculateDistance, dists)))
        #     dists = dists < self.rs[level]
        #     meshsets = np.array(meshsets_current_level)[dists]
        #     if len(meshsets) > 0:
        #         for m in meshsets:
        #             if m in ms02:
        #                 continue
        #             elems = self.gen_lv.mb.get_entities_by_handle(m)
        #             if set(elems).issubset(set(self.elems_with_level)):
        #                 continue
        #             else:
        #                 self.gen_lv.mb.tag_set_data(self.level_tag, elems, np.repeat(level, len(elems)))
        #                 self.elems_with_level = rng.unite(self.elems_with_level, elems)
        #                 ms02.add(m)
        # #################################################################################
        # n = len(vols_lv_0)
        # self.gen_lv.mb.tag_set_data(tag_lv_id, vols_lv_0, np.arange(0, n))
        ################################################################

        if level == len(self.gen_lv.tags_primals) and len(self.elems_with_level) < len(self.gen_lv.all_volumes):

            for m in ms0:
                elements = self.gen_lv.mb.tag_get_handle(m)
                self.gen_lv.mb.tag_set_data(self.level_tag, elements, np.repeat(level, len(elements)))
                self.elems_with_level = rng.unite(self.elems_with_level, elements)

        ms = set()
        n = 0

        for elem in self.elems_with_level:
            if elem in ms:
                continue
            lev = self.gen_lv.mb.tag_get_data(self.level_tag, elem, flat=True)[0]
            if lev == 0:
                self.gen_lv.mb.tag_set_data(tag_lv_id, elem, n)
                ms.add(elem)
                # n+=1
            else:
                nc = self.gen_lv.mb.tag_get_data(self.gen_lv.all_fine_to_primal_tags[lev-1], elem, flat=True)[0]
                elems = self.gen_lv.mb.get_entities_by_type_and_tag(
                self.gen_lv.mb.get_root_set(), types.MBHEX, np.array([self.gen_lv.all_fine_to_primal_tags[lev-1]]),
                np.array([nc]))
                self.gen_lv.mb.tag_set_data(tag_lv_id, elems, np.repeat(n, len(elems)))
                ms = ms | set(elems)
                # n+=1

            n+=1


    def SetLv0(self):
        wells = self.dfw.get_all_wells()
        cent_wells = np.array([self.gen_lv.mtu.get_average_position([v]) for v in wells])
        vols_lv_0 = rng.Range()
        vols_lv_0 = rng.unite(vols_lv_0, wells)

        for well, cent in zip(wells, cent_wells):
            dists = self.gen_lv.vols_centroids - cent
            dists = np.array(list(map(CalculateDistance, dists)))
            dists = dists < self.rs[0]
            vols = np.array(self.gen_lv.all_volumes)[dists]
            if len(vols) > 0:
                vols_lv_0 = rng.unite(vols_lv_0, rng.Range(vols))

        self.gen_lv.mb.tag_set_data(self.level_tag, vols_lv_0, np.repeat(0, len(vols_lv_0)))

    def Run(self):
        self.GenerateAdmMesh()
