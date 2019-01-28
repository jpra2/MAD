import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import time
import os
import yaml

key_lv0 = 'Level0'
key_lvs = 'Levels'
key_ls = 'L'
key_delta = 'delta'
key_gravity = 'gravity'
# parent_dir = os.path.dirname(__file__)
kkeys = [key_lv0, key_lvs, key_ls, key_delta, key_gravity]


class LevelsGenerator:

    def __init__(self, inputname_dict, inputname_h5m):
        global kkeys
        self.__verif = False
        self.kkeys = kkeys
        os.chdir('configs')
        with open(inputname_dict, 'r') as stream:
            self.info = yaml.load(stream)
        os.chdir('..')
        self.Ltot = self.info[self.kkeys[0]]['L']

        self.mb = core.Core()
        self.mb.load_file(inputname_h5m)
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        self.sk = skinner.Skinner(self.mb)
        self.root_set = self.mb.get_root_set()
        self.all_volumes = self.mb.get_entities_by_dimension(self.root_set, 3)
        # e_tags = self.mb.tag_get_tags_on_entity(self.all_volumes[0])
        # print(e_tags)
        # import pdb; pdb.set_trace()
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        self.mtu.construct_aentities(self.all_nodes)
        self.all_faces = self.mb.get_entities_by_dimension(self.root_set, 2)
        self.all_edges = self.mb.get_entities_by_dimension(self.root_set, 1)
        self.vols_centroids = np.array([self.mtu.get_average_position([v]) for v in self.all_volumes])
        self.ns = (len(self.all_volumes))
        self.gravity = self.info[key_gravity]

        self.DefineLenghtNv0()
        self.__verif = True

    def DefineLenghtNv0(self):
        key = self.kkeys[0]
        h = []
        for i in range(3):
            h.append(self.info[key]['L'][i]/self.info[key]['nb'][i])
        self.h = h

    def DefineBoxesLvs(self):

        key_lvs = self.kkeys[1]
        key_ls = self.kkeys[2]
        Ls = []
        n = len(self.info[key_lvs])
        comps = []

        for k in self.info[key_lvs]:
            l = self.info[key_lvs][k][key_ls]
            Ls.append(l[:])

        self.Ls = Ls

        name_tag_primals = 'PRIMALS_NV_'
        name_tag_fine_to_primal = 'FINE_TO_PRIMAL_NV_'
        name_centroid_tag = 'CENTROID_NV_'
        all_names_primals_tags = []
        all_primals_tags = []
        all_fine_to_primal_tags = []
        all_centroids_tag = []


        for i in range(n):

            name_tag = name_tag_primals+str(i+1)
            name_fine_to_primal = name_tag_fine_to_primal+str(i+1)
            name_centroid = name_centroid_tag+str(i+1)
            all_names_primals_tags.append(name_tag)
            primal_tag = self.mb.tag_get_handle(name_tag, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
            fine_to_primal_tag = self.mb.tag_get_handle(name_fine_to_primal, 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
            centroid_tag = self.mb.tag_get_handle(name_centroid, 3, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
            all_centroids_tag.append(centroid_tag)
            all_fine_to_primal_tags.append(fine_to_primal_tag)
            all_primals_tags.append(primal_tag)
            boxes = self.CreateBoxes(Ls[i])
            self.CreateMeshsetsPrimals(boxes, primal_tag, fine_to_primal_tag, centroid_tag)



            ###########################################################
            ### esse passo so eh necessario para malhas nao estruturadas
            ### verifica se existe elementos da malha fina em mais de um volume da malha grossa
            ######## meshsets = self.mb.get_entities_by_type_and_tag(
            ########     self.root_set, types.MBENTITYSET, np.array([primal_tag]),
            ########     np.array([None]))
            ########
            ######## boundary = self.sk.find_geometric_skin(self.all_volumes[0])
            ######## print(boundary)
            ######## import pdb; pdb.set_trace()
            # m0 = set()
            # for m in meshsets:
            #     elems0 = self.mb.get_entities_by_handle(m)
            #     for m2 in meshsets:
            #         if m2 in m0:
            #             continue
            #         elems2 = self.mb.get_entities_by_handle(m2)
            #         inter = rng.intersect(elems0, elems2)
            #
            #         if len(inter) > 0:
            #             for elem in inter:
            #                 ind = np.where(np.array(self.all_volumes) == elem)[0]
            #                 cent_elem = self.vols_centroids[ind]
            #                 cent_m = self.mb.tag_get_data(centroid_tag, m, flat=True)
            #                 cent_m2 = self.mb.tag_get_data(centroid_tag, m2, flat=True)
            #                 dist_m = np.linalg.norm(cent_elem - cent_m)
            #                 dist_m2 = np.linalg.norm(cent_elem - cent_m2)
            #                 if dist_m < dist_m2:
            #                     self.mb.remove_entities(m2, elem)
            #                 else:
            #                     self.mb.remove_entities(m, elem)
            #
            #     m0.add(m)
            #################################################################################


        self.tags_primals = all_primals_tags
        # self.all_names_primals_tags = all_names_primals_tags
        self.all_centroids_tag = all_centroids_tag
        self.all_fine_to_primal_tags = all_fine_to_primal_tags
        all_primals_tags.reverse()# cls.mtu = topo_util.MeshTopoUtil(mb)
        all_names_primals_tags.reverse()

        for i in range(n-1):
            primal_tag = all_primals_tags[i]
            name_primal_tag = all_names_primals_tags[i]
            primal_tag_nv0 = all_primals_tags[i-1]
            name_primal_tag_nv0 = all_names_primals_tags[i-1]

            meshsets = self.mb.get_entities_by_type_and_tag(
                self.root_set, types.MBENTITYSET, np.array([primal_tag]),
                np.array([None]))

            meshsets_nv0 = self.mb.get_entities_by_type_and_tag(
                self.root_set, types.MBENTITYSET, np.array([primal_tag_nv0]),
                np.array([None]))

            ms0 = set()
            for m1 in meshsets:
                cont = 0
                for m0 in meshsets_nv0:
                    if m0 in ms0:
                        continue
                    elems_m1 = self.mb.get_entities_by_handle(m1)
                    elems_m0 = self.mb.get_entities_by_handle(m0)
                    inter = rng.intersect(elems_m1, elems_m0)

                    if len(inter) > 0:
                        if np.allclose(np.array(elems_m0), np.array(inter)):
                            ms0.add(m0)
                            self.mb.add_child_meshset(m1,m0)
                        else:
                            raise RuntimeError('Tamanho das malhas nao compativel')
        #

    def CreateBoxes(self, ll):
        lim =1e-10
        comp = []
        Ltot = self.info[self.kkeys[0]]['L']
        for i in range(3):
            r = int(Ltot[i]//(ll[i]))
            m = Ltot[i]%ll[i]
            points = [k*ll[i] for k in range(r+1)]
            if m > lim:
                points[-1] = Ltot[i]
            comp.append(points[:])

        ids = 0
        boxes = []
        all_points_boxes = []
        centroid_boxes = []
        for i in range(len(comp[0])-1):
            for j in range(len(comp[1])-1):
                for k in range(len(comp[2])-1):
                    points = np.array([
                    np.array([comp[0][i], comp[1][j], comp[2][k]]),
                    np.array([comp[0][i+1], comp[1][j], comp[2][k]]),
                    np.array([comp[0][i], comp[1][j+1], comp[2][k]]),
                    np.array([comp[0][i+1], comp[1][j+1], comp[2][k]]),
                    np.array([comp[0][i], comp[1][j], comp[2][k+1]]),
                    np.array([comp[0][i+1], comp[1][j], comp[2][k+1]]),
                    np.array([comp[0][i], comp[1][j+1], comp[2][k+1]]),
                    np.array([comp[0][i+1], comp[1][j+1], comp[2][k+1]])
                    ])

                    all_points_boxes.append(points)
                    mins = np.array([points[:,0].min(), points[:,1].min(), points[:,2].min()])
                    maxs = np.array([points[:,0].max(), points[:,1].max(), points[:,2].max()])
                    # centroid = np.array([(mins[0]+maxs[0])/2, (mins[1]+maxs[1])/2, (mins[1]+maxs[1])/2])
                    centroid = np.mean(np.array([mins,maxs]), axis=0)
                    boxes.append(np.array([mins,maxs]))
                    centroid_boxes.append(centroid)
                    ids+=1

        all_points_boxes = np.array(all_points_boxes)
        boxes = np.array(boxes)
        centroid_boxes = np.array(centroid_boxes)
        return boxes

    def CreateMeshsetsPrimals(self, boxes, tag, fine_to_primal_tag, centroid_tag):
        id = 0
        global key_delta
        delta = float(self.info[key_delta])

        for box in boxes:
            meshset = self.mb.create_meshset()
            self.mb.tag_set_data(tag, meshset, id)
            indsx = np.where(self.vols_centroids[:,0] > box[0][0]-delta)[0] #min
            indsy = np.where(self.vols_centroids[:,1] > box[0][1]-delta)[0] #min
            indsz = np.where(self.vols_centroids[:,2] > box[0][2]-delta)[0] #min
            c1 = set(indsx) & set(indsy) & set(indsz)
            indsx = np.where(self.vols_centroids[:,0] < box[1][0]+delta)[0] #max
            indsy = np.where(self.vols_centroids[:,1] < box[1][1]+delta)[0] #max
            indsz = np.where(self.vols_centroids[:,2] < box[1][2]+delta)[0] #max
            c2 = set(indsx) & set(indsy) & set(indsz)
            inds = np.array(list(c1 & c2))

            vols = rng.Range(np.array(self.all_volumes)[inds])
            centroid = np.mean(box, axis=0)
            self.mb.tag_set_data(fine_to_primal_tag, vols, np.repeat(id, len(vols)))
            self.mb.add_entities(meshset,vols)
            self.mb.tag_set_data(centroid_tag, meshset, centroid)
            id+=1

    def Run(self):
        self.DefineBoxesLvs()

    def h():
        doc = "The h property."
        def fget(self):
            return self.__h
        def fset(self, value):
            if self.__verif == True:
                raise RuntimeError("You can't to modify the property h")
            self.__h = value
        def fdel(self):
            del self.__h
        return locals()
    h = property(**h())






# inputname_dict = 'ident.yaml'
# inputname = 'test1'
# h5m_file = inputname + '.h5m'
# msh_file = inputname + '.msh'
# vtk_file = inputname + '.vtk'
# generator = LevelsGenerator(inputname_dict, msh_file)
# generator.Run()
# generator.mb.delete_entities(generator.all_faces)
# generator.mb.delete_entities(generator.all_edges)
# generator.mb.write_file(vtk_file)
