import numpy as np
# from math import pi, sqrt
from pymoab import core, types, rng, topo_util, skinner
# import time
# import pyximport; pyximport.install()
# from PyTrilinos import Epetra, AztecOO, EpetraExt  # , Amesos
# import math
# import os
# import shutil
# import random
import sys
# import configparser
import io
import yaml
from pymoab_utils import UtilsPymoab as utpy


class OtherUtils:
    gravity = False
    name_keq_tag = 'K_EQ'
    name_s_grav = 'S_GRAV'
    name_dirichlet_tag = 'P'
    name_neumann_tag = 'Q'
    list_names_tags = [name_keq_tag, name_s_grav, name_dirichlet_tag, name_neumann_tag]

    def __new__(cls, mb=None, comm=None):
        if mb == None:
            raise RuntimeError('Defina o mb core do pymoab')
        elif: comm == None:
            raise RuntimeError('Defina o comm do trilinos')

        cls.list_tags = [mb.tag_get_handle(name) for name in list_names_tags]
        cls.mb = mb
        cls.comm = comm
        # cls.mtu = topo_util.MeshTopoUtil(mb)
        return super().__new__(cls)

    @staticmethod
    def mount_local_problem(map_local, faces_in=None):
        """
        retorna os indices: linha, coluna, valor, sz
        sz = tamanho da matriz
        input:
            map_local: dicionario onde keys = volumes e values = id local
            faces_in: faces internas do conjunto de volumes
        output:
            inds:
        """
        elements = utpy.get_elements(OtherUtils.mb, map_local.keys())
        n = len(elements)

        if faces_in == None:
            # faces = utpy.get_all_faces(OtherUtils.mb, rng.Range(elements))
            # bound_faces = utpy.get_boundary_of_volumes(OtherUtils.mb, elements)
            # faces = rng.subtract(faces, bound_faces)
            faces = rng.subtract(utpy.get_all_faces(OtherUtils.mb, rng.Range(elements)), utpy.get_boundary_of_volumes(OtherUtils.mb, elements))
        else:
            faces = faces_in

        keqs = mb.tag_get_data(OtherUtils.list_tags[0], faces, flat=True)
        elems = [OtherUtils.mb.get_adjacencies(face) for face in faces]
        s_gravs = OtherUtils.mb.tag_get_data(OtherUtils.list_tags[1], faces, flat=True)
        dict_keq = dict(zip(faces, keqs))
        dict_elems = dict(zip(faces, elems))
        dict_s_grav = dict(zip(faces, s_gravs))

        linesM = np.array([])
        colsM = np.array([])
        valuesM = np.array([])
        linesM2 = np.array([])
        valuesM2 = np.array([])
        szM = [n, n]

        b = np.zeros(n)
        s = np.zeros(n)

        for face in faces:
            elems = dict_elems[face]
            keq = dict_keq[face]
            s_grav = dict_s_grav[face]

            linesM = np.append(linesM, [map_local[elems[0]], map_local[elems[1]]])
            colsM = np.append(colsM, [map_local[elems[1]], map_local[elems[0]]])
            valuesM = np.append(valuesM, [-keq, -keq])

            ind0 = np.where(linesM2 == map_local[elems[0]])
            if len(ind0[0]) == 0:
                linesM2 = np.append(linesM2, map_local[elems[0]])
                valuesM2 = np.append(valuesM2, [keq])
            else:
                valuesM2[ind0[0]] += keq

            ind1 = np.where(linesM2 == map_local[elems[1]])
            if len(ind1[0]) == 0:
                linesM2 = np.append(linesM2, map_local[elems[1]])
                valuesM2 = np.append(valuesM2, [keq])
            else:
                valuesM2[ind1[0]] += keq

            s[map_local[elems[0]]] += s_grav
            s[map_local[elems[1]]] -= s_grav

        linesM = np.append(linesM, linesM2)
        colsM = np.append(colsM, linesM2)
        valuesM = np.append(valuesM, valuesM2)

        linesM = linesM.astype(np.int32)
        colsM = colsM.astype(np.int32)

        inds = np.array([linesM, colsM, valuesM, szM, False, False])

        if OtherUtils.gravity == True:
            return inds, s
        else:
            return inds, b

















        # RuntimeError

    @staticmethod
    def get_flow_on_boundary(elements, p_tag, bound_faces=None):
        """
        input:
            elements: volumes ou meshset de volumes
            p_tag = tag da pressao
            bound_faces = faces no contorno do volume
        output:
            dict com keys = elements e values = fluxo na face do contorno
        """
        dict_flux = {}
        elements = utpy.get_elements(elements)
        if bound_faces  == None:
            faces = utpy.get_boundary_of_volumes(OtherUtils.mb, elements)
        else:
            faces = bound_faces

        keqs = OtherUtils.mb.tag_get_data(OtherUtils.list_tags[0], faces, flat=True)
        elems = [OtherUtils.mb.get_adjacencies(face) for face in faces]
        s_gravs = OtherUtils.mb.tag_get_data(OtherUtils.list_tags[1], faces, flat=True)
        dict_keq = dict(zip(faces, keqs))
        dict_elems = dict(zip(faces, elems))
        dict_s_grav = dict(zip(faces, s_gravs))

        for face in faces:
            elems = dict_elems[face]
            keq = dict_keq[face]
            s_grav = dict_s_grav[face]

            p = OtherUtils.mb.tag_get_data(p_tag, elems, flat=True)
            flux = (p[1] - p[0])
            if OtherUtils.gravity == True:
                flux += s_grav

            if elems[0] in elements:
                dict_flux[elems[0]] = flux
            else:
                dict_flux[elems[1]] = -flux

        return dict_flux

    @staticmethod
    def get_slice_by_inds(inds, slice_row, slice_col):
        """
        retorna um slice a partir dos inds
        input:
            inds
            slice_row
            slice_col
        output:
            inds2:
        """

        lines2 = np.array([])
        cols2 = np.array([])
        values2 = np.array([], dtype=np.float64)

        map_l = dict(zip(slice_rows, range(len(slice_rows))))
        map_c = dict(zip(slice_cols, range(len(slice_cols))))
        sz = [len(slice_rows), len(slice_cols)]

        for i in slice_rows:
            assert i in inds[0]
            indices = np.where(inds[0] == i)[0]
            cols = np.array([inds[1][j] for j in indices if inds[1][j] in slice_cols])
            vals = np.array([inds[2][j] for j in indices if inds[1][j] in slice_cols])
            lines = np.repeat(i, len(cols))

            lines2 = np.append(lines2, lines)
            cols2 = np.append(cols2, cols)
            values2 = np.append(values2, vals)

        lines2 = lines2.astype(np.int32)
        cols2 = cols2.astype(np.int32)
        local_inds_l = np.array([map_l[j] for j in lines2]).astype(np.int32)
        local_inds_c = np.array([map_c[j] for j in cols2]).astype(np.int32)

        inds2 = np.array([lines2, cols2, values2, sz, local_inds_l, local_inds_c])
        return inds2
