import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import time
import os
import yaml


class CreateDuals:

    def __new__(cls, mb=None, tags_primals=None, Ltot=None, Ls=None, h=None):
        if mb == None:
            raise RuntimeError('Defina o mb core do pymoab')
        elif tags_primals == None:
            raise RuntimeError('Defina tags_primals ordenada do menor nivel para o maior')
        elif Ltot == None:
            raise RuntimeError('Defina Ltot: comprimento total do reservatorio')
        elif Ls == None:
            raise RuntimeError('Defina Ls: comprimentos dos bolcos dos niveis grossos')
        elif h == None:
            raise RuntimeError('Defina h: comprimentos dos blocos da malha fina')

        cls.tags_primals = tags_primals
        cls.mb = mb
        cls.Ltot = Ltot
        cls.mtu = topo_util.MeshTopoUtil(mb)
        cls.Ls = Ls
        cls.h = h
        return super().__new__(cls)
