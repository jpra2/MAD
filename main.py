import numpy as np
from pymoab import core, types, rng, topo_util, skinner
import time
import os
import yaml
from preprocess.levels_generator import LevelsGenerator as lv
from preprocess.create_duals import CreateDuals as crd
from preprocess.create_adm_mesh import CreateAdmMesh as cradm
from preprocess.define_wells import DefineWells as dfw

parent_dir = os.path.dirname(__file__)
# parent_parent_dir = os.path.join(parent_dir, 'output')
# parent_parent_dir = os.path.dirname(parent_dir)
input_dir = os.path.join(parent_dir, 'input')
output_dir = os.path.join(parent_dir, 'output')
config_dir = os.path.join(input_dir, 'configs')
pre_dir = os.path.join(parent_dir, 'preprocess')

os.chdir(input_dir)
generator_lvs = lv('mesh.yml', '27x27x27.msh')
generator_lvs.Run()
os.chdir(config_dir)
dfw.wells_mono(generator_lvs.mb, 'wells_mono.yml', generator_lvs.vols_centroids, generator_lvs.all_volumes, generator_lvs.gravity, generator_lvs.Ltot)
adm_mesh = cradm(gen_lv=generator_lvs, dfw=dfw)
adm_mesh.GenerateAdmMesh()
os.chdir(parent_dir)
generator_lvs.mb.delete_entities(generator_lvs.all_faces)
generator_lvs.mb.delete_entities(generator_lvs.all_edges)
print('writting_vtk')
generator_lvs.mb.write_file('out.vtk')
# adm_mesh2 = cradm(mb=generator_lvs.mb, tags_primals=generator_lvs.tags_primals, Ltot=generator_lvs.Ltot, Ls=generator_lvs.Ls, h=generator_lvs.h, vols_centroids=1)



# import pdb; pdb.set_trace()
