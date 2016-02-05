#!/usr/bin/env python

import unittest

from glob import glob
from restore_order import getPredConfs, getOrderedSdf

my_id = '1a07C1'
vina_ifn = '/work/jaydy/working/InitialSdf/' + my_id + '_init.pdbqt'
init_sdf_fn = "/work/jaydy/working/InitialSdf/" + my_id + "_init.sdf"
vina_out_dir = '/home/jaydy/work/working/VinaResult/' + my_id
vina_out_ofns = glob(vina_out_dir + '/*-out.pdbqt')

init_sdf = [l.rstrip() for l in file(init_sdf_fn)]
vina_in = [l.rstrip() for l in file(vina_ifn)]


class TestCentroid(unittest.TestCase):
    def test_case1(self):
        all_atom_zones = []
        for vina_ofn in vina_out_ofns:
            vina_out = [l.rstrip() for l in file(vina_ofn)]
            pred_confs = getPredConfs(vina_out)
            all_atom_zones += pred_confs

    def test_case2(self):
        all_atom_zones = []
        for vina_ofn in vina_out_ofns:
            vina_out = [l.rstrip() for l in file(vina_ofn)]
            pred_confs = getPredConfs(vina_out)
            all_atom_zones += pred_confs

        ordered_sdfs = []
        for atom_zone in all_atom_zones:
            ordered_sdfs.append(getOrderedSdf(init_sdf, atom_zone, vina_in))

        for idx, sdf in enumerate(ordered_sdfs):
            vina_sdf_ofn = vina_out_dir + '/' + my_id + '_' + str(idx) + '.sdf'
            with open(vina_sdf_ofn, 'w') as f:
                for l in sdf:
                    f.write(l + "\n")
            print "writing to", vina_sdf_ofn


if __name__ == "__main__":
    unittest.main(exit=False)
