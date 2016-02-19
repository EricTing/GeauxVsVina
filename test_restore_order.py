#!/usr/bin/env python

import unittest

from restore_order import getPredConfs, getOrderedSdf

VINA_IFN = "./1h6v_NDP_B_1.pdb.pdbqt"
INIT_SDF_IFN = "./1h6v_NDP_B_1.pdb.sdf"
VINA_PRED_PDBQT = "./1h6v_NDP_B_1.pdb.pdbqt.vina.pdbqt"
VINA_PRED_SDF = "./1h6v_NDP_B_1.pdb.pdbqt.vina.sdf"


class TestCentroid(unittest.TestCase):
    def test_case(self):
        init_sdf = [l.rstrip() for l in file(INIT_SDF_IFN)]
        vina_in = [l.rstrip() for l in file(VINA_IFN)]
        all_atom_zones = []
        vina_out = [l.rstrip() for l in file(VINA_PRED_PDBQT)]
        pred_confs = getPredConfs(vina_out)
        all_atom_zones += pred_confs

        ordered_sdfs = []
        for atom_zone in all_atom_zones:
            ordered_sdfs.append(getOrderedSdf(init_sdf, atom_zone, vina_in))

        with open(VINA_PRED_SDF, 'w') as ofs:
            for sdf in ordered_sdfs:
                for line in sdf:
                    ofs.write(line)
                    ofs.write("\n")


if __name__ == "__main__":
    unittest.main(exit=False)
