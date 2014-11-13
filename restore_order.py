from glob import glob
import decimal


def vinaLowestEner(vina_ofn):
    """get the lowest energy in the vina output
    """
    vina_out = [l.rstrip() for l in file(vina_ofn)]
    for l in vina_out:
        if "REMARK VINA RESULT" in l:
            splitted = l.split()
            return float(splitted[3])


def chooseBestVinaResult(vina_out_dir):
    """choose the prediction with the lowest energy among the conformations
    """
    ofns = glob(vina_out_dir + '/*-out.pdbqt')
    best_ofn = ''
    lowest_ener = float('infinity')
    for ofn in ofns:
        my_ener = vinaLowestEner(ofn)
        if my_ener < lowest_ener:
            best_ofn = ofn
    return best_ofn


def pdb_atom_type(one_pdbqt_line):
    """str -> [str]
    
    get the atom characteristics of one line in pdbqt file
    """
    type_zone = one_pdbqt_line[:29] + one_pdbqt_line[55:]
    return type_zone.split()


def pdb_atom_coords(one_pdbqt_line):
    """str -> [str]
    
    get the atom coords of one line in pdbqt file
    """
    coords_zone = one_pdbqt_line[30:55]
    return coords_zone.split()


def toRequiredStr(num_str, precision='0.001'):
    """convert the coord str from format 0.0001 to 0.001
    """
    return str(decimal.Decimal(num_str).quantize(decimal.Decimal(precision)))


def getPredConfs(vina_out):
    """[str] -> [[str]]
    Keyword Arguments:
    vina_out -- get the atom zones for each conf in pdbqt
    """
    pred_confs = []
    predicted_coords = []
    for l in vina_out:
        if 'MODEL' in l:
            predicted_coords = []  # clear
        if 'ATOM' in l:
            predicted_coords.append(l)
        if 'ENDMDL' in l:
            pred_confs.append(predicted_coords)
    
    return pred_confs


def getOrderedSdf(sdf, atom_zone, vina_in):
    init_sdf = list(sdf)        # need a copy since list slicing is used
    native_coords = []
    for l in vina_in:
        if "ATOM" in l:
            native_coords.append(l)

    assert len(atom_zone) == len(native_coords)
    for tup in zip(atom_zone, native_coords):
        assert pdb_atom_type(tup[0]) == pdb_atom_type(tup[1])

    init_sdf_coords = []
    tot_atoms = -1
    marker_idx = -1
    for idx, line in enumerate(init_sdf):
        if "OpenBabel" in line:
            marker_idx = idx
            tot_atoms = int(init_sdf[idx + 2].split()[0])
            init_sdf_coords = init_sdf[idx + 3: idx + 3 + tot_atoms]
            break

    rounded_init_sdf_coords = []
    for line in init_sdf_coords:
        xyz_str = line.split()[:3]
        rounded_xyz_str = [toRequiredStr(float(s)) for s in xyz_str]
        rounded_init_sdf_coords.append(rounded_xyz_str)

    def search4NewCoords(old_xyz_str):
        for idx, line in enumerate(native_coords):
            if old_xyz_str == pdb_atom_coords(line):
                new_xyz_str = pdb_atom_coords(atom_zone[idx])
                return new_xyz_str

    predicted_sdf_coords = [search4NewCoords(c) for
                            c in rounded_init_sdf_coords]
    assert len(predicted_sdf_coords) == len(rounded_init_sdf_coords)

    width = 10
    formated_coords = []
    for idx, old_line in enumerate(init_sdf_coords):
        xyz = predicted_sdf_coords[idx]
        xyz = [toRequiredStr(c, '0.0001') for c in xyz]
        formated_xyz_str = ''.join([c.rjust(width) for c in xyz])
        new_line = formated_xyz_str + old_line[30:]
        formated_coords.append(new_line)

    init_sdf[marker_idx + 3: marker_idx + 3 + tot_atoms] = formated_coords

    return init_sdf


if __name__ == "__main__":
    import sys
    my_id = sys.argv[1]

    vina_ifn = '/work/jaydy/working/InitialSdf/' + my_id + '_init.pdbqt'
    init_sdf_fn = "/work/jaydy/working/InitialSdf/" + my_id + "_init.sdf"
    vina_out_dir = '/home/jaydy/work/working/VinaResult/' + my_id
    vina_out_ofns = glob(vina_out_dir + '/*-out.pdbqt')

    init_sdf = [l.rstrip() for l in file(init_sdf_fn)]
    vina_in = [l.rstrip() for l in file(vina_ifn)]

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
