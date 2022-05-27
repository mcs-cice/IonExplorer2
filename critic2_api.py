from structure import *
from read_write import read_file, write_file
import shlex
import subprocess


def write_input(name, start_p, end_p):

    input = ("CRYSTAL" + ' ' + name + ".cif\n\n"
             + "AUTO SEED WS DEPTH 2\n"
             + "FLUXPRINT\n"
             + "    GRAPH 1\n"
             + "    TEXT\n"
             + "ENDFLUXPRINT\n"
             + "POINT {} {} {}\n".format(*start_p)
             + "POINT {} {} {}".format(*end_p))
    write_file(name + ".incritic", input)
    return None


def run_critic2(path):

    args = shlex.split("critic2 " + path + ".incritic " + path + "_critic2_out.txt")
    p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.communicate()
    return None


class Gradient_path:

    def __init__(self, fract_coords, electron_density):

        self.cp1 = None
        self.cp2 = None
        self.fract_coords = np.array(fract_coords)
        self.cart_coords = None
        self.eldens = np.array(electron_density)
        self.translation = get_path_translation(fract_coords)
        self.dist = None
        self.get_cp = None
        self.get_cart_coord = None
        return None

    def split(self, i):

        h_1 = Gradient_path(self.fract_coords[:i + 1], self.eldens[:i + 1])
        h_1.cp1 = self.cp1
        h_1.get_cart_coord = self.get_cart_coord
        if self.cart_coords is not None:
            h_1.cart_coords = self.cart_coords[:i + 1]
        h_2 = Gradient_path(self.fract_coords[i + 1:], self.eldens[i + 1:])
        h_2.cp2 = self.cp2
        h_2.get_cart_coord = self.get_cart_coord
        if self.cart_coords is not None:
            h_2.cart_coords = self.cart_coords[i + 1:]
        return h_1, h_2

    def calc_cart_coord(self):

        self.cart_coords = [self.get_cart_coord(c) for c in self.fract_coords]
        return None

    def slice(self, ind):

        self.fract_coords = self.fract_coords[ind:]
        self.cart_coords = self.cart_coords[ind:]
        self.eldens = self.eldens[ind:]
        self.translation = get_path_translation(self.fract_coords)
        self.dist = None
        return self

    def calc_dist(self, get_cart_coord):

        self.dist = np.linalg.norm(self.cart_coords[0] - get_cart_coord(self.fract_coords[-1] + self.translation))
        return self.dist

    def check_translation(self, tol=0.2):

        if abs(np.linalg.norm(self.cp1.cart_coord - self.cp2.cart_coord) - self.dist) > tol:
            #print(str(self.cp1) + str(self.cp2) + ' ' + str(self.dist) + " - Wrong translation!")
            raise IOError
        return None

    def flip(self, get_cp):

        self.cp1, self.cp2 = self.cp2, self.cp1
        self.fract_coords = list(reversed(self.fract_coords))
        self.eldens = list(reversed(self.eldens))
        if self.cart_coords is not None:
            self.cart_coords = list(reversed(self.cart_coords))
        # self.translation = - self.translation
        self.retranslate(get_cp)
        return self

    def retranslate(self, get_cp):

        self.cp2 = get_cp((self.cp2.name[0], tuple(self.cp2.site - self.cp1.site + self.translation)))
        self.cp1 = get_cp((self.cp1.name[0], (0, 0, 0)))
        return None

    def copy(self):

        copy = Gradient_path(self.fract_coords, self.translation)
        copy.eldens = np.array(self.eldens)
        copy.translation = np.array(self.translation)
        copy.cart_coords = np.array(self.cart_coords)
        copy.dist = self.dist
        copy.cp1 = self.cp1
        copy.cp2 = self.cp2
        return copy

    def concatenate(self, path):

        self.fract_coords = np.concatenate((self.fract_coords, path.fract_coords), axis=0)
        self.cart_coords = np.concatenate((self.cart_coords, path.cart_coords), axis=0)
        self.eldens = np.concatenate((self.eldens, path.eldens), axis=0)
        self.translation = get_path_translation(self.fract_coords)
        self.dist = None
        return self

    def find_correct_translation(self, get_cp, tol=0.2):

        translations = [(i, j, k) for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)]
        for t in translations:
            cp2 = get_cp((self.cp2.name[0], t))
            if abs(np.linalg.norm(self.cp1.cart_coord - cp2.cart_coord) - self.dist) < tol:
                self.cp2 = cp2
                return True
        return False


def read_flux(file_name):

    data = read_file(file_name)
    gradient_paths = []
    path_coord = []
    electron_density = []
    for line in data:
        if "# End gradient path\n" == line:
            gradient_paths.append(Gradient_path(path_coord, electron_density))
            path_coord, electron_density = [], []
        else:
            s_line = line.split()
            if len(s_line) > 3 and s_line[0] != "#":
                path_coord.append(np.array([float(c) for c in s_line[:3]]))
                electron_density.append(float(s_line[3]))
    return gradient_paths


def get_path_translation(path):

    translation = np.array([0,0,0])
    for i in range(len(path)-1):
        delta = path[i] - path[i+1]
        transl = []
        for j, c in enumerate(delta):
            if c > 0.7 :
                transl += [1]
            elif c < -0.7:
                transl += [-1]
            else:
                transl += [0]
        translation += np.array(transl)
    return translation


def extract_eldens(file_name):

    eldens = []
    data = read_file(file_name)
    # Extracting non equivalent critical point list
    for i, line in enumerate(data):
        if "%% POINT" in line:
            for j in range(i, len(data)):
                if "Field value (f):" in data[j]:
                    eldens += [float(data[j].split(":")[1])]
                    break
    return eldens


def extract_cps(file_name):

    cps = []
    r_non_eq_cp = re.compile(("\s*(\d+)\s+\w+\s+\(-?\d\,-?\d\s*\)\s+(\w+)\s+"
                              + "(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)\s+\w+\s+"
                              + "(-?\d+\.\d+E\+?\-?\d+)\s+(-?\d+\.\d+E\+?\-?\d+)\s+(-?\d+\.\d+E\+?\-?\d+)\s*"))
    data = read_file(file_name)
    # Extracting non equivalent critical point list
    for i, line in enumerate(data):
        if "* Critical point list, final report (non-equivalent cps)\n" == line:
            for j in range(i + 4, len(data)):
                m = r_non_eq_cp.match(data[j])
                if m:
                    cps.append([int(m.group(1)), m.group(2)[0],
                                np.array([float(m.group(3)), float(m.group(4)), float(m.group(5))]),
                                float(m.group(7)), float(m.group(8)), float(m.group(9))])
                else:
                    break
            break
    return cps


def get_equiv_position(point, points, tol=0.1):

    dists = cdist(points, [point])
    min_dist = np.min(dists)
    if min_dist < tol:
        return np.where(dists == min_dist)[0][0]
    return None


def add_cps(structure, cps, atom_numb_1=-1, atom_numb_2=2):

    for cp in cps:
        if cp[1] == 'r':
            structure.add_atom(atom_numb_1, cp[2], properties=cp[3])
        elif cp[1] == 'c':
            structure.add_atom(atom_numb_2, cp[2], properties=cp[3])
    structure.multiply_cell()
    return None


def add_paths(structure, gradien_paths, atom_numb_1=-1, atom_numb_2=2):

    unrecognized_paths = []
    for path in gradien_paths:
        path.cart_coords = [structure.cell.get_cart_coord(fc) for fc in path.fract_coords]
        path.calc_dist(structure.cell.get_cart_coord)
    cps = [a for a in structure.atoms.values() if a.element.number == atom_numb_1 or a.element.number == atom_numb_2]
    cps_coord = [cp.cart_coord for cp in cps]
    paths_points = [path.cart_coords[0] for path in gradien_paths]
    paths_points.extend([path.cart_coords[-1] for path in gradien_paths])
    inds = find_equal_ps(paths_points, cps_coord, tol=0.1)

    for i, path in enumerate(gradien_paths):
        cp1_i = inds[i]
        cp2_i = inds[i + len(gradien_paths)]
        if cp1_i is None:
            raise IOError
        if cp2_i is None:
            path.cp1 = cps[cp1_i]
            unrecognized_paths += [path]
        else:
            path.cp1 = cps[cp1_i]
            path.cp2 = cps[cp2_i]
            path.retranslate(structure.get_atom)
            if path.cp1.element.number == atom_numb_2: ##########
                path.flip(structure.get_atom)
            try:
                path.check_translation()
            except:
                path.find_correct_translation(structure.get_atom)
                path.check_translation()
            structure.add_bond(path.cp1, path.cp2, type='S', weight=path)
    if len(unrecognized_paths) != 0:
        find_lost_cps(structure, unrecognized_paths, atom_numb_2)
    return unrecognized_paths


def find_lost_cps(structure, unrecognized_path, atom_numb_2=2):

    common_paths = []
    recognized = np.zeros(len(unrecognized_path), dtype=np.int8)
    for i, path in enumerate(unrecognized_path):
        if recognized[i] == 0:
            common = [path]
            for j in range(i + 1, len(unrecognized_path)):
                if (np.linalg.norm(path.cart_coords[-1] - unrecognized_path[j].cart_coords[-1]) < 0.1):
                    recognized[i], recognized[j] = 1, 1
                    if path.cp1.name[0] != unrecognized_path[j].cp1.name[0]:
                        common.append(unrecognized_path[j])
        if len(common) > 1:
            common_paths.append(common)
    for paths in common_paths:
        eldens = [p.eldens[-1] for p in paths]
        ind = np.where(eldens == min(eldens))[0][0]
        cp2 = structure.add_atom(2, paths[ind].fract_coords[-1], properties=paths[ind].eldens[-1])
        for path in paths:
            path.cp2 = cp2
            path.retranslate(structure.get_atom)
            try:
                path.check_translation()
            except:
                path.find_correct_translation(structure.get_atom)
                path.check_translation()
            if path.cp1.element.number == atom_numb_2:
                path.flip(structure.get_atom)
            #print('lps', path.cp1.name, path.cp2.name, path.dist)
            structure.add_bond(path.cp1, path.cp2, type='S', weight=path)
    if len(unrecognized_path) != sum(recognized):
        unrecognized_path = np.array(unrecognized_path)
        unrecognized_path = unrecognized_path[np.where(recognized == 0)[0]]
        find_lost_paths(structure, unrecognized_path, atom_numb_2)
    return None

def find_lost_paths(structure, unrecognized_path, atom_numb_2=2):

    recognized = np.zeros(len(unrecognized_path), dtype=np.int8)
    for i, path in enumerate(unrecognized_path):
        for bond in structure.nonequiv_bonds.values():
            if bond.weight is not None:
                ind = get_equiv_position(path.cart_coords[-1], bond.weight.cart_coords, tol=0.02)
                if ind is not None:
                    recognized[i] = 1
                    lost_path = path.copy()
                    lost_path.concatenate(bond.weight.copy().slice(ind))
                    lost_path.cp2 = bond.weight.cp2
                    lost_path.retranslate(structure.get_atom)
                    lost_path.calc_dist(structure.cell.get_cart_coord)
                    lost_path.find_correct_translation(structure.get_atom)
                    if path.cp1.element.number == atom_numb_2:
                        path.flip(structure.get_atom)
                    #print('lps', lost_path.cp1.name, lost_path.cp2.name, lost_path.dist)
                    lost_path.check_translation()
                    structure.add_bond(lost_path.cp1, lost_path.cp2, type='S', weight=path)
                    break
    if len(unrecognized_path) != sum(recognized):
        print('#####LOST CP ARE NOT FOUND!!!#####')
    return None


def write_cps_data(name, structure, atom_numb_1=-1, atom_numb_2=2):

    cps_data = "Electron density\n"
    cps = Structure.sort_atoms(list(structure.nonequiv_atoms.values()))
    for cp in cps:
        if cp.element.number == atom_numb_1 or cp.element.number == atom_numb_2:
            cps_data += cp.element.symbol + str(cp.index) + ' ' + str(cp.properties) + '\n'
    write_file(name + "_cps_data.txt", cps_data)
    return None


def write_path_data(name, structure):

    paths_data = ''
    paths = [b for b in structure.nonequiv_bonds.values() if b.weight is not None]
    paths = Structure.sort_bonds(paths)
    for p in paths:
        paths_data += (p.atom_1.element.symbol + str(p.atom_1.index) + ' '
                       + p.atom_2.element.symbol + str(p.atom_2.index) + ' '
                       + str(p.atom_2.site[0]) + ' '
                       + str(p.atom_2.site[1]) + ' '
                       + str(p.atom_2.site[2]) + '\n')
        for gp in p.weight.fract_coords:
            paths_data += str(gp[0]) + ' ' + str(gp[1]) + ' ' + str(gp[2]) + '\n'
    write_file(name + "g_path_data.txt", paths_data)
    return None


def extract_gradient_map(structure, out_name, flux_name):

    gradient_map = structure.copy()
    add_cps(gradient_map, extract_cps(out_name))
    add_paths(gradient_map, read_flux(flux_name))
    return gradient_map

