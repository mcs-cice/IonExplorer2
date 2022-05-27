
import os
import argparse
from structure import *
import read_write as rw
import critic2_api as cr2
from read_write import CIFReader
from multiprocessing import Pool
from itertools import repeat


class IEArgParser:

    def __init__(self):

        self.input = None
        self.out = None
        self.num_proc = None
        self.max_wind = None
        self.min_dist = None
        self.supercell = None
        self.min_hop_dist = None
        self.max_hop_dist = None
        self.migrating_ions = None
        self.write_paths = None
        self.neighbor_mod = None
        self.parse_args()

    @staticmethod
    def init_parser():

        parser = argparse.ArgumentParser(
            description="IonExplorer (version 2.0) is a tool for calculating ion migration trajectories and " +
                        "corresponding barriers in crystalline ion conductors.")

        parser.add_argument("-i", "--path_to_input", nargs="?", type=str, default=None,
                            help="path_to_input: path to cif file.")

        parser.add_argument("-o", "--output_directory", nargs="?", type=str, default=os.getcwd(),
                            help="output_directory: path to a directory for output files. " +
                                 "Default current working directory.")

        parser.add_argument("-m", "--migrating_ion", nargs="?", type=str, default=None,
                            help="migrating_ion: Atomic symbol of the migrating ions.")

        parser.add_argument("-r", "--distance_range", nargs=2, type=float, default=(0, float("inf")),
                            help="distance_range: The minimal and maximal ion hop distance.")

        parser.add_argument("-s", "--search_for_neighbors", nargs=1, type=int, default=1,
                            help="search_for_neighbors: The method of searching for neighboring sites of migrating " +
                                 "ions for potential hops by criteria of 0 - distance cutoff, " +
                                 "and adjacency of Voronoi polyhedra by 1 - vertices, 2 - edges, and 3 - faces." +
                                 "An example 2, default=1.")

        parser.add_argument("-d", "--min_dist", nargs="?", type=float, default=3.,
                            help="min_dist: The minimal distances between the translation equivalent paths." +
                                 "Default=3.0 angstrom. Valid if super_cell is not specified.")

        parser.add_argument("-p", "--num_proc", nargs="?", type=int, default=2,
                            help="num_proc: The number of processes.")

        parser.add_argument("-w", "--max_wind", nargs="?", type=int, default=10,
                            help="max_wind: The maximal number of windows that a migration path can cross. " +
                                 "The criteria apply in analysis of the gradient paths")

        parser.add_argument("-c", "--super_cell", nargs="?", type=int, default=None,
                            help="super_cell: Supercell coefficients. An example 2 2 2, default=auto.")

        parser.add_argument("-wp", "--write_paths", nargs="?", type=bool, default=True,
                            help="write_paths: Write output for migration paths individually."
                                 "An example False, default=True.")

        return parser

    def parse_args(self):

        args = self.init_parser().parse_args()
        input = args.path_to_input
        if input is None:
            raise ValueError("The path to the input cif file is not specified!")
        if not os.path.isfile(input):
            raise ValueError("The input cif file is not found!")
        self.input = input
        self.out = args.output_directory

        if args.migrating_ion is None:
            raise ValueError("The type of migrating ions are not specified!")
        self.migrating_ions = get_number(args.migrating_ion)
        if self.migrating_ions == 0:
            raise ValueError("The atomic symbol {} is not recognized!".format(args.migrating_ion))
        self.min_hop_dist, self.max_hop_dist = args.distance_range
        self.neighbor_mod = args.search_for_neighbors
        self.min_dist = args.min_dist
        self.num_proc = args.num_proc
        self.max_wind = args.max_wind
        self.write_paths = args.write_paths
        self.supercell = args.super_cell
        if self.supercell is not None and len(self.supercell) != 3:
            raise ValueError("Invalid format of super_cell argument!")
        return None


class IonExplorer:

    def __init__(self):

        self.args = IEArgParser()
        self.structures = CIFReader(self.args.input).read_cif()

    def calc(self):

        profiles = []
        trajectories = []
        migration_maps = []
        for structures in self.enumerate_structures():
            mms = self.calc_potential_mms(structures)
            migration_maps.extend(mms)
            for trajs, profs in self.calc_trajectories(mms):
                profiles.append(profs)
                trajectories.append(trajs)
        count_paths = 0
        for mm in migration_maps:
            out_dir = self.args.out + "/" + mm.structure.name
            n_paths = len(mm.potential_hops)
            mm.field_values = profiles[count_paths: n_paths]
            mm.ion_trajectories = trajectories[count_paths: n_paths]
            count_paths += n_paths
            rw.write_file(out_dir + "/PathData.txt", mm.get_map_info())
            mm.plot_profiles().savefig(out_dir + "/MapProfile.png")
            #rw.write_file(out_dir + "/MigrationMap.cif", str(mm.trajectories_to_cif()))
            for i, paths in enumerate(mm.find_lowest_barrier_maps()):
                if paths is not None:
                    rw.write_file("{}/{}_PeriodicMap.cif".format(out_dir, i + 1), str(mm.trajectories_to_cif(paths)))
                    mm.plot_profiles(paths).savefig("{}/{}_MapProfile.png".format(out_dir, i + 1))

    def enumerate_structures(self):

        count = 0
        structures = []
        for i, structure in enumerate(self.structures):
            if structure is not None:
                structure.name = str(i)
                structures.append(structure)
                count += 1
            if count == self.args.num_proc:
                yield structures
        yield structures

    def calc_trajectories(self, migration_maps):

        args = []
        for mm in migration_maps:
            for i, (a_1, a_2) in enumerate(mm.potential_hops):
                args.append([mm.structure, a_1, a_2, self.args.max_wind, self.args.supercell, self.args.min_dist,
                             self.args.out + "/" + mm.structure.name, str(i), self.args.write_paths])
        return self.run_multiproc(IonExplorer.calc_trajectory, args, self.args.num_proc)

    def calc_potential_mms(self, structures):

        n = len(structures)
        args = zip(structures, repeat(self.args.migrating_ions, n), repeat(self.args.neighbor_mod, n),
                   repeat(self.args.min_hop_dist, n), repeat(self.args.max_hop_dist, n))
        return self.run_multiproc(IonExplorer.calc_potential_mm, args)

    @staticmethod
    def calc_potential_mm(structure, ion, mod, min_hop_dist, max_hop_dist):

        mm = MigrationMap(structure)
        mm.find_potential_hops(ion, mod, min_hop_dist, max_hop_dist)
        return mm

    @staticmethod
    def calc_trajectory(structure, atom_1, atom_2, max_wind,
                        supercell, min_dist, out_dir, name, write_paths):

        pf = PathFinder(structure, atom_1, atom_2)
        return pf.calc_trajectory(max_wind, supercell, min_dist, out_dir, name, write_paths)

    @staticmethod
    def run_multiproc(func, args, n_proc=1):

        pool = Pool(n_proc)
        res = pool.starmap(func, args)
        pool.close()
        pool.join()
        return res


class PathFinder:

    def __init__(self, structure, atom_1, atom_2):

        self.structure = structure.copy()
        self.structure.remove_connectivity()
        self.structure.atoms = {k: v for k, v in self.structure.equiv_atoms.items()}
        self.atom_1 = self.structure.get_atom((atom_1.name[0], (0, 0, 0)))
        self.atom_2 = self.structure.get_atom((atom_2.name[0], tuple(atom_2.site - atom_1.site)))
        self.field_1 = None
        self.field_2 = None
        self.transformation = [[0, 0, 0], [1, 1, 1]]
        self.gradient_map = None
        self.optimal_path = []
        self.trajectory = []
        self.field_values = []

    def calc_trajectory(self, max_wind, supercell_size, eq_hop_min_dist, out_dir, name, write_paths=True):

        self.transform_to_supercell(supercell_size, eq_hop_min_dist)
        self.calc_gradient_map(out_dir, name)
        #rw.write_file(out_dir + "/" + str(name) + "_cps.txt", self.get_cps_data())
        #rw.write_file(out_dir + "/" + str(name) + "_gradient_map.cif", str(self.gradient_map))
        self.find_path(max_wind)
        if write_paths:
            rw.write_file(out_dir + "/" + str(name) + "_path.cif", str(self.structure))
            rw.write_file(out_dir + "/" + str(name) + "_profile.txt", self.get_profile())
            self.write_vasp_inp(out_dir, name)
        return self.trajectory, self.field_values

    def transform_to_supercell(self, supercell_size=None, eq_hop_min_dist=3.0, center_path=True):

        supercell = self.structure
        atom_1, atom_2 = self.atom_1, self.atom_2
        cell = [supercell.cell.a, supercell.cell.b, supercell.cell.c]
        hop_vector = atom_1.fract_coord - atom_2.fract_coord - atom_2.site
        if supercell_size is None:
            supercell_size = np.int_(np.ceil(np.abs(hop_vector) + np.full((3, 1), eq_hop_min_dist)[0] / cell))
        elif any(np.abs(hop_vector) >= supercell_size):
            raise Exception(
                "The ion hop is out of the unit cell! Use --super_cell auto or increase supercell size.")
        supercell.transform_to_supercell(*supercell_size)
        supercell.remove_symmetry()
        shift = [0, 0, 0]
        if center_path:
            shift = [0.5, 0.5, 0.5] - (atom_1.fract_coord + atom_2.fract_coord + atom_2.site) / 2.
            self.atom_2 = supercell.get_atom((atom_2.name[0], (0, 0, 0)))
            for a in supercell.nonequiv_atoms.values():
                a.fract_coord = Symmetry.calc_coord_and_translation(a.fract_coord + shift)[0]
                a.cart_coord = supercell.cell.get_cart_coord(a.fract_coord)
        self.transformation = [shift, supercell_size]
        return self.structure, self.transformation

    def calc_gradient_map(self, out_dir, name):

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        path = out_dir + "/" + name
        structure = self.structure.copy()
        structure.remove_atoms([structure.get_atom(self.atom_1.name), structure.get_atom(self.atom_2.name)])
        rw.write_file(path + ".cif", str(structure))
        cr2.write_input(path, self.atom_1.fract_coord, self.atom_2.fract_coord)
        cr2.run_critic2(path)
        gradient_map = cr2.extract_gradient_map(structure, path + "_critic2_out.txt", path + "_flux.txt")
        for b in gradient_map.equiv_bonds.values():
            if b.weight.cp2 is None:
                b.weight.cp2 = b.atom_2
            for c in b.weight.fract_coords:
                b.weight.cart_coords.append(gradient_map.cell.get_cart_coord(c))
        self.gradient_map = gradient_map
        self.field_1, self.field_2 = cr2.extract_eldens(path + "_critic2_out.txt")
        os.remove(path + ".cif")
        os.remove(path + ".incritic")
        os.remove(path + "_flux.txt")
        os.remove(path + "_critic2_out.txt")
        return None

    def find_path(self, length_limit=100, tol=1e-5):

        map = self.gradient_map
        atom_1, atom_2 = self.find_path_ends()
        graph, a_dict, b_dict = map.to_graph()
        weights = [a.properties for a in graph.vertex_attributes["atoms"]]
        path = graph.find_paths((a_dict[atom_1.name[0]], (0, 0, 0)),
                                (a_dict[atom_2.name[0]], (0, 0, 0)), weights, length_limit, tol)[0]
        cps = graph.vertex_attributes["atoms"]
        path = [map.get_atom((cps[i].name[0], tuple(s))) for i, s in path]
        self.optimal_path = path
        self.field_values = [self.field_1, path[0].properties]
        trajectory = [self.atom_1.fract_coord, path[0].fract_coord]
        for i, a_1 in enumerate(path[:-1]):
            bond = map.get_bond((a_1.name, path[i + 1].name, "S"))
            gradient = bond.weight
            if a_1.name == bond.atom_1.name:
                tr = a_1.site - gradient.cp1.name[1]
                trajectory.extend([*gradient.fract_coords + tr, path[i + 1].fract_coord + path[i + 1].site])
                self.field_values.extend([*gradient.eldens, path[i + 1].properties])
            else:
                tr = a_1.site - gradient.cp2.name[1]
                trajectory.extend([*gradient.fract_coords[::-1] + tr, path[i + 1].fract_coord + path[i + 1].site])
                self.field_values.extend([*gradient.eldens[::-1], path[i + 1].properties])
        trajectory.append(self.atom_2.fract_coord)
        self.field_values.append(self.field_2)
        self.field_values = np.array(self.field_values)
        for c in trajectory:
            self.structure.add_atom(0, c, tol=5e-3)
        self.trajectory = (np.array(trajectory) - self.transformation[0]) * self.transformation[1]
        self.remove_duplicate_points()
        return self.trajectory, self.field_values

    def remove_duplicate_points(self, tol=1e-6):

        inds = [0]
        unique_points = [self.trajectory[0]]
        for i, p in enumerate(self.trajectory):
            for up in unique_points:
                if np.linalg.norm(p - up) < tol:
                    break
            else:
                inds.append(i)
                unique_points.append(p)
        self.trajectory = self.trajectory[inds]
        self.field_values = self.field_values[inds]
        return None

    def find_path_ends(self):

        ends = []
        map = self.gradient_map
        s_1 = map.cell.get_cart_coord(self.atom_1.fract_coord)
        s_2 = map.cell.get_cart_coord(self.atom_2.fract_coord)
        paths = list(map.equiv_bonds.values())
        pts = [c for p in paths for c in [p.atom_1.cart_coord, p.atom_2.cart_coord, *p.weight.cart_coords]]
        pt_inds = [i for p in paths for i in [-1, -2, *range(len(p.weight.cart_coords))]]
        gpath_inds = [p for p in paths for _ in range(len(p.weight.cart_coords) + 2)]
        for ind in (find_closest_points(s_1, pts, tol=1e-8)[0], find_closest_points(s_2, pts, tol=1e-8)[0]):
            gp_ind = pt_inds[ind]
            gpath = gpath_inds[ind]
            if gp_ind == -1:
                ends.append(gpath.atom_1)
            elif gp_ind == -2:
                ends.append(gpath.atom_2)
            else:
                atom = map.add_atom(0, gpath.weight.fract_coords[gp_ind], properties=gpath.weight.eldens[gp_ind])
                g_1, g_2 = gpath.weight.split(gp_ind)
                g_1.cp2, g_2.cp1 = atom, atom
                map.add_bond(gpath.atom_1, map.get_atom((atom.name[0], tuple(g_1.translation))), weight=g_1)
                map.add_bond(atom, map.get_atom((gpath.atom_2.name[0], tuple(g_2.translation))), weight=g_2)
                map.remove_bonds([gpath])
                ends.append(atom)
        return ends

    def get_cps_data(self):

        cps_data = "{}\t{}\t{}\t{}\t{}\n".format("CPs", 'X', "Y", "Z", "FIELD_VALUE")
        for a in self.gradient_map.nonequiv_atoms.values():
            if a.properties is not None:
                cps_data += ("{}\t".format(a.element.symbol + str(a.index)) +
                             "{:11.10f}\t{:11.10f}\t{:11.10f}\t{:11.10f}\n".format(*a.fract_coord, a.properties))
        cps_data += "\n"
        for b in self.gradient_map.nonequiv_bonds.values():
            if b.type == "S":
                cps_data += (b.atom_1.element.symbol + str(b.atom_1.index) + "\t" +
                             b.atom_2.element.symbol + str(b.atom_2.index) + "\n")
                for i, c in enumerate(b.weight.fract_coords):
                    cps_data += ("{}\t".format("X" + str(i + 1)) +
                                 "{:11.10f}\t{:11.10f}\t{:11.10f}\t{:11.10f}\n".format(*c, b.weight.eldens[i]))
        return cps_data

    def get_profile(self):

        data = "{}\t{}\t{}\t{}\n".format("X", "Y", "Z", "FIELD_VALUE")
        for i, p in enumerate(self.trajectory):
            data += "{:11.10f}\t{:11.10f}\t{:11.10f}\t{:11.10f}\n".format(*p, self.field_values[i])
        return data

    def write_vasp_inp(self, out_dir, name, image_dist=0.2):

        el_number = self.atom_1.element.number
        structure = self.structure.copy()
        structure.remove_atoms([a for a in structure.nonequiv_atoms.values() if a.element.number == 0])
        structure.remove_atoms([self.atom_1, self.atom_2])
        out_dir = out_dir + "/vasp"
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        out_dir = out_dir + "/" + name
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        points = [self.optimal_path[0]]
        for p in self.optimal_path:
            if np.linalg.norm(points[-1].cart_coord - p.cart_coord) > image_dist:
                points.append(p)
        for i, p in enumerate(points):
            path = out_dir + "/" + str(i)
            if not os.path.isdir(path):
                os.mkdir(path)
            img = structure.copy()
            img.add_atom(el_number, p.fract_coord)
            rw.write_file(path + "/" + "POSCAR", img.to_poscar())
        for i, p in enumerate(points):
            structure.add_atom(el_number, p.fract_coord)
        rw.write_file(out_dir + "/" + "POSCAR", structure.to_poscar())


class MigrationMap:

    def __init__(self, structure):

        self.structure = structure
        self.migration_map = None
        self.potential_hops = None
        self.migration_paths = None
        self.field_values = []
        self.ion_trajectories = []
        self.lb_maps = [None, None, None]

    def get_map_info(self, paths_ind=None):

        if paths_ind is None:
            paths_ind = range(len(self.field_values))
        data = "N\tX\tY\tZ\tFIELD_VALUE\n"
        for i in paths_ind:
            field = self.field_values[i]
            for j, p in enumerate(self.ion_trajectories[i]):
                data += "{}\t{:11.10f}\t{:11.10f}\t{:11.10f}\t{:11.10f}\n".format(i, *p, field[j])
        return data

    def trajectories_to_cif(self, paths_ind=None, tol=-1):

        migration_map = self.structure.copy()
        if paths_ind is None:
            trajectories = np.array(self.ion_trajectories, dtype=object)
        else:
            trajectories = np.array(self.ion_trajectories, dtype=object)[paths_ind]
        for trajectory in trajectories:
            for c in trajectory:
                migration_map.add_atom(0, c, tol)
        return str(migration_map)

    def find_lowest_barrier_maps(self, tol=1e-6):

        hop_inds = {}
        barriers = [None] * len(self.migration_map.equiv_bonds)
        mm_graph, atoms_dict, edge_inds = self.migration_map.to_graph()
        mm_graph.edges_attributes["hop_inds"] = [None] * len(self.migration_map.equiv_bonds)
        for i, (a_1, a_2) in enumerate(self.potential_hops):
            hop_inds[self.migration_map.get_bond((a_1.name, a_2.name, "S")).equal.name] = i
        for bond in self.migration_map.equiv_bonds.values():
            hop_ind = hop_inds[bond.equal.name]
            edge_ind = edge_inds[bond.name]
            field_profile = self.field_values[hop_ind]
            barriers[edge_ind] = max(field_profile) - min(field_profile)
            mm_graph.edges_attributes["hop_inds"][edge_ind] = hop_ind
        mm_graph.edges_attributes["hop_inds"] = np.array(mm_graph.edges_attributes["hop_inds"])
        subgraphs = mm_graph.find_weightless_periodic_subgraphs(barriers, tol)
        self.lb_maps = [list(set(sgr.edges_attributes["hop_inds"])) if sgr is not None else None for sgr in subgraphs]
        return self.lb_maps

    def find_potential_hops(self, migrating_ion, mod=1, min_dist=0, max_dist=float('inf')):

        migration_map = self.structure.copy()
        migration_map.remove_connectivity()
        to_remove = [a for a in migration_map.nonequiv_atoms.values() if a.element.number != migrating_ion]
        migration_map.remove_atoms(to_remove)
        migration_map.multiply_cell()
        atoms = list(migration_map.atoms.values())
        points = np.array([a.cart_coord for a in atoms])
        central_points = [i for i, a in enumerate(atoms) if a.symmetry == 1 and (a.translation == [0, 0, 0]).all()]
        if mod == 0:
            neighbor_sites = find_neighbor_ps(points[central_points], points, min_dist, max_dist)
            neighbors = [(atoms[central_points[i]], atoms[j]) for i, inds in enumerate(neighbor_sites) for j in inds]
        else:
            neighbor_sites = find_vneighbors(points, central_points, mod, min_dist, max_dist).items()
            neighbors = [(atoms[i], atoms[j]) for i, inds in neighbor_sites for j in inds]
        for a1, a2 in neighbors:
            migration_map.add_bond(a1, a2)
        non_equiv_neighbors = [(b.atom_1, b.atom_2) for b in migration_map.nonequiv_bonds.values()]
        non_equiv_neighbors.sort(key=lambda x: np.linalg.norm(x[0].cart_coord - x[1].cart_coord))
        self.migration_map = migration_map
        self.potential_hops = np.array(non_equiv_neighbors)
        return self.potential_hops

    def plot_profiles(self, paths_list=None, title=""):

        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker

        if paths_list is None:
            paths_list = np.arange(len(self.field_values))
        values = np.array(self.field_values, dtype=object)[paths_list]
        fig, ax = plt.subplots(1, figsize=(8, 7.2))
        matplotlib.rcParams.update({'font.size': 13})
        fig.suptitle(r"$\mathrm{" + title + "}$", fontsize=25)
        plt.xlabel("Path (%)", fontsize=20)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(list(range(0, 101, 10)))
        ax.set_ylabel("Electron density (a.u) 10^3", fontsize=20)
        for i, vs in enumerate(values):
            vs = (np.array(vs) - min(vs)) * 1000
            xs = np.array([i * 100 / (len(vs) - 1) for i in range(len(vs))])
            ax.plot(xs, vs, linewidth=4, label=str(paths_list[i]))
        plt.legend()
        return plt


def main():

    ie = IonExplorer()
    ie.calc()
    

if __name__ == '__main__':
    main()



