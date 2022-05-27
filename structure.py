from symmetry import*
from geometry import*
from constants import*
import numpy as np
from functions import remove_items
from graph import PeriodicGraph


class Atom:

    def __init__(
            self,
            atomic_number,
            index,
            fract_coord,
            cart_coord,
            symmetry=1,
            site=(0, 0, 0),
            occupancy=1.0,
            multiplicity=None,
            properties=None
    ):
        self.name = None
        self.element = PERIODIC_TABLE[atomic_number]
        self.index = index
        self.fract_coord = np.array(fract_coord)
        self.cart_coord = np.array(cart_coord)
        self.symmetry = symmetry
        self.site = np.array(site, dtype=np.int8)
        self.translation = np.array([0, 0, 0], dtype=np.int8)
        self.occupancy = occupancy
        self.multiplicity = multiplicity
        self.properties = properties
        self.equal = self
        self.assign_name()

    def copy(self):

        atom_copy = Atom(
            self.element.number,
            self.index,
            self.fract_coord,
            self.cart_coord,
            self.symmetry,
            self.site,
            self.occupancy,
            self.multiplicity,
            self.properties
        )
        atom_copy.translation = np.array(self.translation)
        atom_copy.equal = self.equal
        return atom_copy

    @staticmethod
    def give_name(atomic_number, index, symmetry, site):

        name = ((atomic_number, index, symmetry), tuple(site))
        return name

    def assign_name(self):

        self.name = self.give_name(self.element.number, self.index, self.symmetry, self.site)
        return self.name

    def __str__(self):

        return str(self.name)


class Bond:

    def __init__(self, atom_1, atom_2, type='V', length=None, multiplicity=None, order=None, weight=None):

        self.name = None
        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.type = type
        self.length = length
        self.multiplicity = multiplicity
        self.order = order
        self.weight = weight
        self.equal = self
        self.assign_name()

    def copy(self):

        bond_copy = Bond(
            self.atom_1,
            self.atom_2,
            self.type,
            self.length,
            self.multiplicity,
            self.order,
            self.weight
        )
        bond_copy.equal = self.equal
        return bond_copy

    def cal_length(self):

        return np.linalg.norm(np.array(self.atom_1.cart_coord) - self.atom_2.cart_coord)

    def assign_name(self):

        self.name = (self.atom_1.name, self.atom_2.name, self.type)
        return self.name

    def __str__(self):

        text = str(self.name)
        return text


class Cell:

    def __init__(self, a, b, c, alpha, beta, gamma):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.vectors = None
        self.inv_vectors = None
        self.volume = None
        self.orthogonal = self.is_orthogonal()
        self.calc_cell_vectors()

    def is_orthogonal(self, tol=1e-2):

        if abs(self.alpha - 90) < tol and abs(self.beta - 90) < tol and abs(self.gamma - 90) < tol:
            return True
        return False

    def calc_cell_vectors(self):

        a, b, c = self.a, self.b, self.c
        if self.orthogonal:
            self.vectors = np.array(
                [[a, 0, 0],
                 [0, b, 0],
                 [0, 0, c]]
            )
        else:
            alpha = np.radians(self.alpha % 180)
            betta = np.radians(self.beta % 180)
            gamma = np.radians(self.gamma % 180)
            c1 = c * np.cos(betta)
            c2 = c * (np.cos(alpha) - np.cos(gamma) * np.cos(betta)) / np.sin(gamma)
            c3 = np.sqrt(c * c - c1 * c1 - c2 * c2)
            self.vectors = np.array([[a, 0., 0.],
                                     [b * np.cos(gamma), b * np.sin(gamma), 0.],
                                     [c1, c2, c3]])
            self.inv_vectors = np.linalg.inv(self.vectors)
        return self.vectors

    def get_cart_coord(self, fract_coord):

        if self.orthogonal:
            cart_coord = fract_coord * [self.a, self.b, self.c]
        else:
            cart_coord = (
                    self.vectors[0] * fract_coord[0]
                    + self.vectors[1] * fract_coord[1]
                    + self.vectors[2] * fract_coord[2]
            )
        return cart_coord

    def get_fract_coord(self, cart_coord):

        if self.orthogonal:
            fract_coord = cart_coord / [self.a, self.b, self.c]
        else:
            fract_coord = (
                    self.inv_vectors[0] * cart_coord[0]
                    + self.inv_vectors[1] * cart_coord[1]
                    + self.inv_vectors[2] * cart_coord[2]
            )
        return fract_coord

    def cal_volume(self):

        self.volume = np.dot(self.vectors[0], np.cross(self.vectors[1], self.vectors[2]))
        return self.volume

    def copy(self):

        copy = Cell(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        return copy

    def extend(self, ka, kb, kc):

        self.a = ka * self.a
        self.b = kb * self.b
        self.c = kc * self.c
        self.calc_cell_vectors()
        return None


class Structure:

    def __init__(self):

        self.name = None
        self.structure_data = None
        self.cell = None
        self.symmetry = None
        self.equiv_positions = {}
        self.elements_count = {}
        self.nonequiv_atoms = {}
        self.equiv_atoms = {}
        self.atoms = {}
        self.nonequiv_bonds = {}
        self.equiv_bonds = {}
        self.bonds = {}
        self.graph = None

    def build_structure(self, structure_data):

        self.structure_data = structure_data
        # Unit cell
        self.cell = Cell(
            structure_data.cell_length_a,
            structure_data.cell_length_b,
            structure_data.cell_length_c,
            structure_data.cell_angle_alpha,
            structure_data.cell_angle_beta,
            structure_data.cell_angle_gamma
        )
        # Symmetry
        self.symmetry = Symmetry(structure_data.symmetry_rotations, structure_data.symmetry_translations)
        # Atoms
        atoms_index = []
        n_atoms = len(self.nonequiv_atoms)
        for i, symbol in enumerate(structure_data.atom_site_type_symbol):
            atomic_number = get_number(symbol)
            self.add_atom(atomic_number, structure_data.atom_site_fract[i])
            if n_atoms != len(self.nonequiv_atoms):
                atoms_index += [self.elements_count[symbol]]
                n_atoms = len(self.nonequiv_atoms)
        # Bonds
        for i in range(len(structure_data.symmetry_1)):
            atom_1 = self.get_equal_atom(
                get_number(structure_data.atom_site_type_symbol[structure_data.indexes_1[i]]),
                atoms_index[structure_data.indexes_1[i]],
                structure_data.symmetry_1[i],
                structure_data.translation_1[i]
            )
            atom_2 = self.get_equal_atom(
                get_number(structure_data.atom_site_type_symbol[structure_data.indexes_2[i]]),
                atoms_index[structure_data.indexes_2[i]],
                structure_data.symmetry_2[i],
                structure_data.translation_2[i]
            )

            bond_type = structure_data.bond_type[i]
            self.add_bond(atom_1, atom_2, bond_type)
        return self

    def copy(self):

        copy = Structure()
        copy.name = self.name
        copy.structure_data = self.structure_data
        copy.cell = self.cell.copy()
        copy.symmetry = self.symmetry.copy()
        copy.equiv_positions = {k: {k_2: v_2 for k_2, v_2 in v.items()} for k, v in self.equiv_positions.items()}
        copy.elements_count = {k: v for k, v in self.elements_count.items()}
        # Copy of atoms
        copy.nonequiv_atoms = {n: a.copy() for n, a in self.nonequiv_atoms.items()}
        copy.equiv_atoms = {n: a for n, a in copy.nonequiv_atoms.items()}
        copy.equiv_atoms.update({n: a.copy() for n, a in self.equiv_atoms.items() if n[0][2] != 1})
        copy.atoms = {n: a for n, a in copy.equiv_atoms.items()}
        copy.atoms.update({n: a.copy() for n, a in self.atoms.items() if (n[1] != np.array((0, 0, 0))).any()})
        for a in copy.nonequiv_atoms.values():
            a.equal = a
        # Copy of bonds
        bonds = list(self.bonds.values())
        bonds.sort(key=lambda x: 0 if x.name == x.equal.name else 1)
        for b in bonds:
            bond_copy = b.copy()
            bond_copy.atom_1 = copy.get_atom(bond_copy.atom_1.name)
            bond_copy.atom_2 = copy.get_atom(bond_copy.atom_2.name)
            if b.name == b.equal.name:
                bond_copy.equal = bond_copy
                copy.nonequiv_bonds[b.name] = bond_copy
                copy.equiv_bonds[b.name] = bond_copy
                copy.bonds[b.name] = bond_copy
            elif (b.atom_1.site == np.array((0, 0, 0))).all():
                bond_copy.equal = copy.get_bond(b.equal.name)
                copy.equiv_bonds[b.name] = bond_copy
            else:
                bond_copy.equal = copy.get_bond(b.equal.name)
            copy.bonds[b.name] = bond_copy
        return copy

    def to_graph(self, bond_types=None):

        atoms = list(self.equiv_atoms.values())
        atoms_dict = {a.name[0]: i for i, a in enumerate(atoms)}
        bonds = list(self.equiv_bonds.values())
        edges = []
        translations = []
        for b in bonds:
            if bond_types is None or b.type in bond_types:
                translations.append(b.atom_2.site - b.atom_1.site)
                edges.append((atoms_dict[b.atom_1.name[0]],  atoms_dict[b.atom_2.name[0]]))
        self.graph = PeriodicGraph(len(atoms), edges, translations, {"atoms": atoms}, {"bonds": bonds})
        return self.graph, atoms_dict, {b.name: i for i, b in enumerate(bonds)}

    def increment_element_index(self, symbol):

        index = self.elements_count.get(symbol)
        if index is None:
            index = 1
        else:
            index += 1
        self.elements_count[symbol] = index
        return index

    def multiply_atom(self, atom, tol=0.01):

        equiv_positions = {1: (1, np.array([0, 0, 0]))}
        equiv_atoms = [atom]
        multiplicity = 1
        for i, symop in enumerate(self.symmetry.symm_operations[1:]):
            new_fract_coord = self.symmetry.apply_symmetry(atom.fract_coord, symop)
            new_fract_coord, translation = self.symmetry.calc_coord_and_translation(new_fract_coord)
            new_cart_coord = self.cell.get_cart_coord(new_fract_coord)
            for j, equiv_atom in enumerate(equiv_atoms):
                if np.linalg.norm(new_cart_coord - equiv_atom.cart_coord) < tol:
                    if equiv_atom.symmetry != i + 2:
                        equiv_positions[i + 2] = (equiv_atoms[j].symmetry, equiv_atoms[j].translation - translation)
                    break
            else:
                equiv_atom = atom.copy()
                equiv_atom.fract_coord = new_fract_coord
                equiv_atom.cart_coord = new_cart_coord
                equiv_atom.symmetry = i + 2
                equiv_atom.translation = translation
                equiv_atom.assign_name()
                equiv_atoms += [equiv_atom]
                equiv_positions[i + 2] = (equiv_atom.symmetry, np.array([0, 0, 0]))
                multiplicity += 1
        atom.multiplicity = multiplicity
        self.equiv_atoms.update([(a.name, a) for a in equiv_atoms])
        self.atoms.update(self.equiv_atoms)
        self.equiv_positions[(atom.element.number, atom.index)] = equiv_positions
        return None

    def find_duplicate(self, fract_coord, tol=0.001):

        fract_coord, translation = Symmetry.calc_coord_and_translation(fract_coord)
        atoms = [a for a in self.equiv_atoms.values()]
        if len(atoms) != 0:
            dist = cdist([self.cell.get_cart_coord(fract_coord)], [a.cart_coord for a in atoms])[0]
            min_dist = min(dist)
            if min_dist < tol:
                #return self.get_atom((atoms[np.where(dist < min_dist + 1e-8)[0][0]].name[0], tuple(-translation)))
                return atoms[np.where(dist < min_dist + 1e-8)[0][0]]

    def add_atom(self, atomic_number, fract_coord, symmetry=1, properties=None, tol=1e-4):

        if tol > 0:
            duplicate = self.find_duplicate(fract_coord, tol)
            if duplicate is not None:
                return duplicate

        fract_coord = Symmetry.calc_coord_and_translation(fract_coord)[0]
        index = self.increment_element_index(PERIODIC_TABLE[atomic_number].symbol)
        cart_coord = self.cell.get_cart_coord(fract_coord)
        atom = Atom(atomic_number, index, fract_coord, cart_coord, symmetry=symmetry, properties=properties)
        self.nonequiv_atoms[atom.name] = atom
        self.multiply_atom(atom)
        return atom

    def get_atom(self, atom_name):

        atom = self.atoms.get(atom_name)
        if atom is None:
            translation = atom_name[1]
            atom_name = (atom_name[0], (0, 0, 0))
            atom = self.atoms.get(atom_name)
            if atom is not None:
                atom = atom.copy()
                atom.cart_coord += (
                        translation[0] * self.cell.vectors[0]
                        + translation[1] * self.cell.vectors[1]
                        + translation[2] * self.cell.vectors[2]
                )
                atom.site += translation
                atom.translation += translation
                atom.assign_name()
                self.atoms[atom.name] = atom
        return atom

    def remove_atoms(self, atoms, cond=lambda x: True):

        to_remove = {atom.name[0][:2]: True for atom in atoms if cond(atom)}
        removed_atoms = remove_items(self.atoms, lambda x: to_remove.get(x[0][0][:2]) is not None)
        remove_items(self.equiv_atoms, lambda x: to_remove.get(x[0][0][:2]) is not None)
        remove_items(self.nonequiv_atoms, lambda x: to_remove.get(x[0][0][:2]) is not None)
        r_bonds = [
            v for k, v in self.nonequiv_bonds.items()
            if to_remove.get(k[0][0][:2]) is not None or to_remove.get(k[1][0][:2]) is not None
        ]
        removed_bonds = self.remove_bonds(r_bonds)
        remove_items(self.equiv_positions, lambda x: to_remove.get(x[0]) is not None)
        return removed_atoms, removed_bonds

    def get_equal_atom(self, atomic_number, index, symmetry, translation):

        eq_symmetry, additional_translation = self.equiv_positions[(atomic_number, index)][symmetry]
        name = Atom.give_name(atomic_number, index, eq_symmetry, (0, 0, 0))
        site = translation + additional_translation - self.get_atom(name).translation
        atom = self.get_atom(Atom.give_name(atomic_number, index, eq_symmetry, site))
        return atom

    @staticmethod
    def calc_compositions(atoms, cond=lambda x: True):

        s_dict = {}
        symbols = [a.element.symbol for a in atoms]
        for s in symbols:
            if cond(s):
                if s_dict.get(s) is None:
                    s_dict[s] = 1
                else:
                    s_dict[s] += 1
        composition = list(s_dict.items())
        composition.sort(key=lambda x: get_number(x[0]))
        return composition

    @staticmethod
    def sort_atoms(atoms):
        atoms.sort(key=lambda x: (x.element.number, x.index, x.symmetry, x.site[0], x.site[1], x.site[2]))
        return atoms

    def add_bond(self, atom_1, atom_2, type='S', length=None, multiplicity=None, order=None, weight=None):

        atom_1, atom_2 = self.sort_atoms([atom_1, atom_2])
        atom_2 = self.get_atom((atom_2.name[0], tuple(atom_2.site - atom_1.site)))
        atom_1 = self.get_atom((atom_1.name[0], (0, 0, 0)))
        if self.bonds.get((atom_1.name, atom_2.name, type)) is not None:
            return None
        bond = Bond(atom_1, atom_2, type, length, multiplicity, order, weight)
        self.nonequiv_bonds[bond.name] = bond
        self.multiply_bond(bond)
        return bond

    def multiply_bond(self, bond):

        s = self.symmetry
        equiv_bonds = {}
        self.equiv_bonds[bond.name] = bond
        self.bonds[bond.name] = bond
        for i, symop in enumerate(s.symm_operations[1:]):
            symmetry_1 = s.symm_operations[bond.atom_1.symmetry - 1]
            symmetry_1 = (symmetry_1[0], symmetry_1[1] + bond.atom_1.translation)
            new_symmetry_1, translations_1 = s.get_equiv_symop(s.summarize_operations(symmetry_1, symop))
            atom_1 = self.get_equal_atom(bond.atom_1.element.number, bond.atom_1.index, new_symmetry_1, translations_1)
            symmetry_2 = s.symm_operations[bond.atom_2.symmetry - 1]
            symmetry_2 = (symmetry_2[0],  symmetry_2[1] + bond.atom_2.translation)
            new_symmetry_2, translations_2 = s.get_equiv_symop(s.summarize_operations(symmetry_2, symop))
            atom_2 = self.get_equal_atom(bond.atom_2.element.number, bond.atom_2.index, new_symmetry_2, translations_2)
            atom_1, atom_2 = self.sort_atoms([atom_1, atom_2])
            atom_2 = self.get_atom((atom_2.name[0], tuple(atom_2.site - atom_1.site)))
            atom_1 = self.get_atom((atom_1.name[0], (0, 0, 0)))
            equiv_bond = Bond(atom_1, atom_2, bond.type, bond.length)
            equiv_bond.weight = bond.weight
            equiv_bond.equal = bond.equal
            equiv_bonds[equiv_bond.name] = equiv_bond
            self.equiv_bonds[equiv_bond.name] = equiv_bond
            self.bonds[equiv_bond.name] = equiv_bond
            if abs(bond.cal_length() - equiv_bond.cal_length()) > 0.01:
                raise RuntimeError("")
        return None

    def get_bond(self, name):

        bond = self.bonds.get(name)
        if bond is None:
            atom_1 = self.get_atom(name[0])
            atom_2 = self.get_atom(name[1])
            type = name[2]
            atom_1, atom_2 = self.sort_atoms([atom_1, atom_2])
            bond = self.equiv_bonds[
                ((atom_1.name[0], (0, 0, 0)), (atom_2.name[0], tuple(atom_2.site - atom_1.site)), type)
            ]
            if atom_1.name == name[0]:
                translation = name[0][1]
            else:
                translation = name[1][1]
            if translation != (0, 0, 0):
                bond = bond.copy()
                bond.atom_1 = self.get_atom((bond.atom_1.name[0], tuple(bond.atom_1.site + translation)))
                bond.atom_2 = self.get_atom((bond.atom_2.name[0], tuple(bond.atom_2.site + translation)))
                bond.assign_name()
                self.bonds[bond.name] = bond
        return bond

    def remove_bonds(self, bonds, cond=lambda x: True):

        to_remove = {bond.equal.name: True for bond in bonds if cond(bond)}
        remove_items(self.nonequiv_bonds, lambda x: to_remove.get(x[1].equal.name))
        remove_items(self.equiv_bonds, lambda x: to_remove.get(x[1].equal.name))
        removed = remove_items(self.bonds, lambda x: to_remove.get(x[1].equal.name))
        return removed

    @staticmethod
    def sort_bonds(bonds):

        bonds.sort(key=lambda x: (
            x.atom_1.element.number, x.atom_1.index,
            x.atom_2.element.number, x.atom_2.index,
            x.atom_2.site[0],
            x.atom_2.site[1],
            x.atom_2.site[2]
        ))
        return bonds

    def remove_connectivity(self):

        self.nonequiv_bonds = {}
        self.equiv_bonds = {}
        self.bonds = {}
        return None

    def multiply_cell(self, a=(-1, 2), b=(-1, 2), c=(-1, 2)):

        atoms = []
        bonds = []
        translations = [(t1, t2, t3)
                        for t1 in range(a[0], a[1])
                        for t2 in range(b[0], b[1])
                        for t3 in range(c[0], c[1])]
        for atom in self.equiv_atoms.values():
            for t in translations:
                name = (atom.name[0], t)
                atoms.append(self.get_atom(name))
        for bond in self.equiv_bonds.values():
            for t in translations:
                name = ((bond.name[0][0], t), (bond.name[1][0], tuple(bond.atom_2.site + t)),  bond.name[2])
                bonds.append(self.get_bond(name))
        return atoms, bonds

    def transform_to_supercell(self, ka, kb, kc):

        atoms_dict = {k: v for k, v in self.atoms.items()}
        bonds_dict = {k: v for k, v in self.bonds.items()}
        if (ka, kb, kc) != (1, 1, 1):
            self.cell.extend(ka, kb, kc)
            atoms, bonds = self.multiply_cell(a=(0, ka), b=(0, kb), c=(0, kc))
            atoms_dict, bonds_dict = self.remove_symmetry()
            element_ind = {}
            for atom in atoms:
                if (atom.site != (0, 0, 0)).any():
                    element_ind[atom.name] = (self.increment_element_index(atom.element.symbol), atom)
                else:
                    atom.equel = atom
                    element_ind[atom.name] = (atom.index, atom)
            for atom in self.atoms.values():
                c, t = Symmetry.calc_coord_and_translation((atom.fract_coord + atom.site) / [ka, kb, kc], prec=1e-5)
                atom.index, atom.equel = element_ind[(atom.name[0], tuple(atom.site + t * [ka, kb, kc]))]
                atom.fract_coord = c
                atom.site = -t
                atom.translation = np.array(atom.site)
            for atom in self.atoms.values():
                atom.assign_name()
            for bond in self.nonequiv_bonds.values():
                bond.equal = bond
            for bond in self.bonds.values():
                bond.assign_name()
            self.nonequiv_atoms = {atom.name: atom for atom in atoms}
            self.equiv_atoms = {atom.name: atom for atom in atoms}
            self.atoms = {atom.name: atom for atom in self.atoms.values()}
            self.nonequiv_bonds = {bond.name: bond for bond in bonds}
            self.equiv_bonds = {bond.name: bond for bond in bonds}
            self.bonds = {bond.name: bond for bond in self.bonds.values()}
            self.equiv_positions = {(a.element.number, a.index): {1: (1, np.array([0, 0, 0]))}
                                    for a in self.nonequiv_atoms.values()}
        return atoms_dict, bonds_dict

    def remove_symmetry(self):

        atoms_dict = {}
        self.symmetry = Symmetry(np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]), np.array([[0., 0., 0.]]))
        self.symmetry.name = "P1"
        self.elements_count = {}
        atoms = self.sort_atoms(list(self.equiv_atoms.values()))
        for atom in atoms:
            atoms_dict[atom.name] = atom
            atom.index = self.increment_element_index(atom.element.symbol)
            atom.translation = np.array([0, 0, 0])
            atom.symmetry = 1
            atom.equal = atom
            atom.assign_name()
        for atom in self.atoms.values():
            if np.any(atom.site != [0, 0, 0]):
                atoms_dict[atom.name] = atom
                atom.index = atoms_dict[(atom.name[0], (0, 0, 0))].index
                atom.symmetry = 1
                atom.translation = atom.site
                atom.assign_name()
        self.atoms = {atom.name: atom for atom in self.atoms.values()}
        self.equiv_atoms = {atom.name: atom for atom in atoms}
        self.nonequiv_atoms = {atom.name: atom for atom in atoms}
        bonds_dict = {bond.name: bond for bond in self.bonds.values()}
        self.bonds = {bond.assign_name(): bond for bond in self.bonds.values()}
        self.equiv_bonds = {bond.assign_name(): bond for bond in self.equiv_bonds.values()}
        self.nonequiv_bonds = {bond.assign_name(): bond for bond in self.equiv_bonds.values()}
        self.equiv_positions = {(a.element.number, a.index): {1: (1, np.array([0, 0, 0]))}
                                for a in self.nonequiv_atoms.values()}
        return atoms_dict, bonds_dict

    def __str__(self):

        data = "data_" + self.structure_data.block_name + '\n'
        data += "_cell_length_a                      " + str(self.cell.a) + '\n'
        data += "_cell_length_b                      " + str(self.cell.b) + '\n'
        data += "_cell_length_c                      " + str(self.cell.c) + '\n'
        data += "_cell_angle_alpha                   " + str(self.cell.alpha) + '\n'
        data += "_cell_angle_beta                    " + str(self.cell.beta) + '\n'
        data += "_cell_angle_gamma                   " + str(self.cell.gamma) + '\n'
        data += "_symmetry_space_group_name_H-M      " + str(self.symmetry.name) + '\n'
        data += ("loop_\n"
                 + "_space_group_symop_id\n"
                 + "_space_group_symop_operation_xyz\n")
        for j, s in enumerate(self.symmetry.symm_operations):
            data += (str(j + 1) + ' ' + self.symmetry.get_symop_as_xyz(*s) + '\n')
        data += ("loop_\n"
                 + "_atom_site_label\n"
                 + "_atom_site_type_symbol\n"
                 + "_atom_site_fract_x\n"
                 + "_atom_site_fract_y\n"
                 + "_atom_site_fract_z\n")
        atoms = list(self.nonequiv_atoms.values())
        atoms = self.sort_atoms(atoms)
        for atom in atoms:
            data += (atom.element.symbol + str(atom.index) + ' '
                     + atom.element.symbol + ' '
                     + '%6.5f' % (atom.fract_coord[0]) + ' '
                     + '%6.5f' % (atom.fract_coord[1]) + ' '
                     + '%6.5f' % (atom.fract_coord[2]) + '\n')
        bonds = list(self.nonequiv_bonds.values())
        if len(bonds) != 0:
            data += ("loop_\n"
                     + "_topol_link.node_label_1\n"
                     + "_topol_link.node_label_2\n"
                     + "_topol_link.site_symmetry_symop_1\n"
                     + "_topol_link.site_symmetry_translation_1_x\n"
                     + "_topol_link.site_symmetry_translation_1_y\n"
                     + "_topol_link.site_symmetry_translation_1_z\n"
                     + "_topol_link.site_symmetry_symop_2\n"
                     + "_topol_link.site_symmetry_translation_2_x\n"
                     + "_topol_link.site_symmetry_translation_2_y\n"
                     + "_topol_link.site_symmetry_translation_2_z\n"
                     + "_topol_link.length\n"
                     + "_topol_link.type\n")
            bonds.sort(key=(lambda x: ((x.atom_1.element.number, x.atom_1.index),
                                       (x.atom_2.element.number, x.atom_2.index))))
            for bond in bonds:
                data += (bond.atom_1.element.symbol + str(bond.atom_1.index) + ' '
                         + bond.atom_2.element.symbol + str(bond.atom_2.index) + ' '
                         + str(bond.atom_1.symmetry) + ' '
                         + str(bond.atom_1.translation[0]) + ' '
                         + str(bond.atom_1.translation[1]) + ' '
                         + str(bond.atom_1.translation[2]) + ' '
                         + str(bond.atom_2.symmetry) + ' '
                         + str(bond.atom_2.translation[0]) + ' '
                         + str(bond.atom_2.translation[1]) + ' '
                         + str(bond.atom_2.translation[2]) + ' '
                         + str(bond.length) + ' '
                         + bond.type + '\n')
        data += "#End of " + self.structure_data.block_name + '\n'
        return data

    def to_poscar(self, fract=True, symbols=True, scale=1.0, title="Generated by IonExplorer 2.0"):

        poscar = title + '\n'
        poscar += "%2.1f" % scale + '\n'
        poscar += ("{:20.10f}{:20.10f}{:20.10f}\n".format(*self.cell.vectors[0]) +
                   "{:20.10f}{:20.10f}{:20.10f}\n".format(*self.cell.vectors[1]) +
                   "{:20.10f}{:20.10f}{:20.10f}\n".format(*self.cell.vectors[2]))
        atoms = list(self.equiv_atoms.values())
        atoms.sort(key=lambda a: (a.element.number, a.element.symbol))
        composition = self.calc_compositions(atoms)
        if symbols:
            poscar += ''.join(["{:>5s}".format(s) for s, n in composition]) + '\n'
        poscar += ''.join(["{:>5d}".format(n) for s, n in composition]) + '\n'
        if fract:
            poscar += "Direct\n"
            poscar += ''.join(["{:16.10f}{:16.10f}{:16.10f}\n".format(*a.fract_coord) for a in atoms])
        else:
            poscar += "Cartesian\n"
            poscar += ''.join(["{:16.10f}{:16.10f}{:16.10f}\n".format(*a.cart_coord) for a in atoms])
        return poscar