import re
import numpy as np
from fractions import Fraction


class Symmetry:

    def __init__(self, rotations, translations, space_group_name=None):

        self.order = len(rotations)
        self.rotations = rotations
        self.translations = translations
        trs = [self.calc_coord_and_translation(t) for t in translations]
        self.std_translations = np.array([t for t, a in trs])
        self.adt_translations = - np.array([a for t, a in trs], dtype=int)
        self.symm_operations = [[rotations[i], translations[i]] for i in range(len(rotations))]
        self.equivalent_positions = {}
        self.name = space_group_name

    def copy(self):

        copy = Symmetry([], [], self.name)
        copy.order = self.order
        copy.rotations = np.array(self.rotations)
        copy.translations = np.array(self.translations)
        copy.std_translations = np.array(self.std_translations)
        copy.adt_translations = - np.array(self.adt_translations, dtype=int)
        copy.symm_operations = [[copy.rotations[i], copy.translations[i]] for i in range(len(copy.rotations))]
        copy.equivalent_positions = {k: np.array(v) for k, v in self.equivalent_positions.items()}
        return copy

    @staticmethod
    def pars_sym_code(sym):

        r = re.compile(("(\-)?\+?(\d\/\d|\d+\.\d+|\d|[x,y,z])"
                        + "(\+|\-)?(\d\/\d|\d+\.\d+|\d|[x,y,z])?"
                        + "(\+|\-)?(\d\/\d|\d+\.\d+|\d|[x,y,z])?"
                        + "(\+|\-)?(\d\/\d|\d+\.\d+|\d|[x,y,z])?"))
        rot = np.zeros((3, 3), dtype='int')
        trans = np.zeros(3)
        for i, s in enumerate(sym.replace(' ', '').lower().split(',')):
            m = r.match(s)
            if m:
                for c in ((m.group(1), m.group(2)), (m.group(3), m.group(4)),
                          (m.group(5), m.group(6)), (m.group(7), m.group(8))):
                    sing = 1
                    if c[0] == '-':
                        sing = -1
                    if c[1] is not None and '/' in c[1]:
                        trans[i] += sing * float(c[1][0]) / float(c[1][2])
                    elif c[1] is not None and '.' in c[1]:
                        trans[i] += sing * float(c[1][0])
                    elif c[1] is not None and c[1].isnumeric():
                        trans[i] += sing * float(c[1])
                    elif c[1] is not None:
                        rot[i][ord(c[1]) - ord('x')] = sing
            else:
                raise RuntimeError("Invalid format of the symmetry code " + '<strong>' + sym + '</strong>' + '!')
        return rot, trans

    @staticmethod
    def get_symop_as_xyz(rotation, translation):

        symop_as_xyz = ""
        xyz = ["x", "y", "z"]
        for i, t in enumerate(translation):
            if abs(t) > 1e-8:
                symop_as_xyz += str(Fraction(str(t)).limit_denominator(100))
            for j, r in enumerate(rotation[i]):
                if r < 0:
                    symop_as_xyz += "-" + xyz[j]
                elif r > 0:
                    if symop_as_xyz == "" or symop_as_xyz[-1] == ",":
                        symop_as_xyz += xyz[j]
                    else:
                        symop_as_xyz += "+" + xyz[j]
            if i != 2:
                symop_as_xyz += ","
        return symop_as_xyz

    @staticmethod
    def summarize_operations(op_1, op_2):

        rotation = np.dot(op_2[0], op_1[0])
        translations = [
            (op_2[1][0] + (op_2[0][0] * op_1[1]).sum()),
            (op_2[1][1] + (op_2[0][1] * op_1[1]).sum()),
            (op_2[1][2] + (op_2[0][2] * op_1[1]).sum())
        ]
        symop = [rotation, translations]
        return symop

    def get_equiv_symop(self, symop):

        rot = symop[0]
        transl, additional_transl = self.calc_coord_and_translation(symop[1])
        symmetry_index = self.get_index([rot, transl], standardized=True)
        additional_transl += self.adt_translations[symmetry_index - 1]
        if symmetry_index is None:
            raise Exception("Parsing of symmetry parameters. The list of symmetry operations is not complete!")
        return symmetry_index, -additional_transl

    @staticmethod
    def apply_symmetry(fract_coord, symop, prec=1e-4):

        new_coord = np.dot(symop[0], fract_coord)
        new_coord = new_coord + symop[1]
        new_coord = np.array([round(c) if abs(c) < prec or abs(abs(c % 1) - 1) < prec else c for c in new_coord])
        return new_coord

    @staticmethod
    def calc_coord_and_translation(fract_coord, prec=1e-4):

        fract_coord = [round(c) if abs(c) < prec or abs(abs(c % 1) - 1) < prec else c for c in fract_coord]
        translation = np.array([0, 0, 0])
        new_coord = np.array([0., 0., 0.])
        for i in range(len(new_coord)):
            new_coord[i], translation[i] = (lambda x: (x % 1, - int(x - (x % 1))))(fract_coord[i])
        return new_coord, translation

    def get_index(self, symm_operation, standardized=False, prec=0.001):

        r_1, t_1 = symm_operation
        if standardized:
            for i, r_2 in enumerate(self.rotations):
                if (r_1 == r_2).all() and np.linalg.norm(t_1 - self.std_translations[i]) < prec:
                    return i + 1
        else:
            for i, r_2 in enumerate(self.rotations):
                if (r_1 == r_2).all() and np.linalg.norm(t_1 - self.translations[i]) < prec:
                    return i + 1
        return None