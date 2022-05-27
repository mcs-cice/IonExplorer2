import re
import numpy as np
import symmetry as sm

########################################################
############            CIF TAGS            ############
########################################################

############           CELL TAGS            ############
length_a = ["_cell_length_a", "_cell.length_a"]
length_b = ["_cell_length_b", "_cell.length_b"]
length_c = ["_cell_length_c", "_cell.length_c"]
angle_alpha = ["_cell_angle_alpha", "_cell.angle_alpha"]
angle_beta = ["_cell_angle_beta", "_cell.angle_beta"]
angle_gamma = ["_cell_angle_gamma", "_cell.angle_gamma"]

###########           SYMMETRY TAGS           ###########
symop = ["_symmetry_equiv_pos_as_xyz",
         "_space_group_symop_operation_xyz",
         "_space_group_symop.operation_xyz"]
sp_group = ["_symmetry_space_group_name_H-M"]

###########             ATOM TAGS             ############
index = ["_atom_site.id", "_atom_site_id "]
symbol = ["_atom_site_type_symbol", "_atom_site.type_symbol"]
label = ["_atom_site_label", "_atom_site.label"]
fract_x = ["_atom_site_fract_x", "_atom_site.fract_x"]
fract_y = ["_atom_site_fract_y", "_atom_site.fract_y"]
fract_z = ["_atom_site_fract_z", "_atom_site.fract_z"]
occupancy = ["_atom_site_occupancy", "_atom_site.occupancy"]

###########          CONNECTIVITIES           ############
label_1 = ["_topol_link.node_label_1"]
label_2 = ["_topol_link.node_label_2"]
symmetry_1 = ["_topol_link.site_symmetry_symop_1"]
symmetry_2 = ["_topol_link.site_symmetry_symop_2"]
translation_1_x = ["_topol_link.site_symmetry_translation_1_x"]
translation_1_y = ["_topol_link.site_symmetry_translation_1_y"]
translation_1_z = ["_topol_link.site_symmetry_translation_1_z"]
translation_2_x = ["_topol_link.site_symmetry_translation_2_x"]
translation_2_y = ["_topol_link.site_symmetry_translation_2_y"]
translation_2_z = ["_topol_link.site_symmetry_translation_2_z"]
bond_type = ["_topol_link.type"]
############               END                ############


class StructureData:

    def __init__(self, block_name, data):

        self.block_name = block_name
        self.data = data
        # Unit cell parameters
        self.cell_length_a = self._pars_num_param(*self._extract(length_a, 'cell_length_a - is not found!'))
        self.cell_length_b = self._pars_num_param(*self._extract(length_b, 'cell_length_b - is not found!'))
        self.cell_length_c = self._pars_num_param(*self._extract(length_c, 'cell_length_c - is not found!'))
        self.cell_angle_alpha = self._pars_num_param(*self._extract(angle_alpha, 'cell_angle_alpha - is not found!'))
        self.cell_angle_beta = self._pars_num_param(*self._extract(angle_beta, 'cell_angle_beta - is not found!'))
        self.cell_angle_gamma = self._pars_num_param(*self._extract(angle_gamma, 'cell_angle_gamma - is not found!'))
        # Symmetry
        self.space_group_name = self._extract_2(sp_group)[0]
        self.symmetry_equiv_pos_as_xyz = self._extract(symop, "Symmetry operations are not found!")[0]
        self.symmetry_rotations = []
        self.symmetry_translations = []
        for sym in self.symmetry_equiv_pos_as_xyz:
            rot, trans = sm.Symmetry.pars_sym_code(sym)
            self.symmetry_rotations.append(rot)
            self.symmetry_translations.append(trans)
        # Atoms
        self.atom_index = []
        self.atom_site_dis = []
        self.atom_site_type_symbol = []
        self.atom_site_label = []
        atom_site_type_symbol = self._extract_2(symbol)[0]
        atom_site_label = self._extract_2(label)[0]
        if atom_site_type_symbol is not None:
            self._pars_symbols(atom_site_type_symbol, "Symbol")
            self._pars_labels(atom_site_label, "Label")
        elif atom_site_label is None:
            raise Exception("The atomic symbols and labels are not found!")
        else:
            self._pars_labels(atom_site_label, "Label")
        x = self._extract(fract_x, '_atom_site_fract_x - is not found!')[0]
        y = self._extract(fract_y, '_atom_site_fract_y - is not found!')[0]
        z = self._extract(fract_z, '_atom_site_fract_z - is not found!')[0]
        x = [self._pars_num_param(x, "_atom_site_fract_x " + str(i + 1)) for i, x in enumerate(x)]
        y = [self._pars_num_param(y, "_atom_site_fract_y " + str(i + 1)) for i, y in enumerate(y)]
        z = [self._pars_num_param(z, "_atom_site_fract_z " + str(i + 1)) for i, z in enumerate(z)]
        self.atom_site_fract = np.array([[x[i], y[i], z[i]] for i in range(len(x))])
        self.atom_site_occupancy = self._get_param(occupancy)
        if self.atom_site_occupancy is None:
            self.atom_site_occupancy = [1.0] * len(self.atom_site_type_symbol)
        # Connectivity
        self.indexes_1 = []
        self.indexes_2 = []
        self.symmetry_1 = []
        self.symmetry_2 = []
        self.translation_1 = []
        self.translation_2 = []
        self.bond_type = []
        try:
            if self.atom_site_label is None:
                raise Exception("The atomic symbols and labels are not found!")
            _label_1 = self._extract(label_1, 'Atomic label_1 are not found!')[0]
            _label_2 = self._extract(label_2, 'Atomic label_2 are not found!')[0]
            dict_labels = {l: i for i, l in enumerate(self.atom_site_label)}
            self.indexes_1 = [dict_labels[l] for l in _label_1]
            self.indexes_2 = [dict_labels[l] for l in _label_2]
            self.symmetry_1 = [int(el) for el in self._extract(symmetry_1, 'symmetry_1 - is not found!')[0]]
            self.symmetry_2 = [int(el) for el in self._extract(symmetry_2, 'symmetry_2 - is not found!')[0]]
            t_1_x = [int(el) for el in self._extract(translation_1_x, 'translation_1_x - is not found!')[0]]
            t_1_y = [int(el) for el in self._extract(translation_1_y, 'translation_1_y - is not found!')[0]]
            t_1_z = [int(el) for el in self._extract(translation_1_z, 'translation_1_z - is not found!')[0]]
            t_2_x = [int(el) for el in self._extract(translation_2_x, 'translation_2_x - is not found!')[0]]
            t_2_y = [int(el) for el in self._extract(translation_2_y, 'translation_2_y - is not found!')[0]]
            t_2_z = [int(el) for el in self._extract(translation_2_z, 'translation_2_z - is not found!')[0]]
            self.translation_1 = [[t_1_x[i], t_1_y[i], t_1_z[i]] for i in range(len(t_1_x))]
            self.translation_2 = [[t_2_x[i], t_2_y[i], t_2_z[i]] for i in range(len(t_2_x))]
            self.bond_type = self._extract(bond_type, 'Bond_type - is not found!')[0]
            if (
                    not
                    (
                            len(_label_1) == len(_label_2) == len(self.symmetry_1) == len(self.symmetry_2)
                            == len(self.translation_1) == len(self.translation_2) == len(self.bond_type)
                    )
            ):
                raise Exception("...!!!!!!!!!!!!!!")
        except Exception as e:
            #print(str(e) + " Impossible to read the atomic connectivity!")
            pass

    def _extract(self, tags, exception=''):

        for tag in tags:
            param = self.data.get(tag)
            if param is not None:
                return param, tag
        raise Exception(exception)

    def _extract_2(self, tags):

        for tag in tags:
            param = self.data.get(tag)
            if param is not None:
                return param, tag
        return None, None

    def _get_param(self, tags):
        
        for tag in tags:
            param = self.data.get(tag)
            if param is not None:
                return param
        return None

    def _pars_num_param(self, param, tag):

        r_number = re.compile("\s*(-?\d*)(\.\d?)?\(?(\d+)?\)?\s*\w*")
        m_number = r_number.match(param)
        if m_number:
            param = float(''.join([g if g is not None else '' for g in m_number.groups()]))
        else:
            raise Exception(tag + " - invalid value : " + param)
        return param

    def _pars_symbols(self, symbols, tag):

        r_symbol = re.compile("([a-zA-Z]+)([^a-zA-Z]+.*)?")
        for i, symbol in enumerate(symbols):
            m_symbol = r_symbol.match(symbol)
            if m_symbol:
                self.atom_site_type_symbol.append(m_symbol.group(1))
            else:
                raise Exception(tag + ' ' + str(i) + " - invalid value : " + symbol)
        return None
        
    def _pars_labels(self, labels, tag):

        r_label = re.compile("\d*([a-zA-Z]+)_?(\d+)?_?([a-zA-Z]+)?")
        extract_symbols = False
        if self.atom_site_type_symbol is None:
            extract_symbols = True
            self.atom_site_type_symbol = []
        for i, l in enumerate(labels):
            m_label = r_label.match(l)
            if m_label:
                self.atom_site_label.append(l)
                self.atom_site_dis.append(m_label.group(3))
                if extract_symbols:
                    self.atom_site_type_symbol.append(m_label.group(1))
            else:
                raise Exception(tag + ' ' + str(i) + " - invalid value : " + l)
        return None

    def __str__(self):

        data = (
                "data_" + self.block_name + '\n' +
                "_cell_length_a\t" + '%8.5f' % self.cell_length_a + '\n' +
                "_cell_length_b\t" + '%8.5f' % self.cell_length_b + '\n' +
                "_cell_length_c\t" + '%8.5f' % self.cell_length_c + '\n' +
                "_cell_angle_alpha\t" + '%5.2f' % self.cell_angle_alpha + '\n' +
                "_cell_angle_beta\t" + '%5.2f' % self.cell_angle_beta + '\n' +
                "_cell_angle_gamma\t" + '%5.2f' % self.cell_angle_gamma + '\n' +
                "loop_\n" +
                "_symmetry_equiv_pos_as_xyz\n" +
                ''.join(["\'" + s + "\'" + "\n" for s in self.symmetry_equiv_pos_as_xyz]) +
                "loop_\n" +
                "_atom_site_label\n" +
                "_atom_site_type_symbol\n" +
                "_atom_site_fract_x\n" +
                "_atom_site_fract_y\n" +
                "_atom_site_fract_z\n" +
                ''.join(
                    [self.atom_site_label[i] + '\t' +
                         self.atom_site_type_symbol[i] +
                         '%10.5f' % self.atom_site_fract[i][0] +
                         '%10.5f' % self.atom_site_fract[i][1] +
                         '%10.5f' % self.atom_site_fract[i][2] + '\n'
                         for i in range(len(self.atom_site_label))]
                )
        )
        return data
