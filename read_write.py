from CifFile import ReadCif
from structure import *
from structure_data import StructureData
import os


class TextWrapper:

    def __init__(self, text):

        self.text = text

    def read(self):

        return self.text


class CIFReader:

    def __init__(self, file_name):

        if not os.path.isfile(file_name):
            raise ValueError("The file is not found!")
        self.file_name = file_name

    def read_cif_blocks(self):

        with open(self.file_name, "r") as open_file:
            cif_block = ""
            for l in open_file:
                if l.startswith("data_"):
                    cif_block = l
                    for l in open_file:
                        if l.startswith("#") or l.strip == "":
                            continue
                        if l.startswith("data_"):
                            yield cif_block
                            cif_block = ""
                        cif_block += l
        yield cif_block

    def read_cif(self):

        for i, cif_block in enumerate(self.read_cif_blocks()):
            structure = None
            try:
                block_name, data = list(ReadCif(TextWrapper(cif_block)).items())[0]
                structure_data = StructureData(block_name, data)
                structure = Structure().build_structure(structure_data)
            except Exception as exc:
                print("Parsing of structure {} is failed!".format(i), str(exc))
            yield structure


def read_file(file_name, out_list=True):

    with open(file_name, 'r') as open_file:
        if out_list:
            text = open_file.readlines()
        else:
            text = open_file.read()
    return text


def write_file(file_name, text, key='w'):

    with open(file_name, key) as new_file:
        if text is list:
            new_file.writelines(text)
        else:
            new_file.write(text)
    return None
