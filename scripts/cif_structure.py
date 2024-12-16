import numpy as np
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from ase.io import read, write
import pandas as pd

class CifStructure:
    def __init__(self, ase_structure=None, df=None, graph=None, distances=None):
        """
        Initialize the MainStructure object.
        
        Parameters:
            ase_structure (Atoms): An ASE Atoms object.
            df (pd.DataFrame): DataFrame of atomic positions/metadata.
            graph (nx.Graph): Graph representing atomic connectivity.
        """
        self.ase_structure = ase_structure
        self.df = df
        self.graph = graph if graph else nx.Graph()
        self.distances = distances

    @classmethod
    def from_cif(cls, cif_path):
        """Initialize MainStructure by reading a CIF file."""
        ase_struc = read(cif_path)
        df = cls.structure_to_df(ase_struc)
        return cls(ase_struc, df)

    @classmethod
    def from_dataframe(cls, df, cell=None, pbc=(True,True,True)):
        """Initialize MainStructure from a DataFrame (with optional cell info)."""
        ase_struc = cls.df_to_structure(df, cell=cell, pbc=pbc)
        return cls(ase_struc, df)

    @staticmethod
    def structure_to_df(ase_structure, res_name='SIO'):
        """Convert an ASE structure to a DataFrame with columns ['x','y','z','atom_name','res_name']."""
        df = pd.DataFrame(ase_structure.get_positions(), columns=['x', 'y', 'z'])
        df['atom_name'] = ase_structure.get_chemical_symbols()
        df['res_name'] = res_name
        return df
    
    def _sort_dataframe(self, key='default'):
        if key != 'default':
            self.df = self.df.sort_values(key)
        self.df['id'] = range(len(self.df))

    @staticmethod
    def df_to_structure(df, cell=None, pbc=(True,True,True)):
        """Create an ASE Atoms object from a DataFrame."""
        symbols = df['atom_name'].tolist()
        positions = df[['x','y','z']].values
        structure = Atoms(symbols=symbols, positions=positions)
        if cell is not None:
            structure.set_cell(cell)
            structure.set_pbc(pbc)
        return structure

    def make_cif_name(self):
        self._sort_dataframe()
        # Use pandas to calculate row numbers (cif_id) for each atom_name group
        self.df['cif_id'] = self.df.groupby('atom_name').cumcount() + 1
        
        # Create the cif_name column
        self.df['cif_name'] = self.df['atom_name'] + self.df['cif_id'].astype(str)

        self.node2cif_dict = self.df.set_index('node_name')['cif_name'].to_dict()

    def write_cif(self, filename):
        """
        Write CIF file with bond information. 
        For simplicity, you might just write the ASE structure and separately handle bond info 
        (like the _geom_bond_atom_site_label lines).
        """
        write(filename, self.ase_structure)


    def write_cif_with_bonds(self, filename):
        """
        Write CIF file with bond information. 
        For simplicity, you might just write the ASE structure and separately handle bond info 
        (like the _geom_bond_atom_site_label lines).
        """
        # 1. Write the structure itself
        write(filename, self.ase_structure)
        cif_string = self.make_cif_string()
        with open(filename, 'a') as f:
            f.write(cif_string)

    def make_cif_string(self):
        self.make_cif_name()
        s = """loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
"""

        for i, j in self.distances:
            s += f'{self.node2cif_dict[i]:10s} {self.node2cif_dict[j]:10s} {self.distances[(i, j)]:10.4f}\n'
        print(s)
        return s        
