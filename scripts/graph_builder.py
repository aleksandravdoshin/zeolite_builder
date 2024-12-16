import numpy as np
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from ase.io import read, write
import pandas as pd
import duckdb

class GraphBuilder:
    """
    Handles the logic of building a graph from an ASE structure and DataFrame using covalent radii,
    ensuring no edges are created between atoms with distance 0.0.
    """

    distances = {}  # Class-level dictionary to store edge distances

    @staticmethod
    def covalent_radii():
        """Returns a dictionary of covalent radii for common elements."""
        return {
            'H': 0.31,
            'C': 0.76,
            'O': 0.66,
            'N': 0.71,
            'Si': 1.11,
            # Add more as needed
        }

    @staticmethod
    def calculate_bond_distance(atom1, atom2, tolerance=0.2):
        """Calculate the bond distance threshold for two atoms."""
        radii = GraphBuilder.covalent_radii()
        r1 = radii.get(atom1, 0.0)  # Default to 0.0 if atom is not found
        r2 = radii.get(atom2, 0.0)
        return r1 + r2 + tolerance

    @classmethod
    def build_graph(cls, df, ase_structure, atom_pairs=None, default_tolerance=0.2):
        """
        Build a graph for atom pairs based on covalent radii and tolerances.

        Parameters:
            df (pd.DataFrame): DataFrame with atomic positions and labels.
            ase_structure (Atoms): ASE structure for periodic boundary conditions.
            atom_pairs (list): List of atom pairs to consider for bonding.
                               Example: [('C', 'O'), ('C', 'C')]
            default_tolerance (float): Default tolerance to add to the bond distance.

        Returns:
            nx.Graph: Graph representing atomic connectivity.
        """
        a, b, c = np.diag(ase_structure.cell)
        radii = cls.covalent_radii()

        if atom_pairs is None:
            atom_pairs = {(a1, a2) for a1 in radii.keys() for a2 in radii.keys()}

        # Assign unique node names to the DataFrame
        if 'node_name' not in df.columns:
            df['node_name'] = [f"{row.atom_name}_{i}" for i, row in df.iterrows()]

        G = nx.Graph()
        cls.distances = {}  # Reset distances

        for atom1, atom2 in atom_pairs:
            # Get coordinates for the atom types
            array1 = df.loc[df.atom_name == atom1, ['x', 'y', 'z']].values
            array2 = df.loc[df.atom_name == atom2, ['x', 'y', 'z']].values
            if len(array1) == 0 or len(array2) == 0:
                continue
            # Compute the bond distance threshold
            bond_distance = cls.calculate_bond_distance(atom1, atom2, default_tolerance)

            # Find neighbors and distances within this distance
            neighbors_list, distances_list = cls.neighbors_within_radius(array1, array2, bond_distance, a, b, c)

            # Add edges to the graph
            for i, (nbr_indices, nbr_distances) in enumerate(zip(neighbors_list, distances_list)):
                source_node = df.loc[df.atom_name == atom1].iloc[i]['node_name']
                for neighbor_idx, dist in zip(nbr_indices, nbr_distances):
                    # Avoid self-loops or edges with distance 0.0
                    if dist == 0.0:
                        continue
                    target_node = df.loc[df.atom_name == atom2].iloc[neighbor_idx]['node_name']
                    G.add_edge(source_node, target_node)
                    cls.distances[(source_node, target_node)] = dist

        return G

    @staticmethod
    def neighbors_within_radius(array1, array2, r_c, a, b, c):
        """
        Find neighbors and distances within a radius r_c, applying 3D periodic boundary conditions.

        Returns:
            list of indices and distances for neighbors.
        """
        translations = [
            np.array([dx * a, dy * b, dz * c])
            for dx in [-1, 0, 1]
            for dy in [-1, 0, 1]
            for dz in [-1, 0, 1]
        ]
        expanded_array2 = np.vstack([array2 + t for t in translations])
        expanded_to_original = np.tile(np.arange(len(array2)), len(translations))

        nbrs = NearestNeighbors(radius=r_c, algorithm='auto').fit(expanded_array2)
        distances, expanded_indices = nbrs.radius_neighbors(array1, return_distance=True)

        original_indices = [np.unique(expanded_to_original[indices]) for indices in expanded_indices]
        original_distances = [
            distances[i][np.unique(expanded_to_original[indices], return_index=True)[1]]
            for i, indices in enumerate(expanded_indices)
        ]

        return original_indices, original_distances
