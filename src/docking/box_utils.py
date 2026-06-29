import numpy as np
from rdkit import Chem
from pathlib import Path


def calculate_centroid(ligand_pdb: Path):
    mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
    if not mol:
        raise ValueError("Não foi possível ler o PDB do ligante para calcular o box.")

    conf = mol.GetConformer()
    coords = conf.GetPositions()
    center = np.mean(coords, axis=0)

    return {
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "center_z": float(center[2]),
        "size_x": 22.0,
        "size_y": 22.0,
        "size_z": 22.0,
    }
