import urllib.request
from pathlib import Path
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
import os

def download_pdb(pdb_id: str, raw_dir: Path) -> Path:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_file = raw_dir / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        urllib.request.urlretrieve(url, pdb_file)
    return pdb_file

def download_pubchem_sdf(cid: str, out_path: Path) -> Path:
    """Baixa um composto do PubChem em formato SDF 3D."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/record/SDF/?record_type=3d"
    if not out_path.exists():
        urllib.request.urlretrieve(url, out_path)
    return out_path

def split_receptor_ligand(pdb_file: Path, processed_dir: Path):
    receptor_path = processed_dir / "receptor.pdb"
    ligand_path = processed_dir / "ligand.pdb"
    
    with open(pdb_file, "r") as fin, \
         open(receptor_path, "w") as frec, \
         open(ligand_path, "w") as flig:
        for line in fin:
            if line.startswith("ATOM"):
                frec.write(line)
            elif line.startswith("HETATM") and "HOH" not in line and "WAT" not in line:
                flig.write(line)
            elif line.startswith("TER") or line.startswith("END"):
                frec.write(line)
    return receptor_path, ligand_path

def prepare_receptor(pdb_path: Path, out_path: Path):
    """
    Prepara o receptor (PDBQT) a partir de um PDB.
    Utiliza o mk_prepare_receptor do Meeko diretamente.
    """
    # Usamos o arquivo PDB original gerado pelo split
    # Adicionamos flags para lidar com altlocs e resíduos problemáticos
    cmd = (
        f"mk_prepare_receptor --read_pdb {pdb_path} "
        f"-o {out_path.with_suffix('')} "
        f"--default_altloc A --allow_bad_res -p"
    )
    
    resultado = os.system(cmd)
    
    if resultado != 0 or not out_path.exists():
        raise RuntimeError(f"Falha na conversão do receptor via mk_prepare_receptor. Comando: {cmd}")

def prepare_ligand_sdf(sdf_path: Path, out_path: Path):
    """Prepara o ligante (PDBQT) a partir de um SDF usando mk_prepare_ligand."""
    cmd = f"mk_prepare_ligand -i {sdf_path} -o {out_path}"
    resultado = os.system(cmd)
    if resultado != 0:
        raise RuntimeError(f"Falha na conversão do ligante via mk_prepare_ligand. Comando: {cmd}")

def prepare_ligand(pdb_path: Path, out_path: Path):
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if not mol:
        raise ValueError("RDKit não conseguiu processar o PDB do ligante.")
    
    frags = Chem.GetMolFrags(mol, asMols=True)
    # Isola o maior fragmento (assume ser o ligante de interesse)
    ligand = max(frags, key=lambda m: m.GetNumAtoms())
    ligand = Chem.AddHs(ligand)

    prep = MoleculePreparation()
    setups = prep.prepare(ligand)
    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(setups[0])
    
    with open(out_path, "w") as f:
        f.write(pdbqt_string)
