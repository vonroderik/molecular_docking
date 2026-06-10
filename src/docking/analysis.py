import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign

def extract_vina_score(log_path: Path):
    try:
        with open(log_path, "r") as f:
            for line in f:
                if line.startswith("   1 "):
                    return float(line.split()[1])
    except Exception:
        return None

def analyze_results(docked_pdbqt: Path, reference_pdb: Path, results_dir: Path):
    sdf_out = results_dir / "docked_poses.sdf"
    
    # Exporta para SDF usando o meeko CLI
    os.system(f"mk_export {docked_pdbqt} -s {sdf_out}")
    
    if not sdf_out.exists():
        return None, "Falha ao gerar o arquivo SDF."

    try:
        ref_mol = Chem.MolFromPDBFile(str(reference_pdb), removeHs=True)
        suppl = Chem.SDMolSupplier(str(sdf_out), removeHs=True)
        best_pose = suppl[0]
        
        if not best_pose:
            return None, "RDKit falhou em ler a pose do SDF."
            
        rmsd = rdMolAlign.GetBestRMS(best_pose, ref_mol)
        return rmsd, None
    except Exception as e:
        return None, str(e)
