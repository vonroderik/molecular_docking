import subprocess
from pathlib import Path

def run_vina(
    vina_path: Path,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box_params: dict,
    out_pdbqt: Path,
    log_path: Path,
    exhaustiveness: int = 16
):
    # Usamos resolve().as_posix() para garantir caminhos absolutos com barras normais (/)
    # Isso evita problemas com escapes de contra-barra (\) no Windows dentro do config do Vina.
    config_content = f"""receptor = {receptor_pdbqt.resolve().as_posix()}
ligand = {ligand_pdbqt.resolve().as_posix()}

center_x = {box_params['center_x']:.3f}
center_y = {box_params['center_y']:.3f}
center_z = {box_params['center_z']:.3f}

size_x = {box_params['size_x']:.1f}
size_y = {box_params['size_y']:.1f}
size_z = {box_params['size_z']:.1f}

exhaustiveness = {exhaustiveness}
num_modes = 9
energy_range = 3
"""
    
    config_file = out_pdbqt.parent / "config.txt"
    with open(config_file, "w") as f:
        f.write(config_content)
    
    vina_cmd = [
        str(vina_path.resolve().as_posix()),
        "--config", str(config_file.resolve().as_posix()),
        "--out", str(out_pdbqt.resolve().as_posix())
    ]
    
    with open(log_path, "w") as log_file:
        subprocess.run(vina_cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    
    return config_file
