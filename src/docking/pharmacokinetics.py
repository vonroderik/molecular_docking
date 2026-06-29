from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


def calculate_admet_descriptors(ligand_sdf: Path) -> dict:
    """
    Calcula descritores farmacocinéticos (ADMET) baseados em filtros moleculares clássicos
    (Regra de Cinco de Lipinski e Regras de Veber) usando o RDKit de forma local.
    Analisa a primeira pose válida contida no arquivo SDF fornecido.
    """
    ligand_sdf = Path(ligand_sdf)
    if not ligand_sdf.exists():
        raise FileNotFoundError(f"Arquivo SDF do ligante não encontrado em: {ligand_sdf}")

    try:
        suppl = Chem.SDMolSupplier(str(ligand_sdf))
        mol = None
        for m in suppl:
            if m is not None:
                mol = m
                break
    except Exception as e:
        raise RuntimeError(f"Erro ao ler o arquivo SDF com o RDKit: {e}")

    if mol is None:
        raise ValueError(f"RDKit não conseguiu identificar uma molécula válida no arquivo: {ligand_sdf}")

    try:
        # Cálculo dos descritores físico-químicos com RDKit
        mw = float(Descriptors.ExactMolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        hbd = int(Descriptors.NumHDonors(mol))
        hba = int(Descriptors.NumHAcceptors(mol))
        tpsa = float(Descriptors.TPSA(mol))
        rotb = int(Descriptors.NumRotatableBonds(mol))
    except Exception as e:
        raise RuntimeError(f"Erro ao calcular os descritores moleculares com RDKit: {e}")

    # Validação da Regra de Cinco de Lipinski
    # Parâmetros: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10. Permite-se no máximo 1 violação.
    lipinski_violations = []
    if mw > 500.0:
        lipinski_violations.append(f"Peso Molecular elevado ({mw:.2f} > 500)")
    if logp > 5.0:
        lipinski_violations.append(f"LogP elevado ({logp:.2f} > 5)")
    if hbd > 5:
        lipinski_violations.append(f"Doadores de H em excesso ({hbd} > 5)")
    if hba > 10:
        lipinski_violations.append(f"Aceitadores de H em excesso ({hba} > 10)")

    lipinski_pass = len(lipinski_violations) <= 1

    # Validação das Regras de Veber
    # Parâmetros: Ligações Rotacionáveis <= 10, TPSA <= 140. Não são permitidas violações.
    veber_violations = []
    if rotb > 10:
        veber_violations.append(f"Ligações rotacionáveis em excesso ({rotb} > 10)")
    if tpsa > 140.0:
        veber_violations.append(f"TPSA elevado ({tpsa:.2f} > 140)")

    veber_pass = len(veber_violations) == 0

    pass_filters = lipinski_pass and veber_pass

    # Estruturação e formatação com duas casas decimais
    return {
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hydrogen_bond_donors": hbd,
        "hydrogen_bond_acceptors": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotb,
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": lipinski_pass,
        "veber_violations": veber_violations,
        "veber_pass": veber_pass,
        "pass_filters": pass_filters,
    }
