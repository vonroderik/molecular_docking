#!/usr/bin/env python

import json
from pathlib import Path

import questionary
import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from docking import analysis, box_utils, preparation, vina_runner, pharmacokinetics
from docking.preparation import get_executable

app = typer.Typer(help="Pipeline de Docking Molecular Automatizado")
console = Console()

# Definição de caminhos base
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
VINA_BIN = None


def get_vina_bin() -> Path:
    import os
    import platform
    import shutil
    import stat
    import urllib.request

    bin_dir = BASE_DIR / "bin"
    bin_dir.mkdir(exist_ok=True)

    if platform.system() == "Windows":
        vina_bin = bin_dir / "vina.exe"
        if not vina_bin.exists():
            console.print("[yellow]Baixando AutoDock Vina para Windows...[/yellow]")
            url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win64.exe"
            urllib.request.urlretrieve(url, vina_bin)
        return vina_bin
    else:
        vina_bin = bin_dir / "vina"
        if not vina_bin.exists():
            # Tenta verificar se já existe vina globalmente no PATH
            system_vina = shutil.which("vina")
            if system_vina:
                return Path(system_vina)

            console.print("[yellow]Baixando AutoDock Vina para Linux...[/yellow]")
            url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
            try:
                urllib.request.urlretrieve(url, vina_bin)
                st = os.stat(vina_bin)
                os.chmod(vina_bin, st.st_mode | stat.S_IEXEC)
            except Exception as e:
                raise RuntimeError(
                    f"AutoDock Vina para Linux não pôde ser baixado e não está instalado no PATH global. "
                    f"Erro original: {e}"
                )
        return vina_bin


def render_interactions_table(interactions: dict):
    """
    Exibe uma tabela no terminal com os aminoácidos do receptor que fizeram
    contatos estáticos (pontes de hidrogênio e contatos hidrofóbicos) na Pose 1.
    """
    table = Table(
        title="[bold magenta]Interações Estáticas Receptor-Ligante (Pose 1)[/bold magenta]",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Resíduo", style="yellow")
    table.add_column("Tipo de Interação", style="green")
    table.add_column("Distância (Å)", justify="right", style="white")

    # Adiciona as pontes de hidrogênio
    hbonds = interactions.get("hydrogen_bonds", [])
    for hb in hbonds:
        res = f"{hb['resname']} {hb['resnr']}"
        table.add_row(res, "Ponte de Hidrogênio", f"{hb['distance']:.2f}")

    # Adiciona os contatos hidrofóbicos
    hcontacts = interactions.get("hydrophobic_contacts", [])
    for hc in hcontacts:
        res = f"{hc['resname']} {hc['resnr']}"
        table.add_row(res, "Contato Hidrofóbico", f"{hc['distance']:.2f}")

    # Mensagem caso não existam interações mapeadas
    if not hbonds and not hcontacts:
        table.add_row("Nenhuma interação mapeada", "-", "-")

    console.print(table)


def render_admet_table(admet: dict):
    """
    Exibe uma tabela no terminal com os descritores farmacocinéticos (ADMET)
    e o status dos filtros moleculares clássicos (Lipinski e Veber).
    """
    if "error" in admet:
        console.print(f"[bold red]Erro ao calcular ADMET:[/bold red] {admet['error']}")
        return

    table = Table(
        title="[bold magenta]Triagem Farmacocinética (ADMET)[/bold magenta]",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Propriedade", style="yellow")
    table.add_column("Valor Calculado", justify="right", style="white")
    table.add_column("Filtro/Limite", style="blue")
    table.add_column("Status", justify="center")

    # Linha Lipinski MW
    mw = admet.get("molecular_weight", 0.0)
    mw_status = "[bold green]OK[/bold green]" if mw <= 500 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Peso Molecular (MW)", f"{mw:.2f} g/mol", "<= 500.00", mw_status)

    # Linha Lipinski LogP
    logp = admet.get("logp", 0.0)
    logp_status = "[bold green]OK[/bold green]" if logp <= 5 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Lipofilicidade (LogP)", f"{logp:.2f}", "<= 5.00", logp_status)

    # Linha Lipinski HBD
    hbd = admet.get("hydrogen_bond_donors", 0)
    hbd_status = "[bold green]OK[/bold green]" if hbd <= 5 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Doadores de H (HBD)", str(hbd), "<= 5", hbd_status)

    # Linha Lipinski HBA
    hba = admet.get("hydrogen_bond_acceptors", 0)
    hba_status = "[bold green]OK[/bold green]" if hba <= 10 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Aceitadores de H (HBA)", str(hba), "<= 10", hba_status)

    # Linha Veber TPSA
    tpsa = admet.get("tpsa", 0.0)
    tpsa_status = "[bold green]OK[/bold green]" if tpsa <= 140 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Superfície Polar (TPSA)", f"{tpsa:.2f} Å²", "<= 140.00", tpsa_status)

    # Linha Veber RotBonds
    rotb = admet.get("rotatable_bonds", 0)
    rotb_status = "[bold green]OK[/bold green]" if rotb <= 10 else "[bold red]VIOLADO[/bold red]"
    table.add_row("Ligações Rotacionáveis", str(rotb), "<= 10", rotb_status)

    console.print(table)

    # Veredito Geral
    pass_filters = admet.get("pass_filters", False)
    if pass_filters:
        veredito = "[bold white on green]  APROVADO  [/bold white on green]"
        console.print(Panel(f"Veredito de Triagem ADMET: {veredito}\nA molécula atende aos critérios clássicos de Lipinski (máximo 1 violação) e Veber (0 violações).", border_style="green"))
    else:
        veredito = "[bold white on red]  REPROVADO  [/bold white on red]"
        violacoes = admet.get("lipinski_violations", []) + admet.get("veber_violations", [])
        violacoes_str = ", ".join(violacoes) if violacoes else "Filtros não atendidos."
        console.print(Panel(f"Veredito de Triagem ADMET: {veredito}\nViolações: [red]{violacoes_str}[/red]", border_style="red"))


@app.command(name="validate")
def validate(
    pdb_id: str = typer.Option(
        "4HG7", "--pdb", help="ID do PDB para baixar e analisar"
    ),
    exhaustiveness: int = typer.Option(
        16, "--ex", help="Exaustividade do Vina (padrão: 16)"
    ),
):
    """
    POSITIVE CONTROL / REDOCKING:
    Baixa um PDB, separa o ligante nativo, prepara e executa o docking para validar o RMSD.
    """
    global VINA_BIN
    if VINA_BIN is None:
        VINA_BIN = get_vina_bin()

    run_dir = DATA_DIR / pdb_id
    raw_dir = run_dir / "raw"
    processed_dir = run_dir / "processed"
    results_dir = run_dir / "results"

    for folder in [raw_dir, processed_dir, results_dir]:
        folder.mkdir(parents=True, exist_ok=True)

    console.print(
        Panel.fit(
            f"[bold blue]Pipeline de Validação (Redocking)[/bold blue]\n"
            f"PDB ID: {pdb_id} | Exhaustiveness: {exhaustiveness}\n"
            f"Output: {run_dir}",
            border_style="blue",
        )
    )

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task1 = progress.add_task(
                description="Baixando e separando PDB...", total=1
            )
            pdb_file = preparation.download_pdb(pdb_id, raw_dir)
            rec_pdb, lig_pdb = preparation.split_receptor_ligand(
                pdb_file, processed_dir
            )
            progress.update(task1, completed=1)

            task2 = progress.add_task(
                description="Preparando receptor e ligante (PDBQT)...", total=1
            )
            rec_pdbqt = processed_dir / "receptor.pdbqt"
            lig_pdbqt = processed_dir / "ligand.pdbqt"
            preparation.prepare_receptor(rec_pdb, rec_pdbqt)
            preparation.prepare_ligand(lig_pdb, lig_pdbqt)
            progress.update(task2, completed=1)

            task3 = progress.add_task(description="Calculando grid box...", total=1)
            box_params = box_utils.calculate_centroid(lig_pdb)
            progress.update(task3, completed=1)

            task4 = progress.add_task(description="Rodando AutoDock Vina...", total=1)
            docked_out = results_dir / "docked.pdbqt"
            vina_log = results_dir / "vina_log.txt"
            vina_runner.run_vina(
                VINA_BIN,
                rec_pdbqt,
                lig_pdbqt,
                box_params,
                docked_out,
                vina_log,
                exhaustiveness,
            )
            progress.update(task4, completed=1)

            task5 = progress.add_task(description="Analisando resultados...", total=1)
            score = analysis.extract_vina_score(vina_log)
            rmsd, error = analysis.analyze_results(docked_out, lig_pdb, results_dir)
            progress.update(task5, completed=1)

            # Executa fluxo do PLIP imediatamente após a análise inicial
            task_plip = progress.add_task(
                description="Executando PLIP (Docker)...", total=1
            )
            complex_pdb = results_dir / "complex.pdb"
            analysis.generate_complex_pdb(
                rec_pdb, results_dir / "docked_poses.sdf", complex_pdb
            )

            plip_ok, plip_msg = analysis.run_plip_docker(complex_pdb, results_dir)
            if not plip_ok:
                raise RuntimeError(plip_msg)

            interactions = analysis.parse_plip_xml(results_dir / "complex_report.xml")

            # Triagem ADMET
            try:
                admet = pharmacokinetics.calculate_admet_descriptors(results_dir / "docked_poses.sdf")
            except Exception as admet_err:
                admet = {"error": str(admet_err), "pass_filters": False}

            interactions["pharmacokinetics"] = admet

            # Salva o arquivo JSON consolidado
            with open(results_dir / "interactions.json", "w") as f:
                json.dump(interactions, f, indent=4)

            progress.update(task_plip, completed=1)

        console.print("\n[bold green]✓ Validação concluída![/bold green]")
        console.print(f"[bold]Energia de Afinidade (Score):[/bold] {score} kcal/mol")

        if rmsd is not None:
            color = "green" if rmsd <= 2.0 else "red"
            veredito = "SUCESSO" if rmsd <= 2.0 else "FRACASSO"
            console.print(
                f"[bold]RMSD (vs Cristal):[/bold] [{color}]{rmsd:.3f} Å[/{color}]"
            )
            console.print(f"[bold]Veredito:[/bold] [{color}]{veredito}[/{color}]")
        else:
            console.print(f"[bold red]Erro na análise:[/bold red] {error}")

        # Renderiza a tabela de contatos estáticos no terminal
        render_interactions_table(interactions)
        render_admet_table(interactions.get("pharmacokinetics", {}))

    except Exception as e:
        console.print(f"\n[bold red]FATAL ERROR:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="screen")
def screen(
    receptor: Path = typer.Option(
        ..., "--receptor", help="Caminho para o receptor preparado (.pdbqt)"
    ),
    ligand: Path = typer.Option(
        ..., "--ligand", help="Caminho para o novo composto preparado (.pdbqt)"
    ),
    cx: float = typer.Option(..., "--cx", help="Coordenada X do centro do sítio ativo"),
    cy: float = typer.Option(..., "--cy", help="Coordenada Y do centro do sítio ativo"),
    cz: float = typer.Option(..., "--cz", help="Coordenada Z do centro do sítio ativo"),
    size: float = typer.Option(22.0, "--size", help="Tamanho da caixa (A)"),
    exhaustiveness: int = typer.Option(16, "--ex", help="Exaustividade do Vina"),
):
    """
    TRIAGEM VIRTUAL (VIRTUAL SCREENING):
    Executa o docking de um novo xenobiótico em um receptor já preparado em coordenadas específicas.
    """
    global VINA_BIN
    if VINA_BIN is None:
        VINA_BIN = get_vina_bin()

    # Isolamento de output pelo nome do ligante
    ligand_name = ligand.stem
    results_dir = DATA_DIR / "screening" / ligand_name
    results_dir.mkdir(parents=True, exist_ok=True)

    console.print(
        Panel.fit(
            f"[bold cyan]Triagem Virtual[/bold cyan]\n"
            f"Receptor: {receptor.name} | Ligante: {ligand_name}\n"
            f"Box: Center({cx}, {cy}, {cz}) | Size({size})",
            border_style="cyan",
        )
    )

    # Montagem dos parâmetros da caixa
    box_params = {
        "center_x": cx,
        "center_y": cy,
        "center_z": cz,
        "size_x": size,
        "size_y": size,
        "size_z": size,
    }

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            # Execução do Vina
            task1 = progress.add_task(
                description=f"Rodando Vina para {ligand_name}...", total=1
            )
            docked_out = results_dir / f"{ligand_name}_docked.pdbqt"
            vina_log = results_dir / f"{ligand_name}_vina.log"

            vina_runner.run_vina(
                VINA_BIN,
                receptor,
                ligand,
                box_params,
                docked_out,
                vina_log,
                exhaustiveness,
            )
            progress.update(task1, completed=1)

            # Extração de score
            task2 = progress.add_task(description="Extraindo score...", total=1)
            score = analysis.extract_vina_score(vina_log)
            progress.update(task2, completed=1)

            # Exporta para SDF e executa o fluxo do PLIP
            task_plip = progress.add_task(
                description="Executando PLIP (Docker)...", total=1
            )
            sdf_out = results_dir / "docked_poses.sdf"
            import subprocess

            exec_name = get_executable("mk_export")
            subprocess.run([exec_name, str(docked_out), "-s", str(sdf_out)], check=True)

            # Resolve o receptor PDB correspondente
            receptor_pdb = receptor.with_suffix(".pdb")
            if not receptor_pdb.exists():
                pdbs = list(receptor.parent.glob("*.pdb"))
                if pdbs:
                    receptor_pdb = pdbs[0]
                else:
                    raise FileNotFoundError(
                        f"Não foi possível encontrar o arquivo receptor PDB em {receptor.parent}"
                    )

            complex_pdb = results_dir / "complex.pdb"
            analysis.generate_complex_pdb(receptor_pdb, sdf_out, complex_pdb)

            plip_ok, plip_msg = analysis.run_plip_docker(complex_pdb, results_dir)
            if not plip_ok:
                raise RuntimeError(plip_msg)

            interactions = analysis.parse_plip_xml(results_dir / "complex_report.xml")

            # Triagem ADMET
            try:
                admet = pharmacokinetics.calculate_admet_descriptors(sdf_out)
            except Exception as admet_err:
                admet = {"error": str(admet_err), "pass_filters": False}

            interactions["pharmacokinetics"] = admet

            # Salva o arquivo JSON consolidado
            with open(results_dir / "interactions.json", "w") as f:
                json.dump(interactions, f, indent=4)

            progress.update(task_plip, completed=1)

        console.print(
            f"\n[bold green]✓ Triagem concluída para {ligand_name}![/bold green]"
        )
        console.print(
            f"[bold]Energia de Afinidade (Score):[/bold] [yellow]{score}[/yellow] kcal/mol"
        )
        console.print(f"[bold]Resultado salvo em:[/bold] {docked_out}")

        # Renderiza a tabela de contatos estáticos no terminal
        render_interactions_table(interactions)
        render_admet_table(interactions.get("pharmacokinetics", {}))

    except Exception as e:
        console.print(f"\n[bold red]FATAL ERROR during screening:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="interactive")
def interactive():
    """Interface interativa (TUI) para facilitar o uso do pipeline."""
    while True:
        choice = questionary.select(
            "O que você deseja fazer?",
            choices=[
                "1. Validação (Redocking)",
                "2. Download de Ligante (PubChem)",
                "3. Preparação de Ligante (SDF -> PDBQT)",
                "4. Triagem Virtual (Screening)",
                "5. Sair",
            ],
        ).ask()

        if choice == "1. Validação (Redocking)":
            pdb_id = questionary.text(
                "Digite o ID do PDB (ex: 4HG7):", default="4HG7"
            ).ask()
            ex = questionary.text("Exaustividade (ex: 16):", default="16").ask()
            validate(pdb_id=pdb_id, exhaustiveness=int(ex))

        elif choice == "2. Download de Ligante (PubChem)":
            cid = questionary.text("Digite o CID do composto no PubChem:").ask()
            name = questionary.text(
                "Digite o nome do arquivo (ex: desoxicolato.sdf):"
            ).ask()
            if not name.endswith(".sdf"):
                name += ".sdf"
            out_path = DATA_DIR / name
            preparation.download_pubchem_sdf(cid, out_path)
            console.print(f"[bold green]✓ Download concluído:[/bold green] {out_path}")

        elif choice == "3. Preparação de Ligante (SDF -> PDBQT)":
            sdf_file = questionary.path("Caminho para o arquivo SDF:").ask()
            name = Path(sdf_file).stem
            out_pdbqt = questionary.text(
                "Caminho de saída (PDBQT):", default=f"data/{name}.pdbqt"
            ).ask()
            preparation.prepare_ligand_sdf(Path(sdf_file), Path(out_pdbqt))
            console.print(
                f"[bold green]✓ Preparação concluída:[/bold green] {out_pdbqt}"
            )

        elif choice == "4. Triagem Virtual (Screening)":
            receptor = questionary.path("Caminho para o receptor (.pdbqt):").ask()
            ligand = questionary.path("Caminho para o ligante (.pdbqt):").ask()
            cx = questionary.text("Coordenada X:").ask()
            cy = questionary.text("Coordenada Y:").ask()
            cz = questionary.text("Coordenada Z:").ask()
            size = questionary.text("Tamanho da caixa (A):", default="22.0").ask()
            ex = questionary.text("Exaustividade (ex: 16):", default="16").ask()

            screen(
                receptor=Path(receptor),
                ligand=Path(ligand),
                cx=float(cx),
                cy=float(cy),
                cz=float(cz),
                size=float(size),
                exhaustiveness=int(ex),
            )

        elif choice == "5. Sair":
            break


if __name__ == "__main__":
    app()
