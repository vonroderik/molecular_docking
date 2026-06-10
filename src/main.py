import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
from docking import preparation, box_utils, vina_runner, analysis
import questionary

app = typer.Typer(help="Pipeline de Docking Molecular Automatizado")
console = Console()

# Definição de caminhos base
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
VINA_BIN = BASE_DIR / "bin" / "vina.exe"


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

        console.print(
            f"\n[bold green]✓ Triagem concluída para {ligand_name}![/bold green]"
        )
        console.print(
            f"[bold]Energia de Afinidade (Score):[/bold] [yellow]{score}[/yellow] kcal/mol"
        )
        console.print(f"[bold]Resultado salvo em:[/bold] {docked_out}")

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
            pdb_id = questionary.text("Digite o ID do PDB (ex: 4HG7):", default="4HG7").ask()
            ex = questionary.text("Exaustividade (ex: 16):", default="16").ask()
            validate(pdb_id=pdb_id, exhaustiveness=int(ex))

        elif choice == "2. Download de Ligante (PubChem)":
            cid = questionary.text("Digite o CID do composto no PubChem:").ask()
            name = questionary.text("Digite o nome do arquivo (ex: desoxicolato.sdf):").ask()
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
            console.print(f"[bold green]✓ Preparação concluída:[/bold green] {out_pdbqt}")

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
