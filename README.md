# Manual do Pipeline de Docking Molecular

Este projeto automatiza fluxos de trabalho de docking molecular utilizando AutoDock Vina, Meeko e RDKit.

O sistema suporta dois fluxos principais: Validação (Redocking) e Triagem Virtual (Screening).

---

## Comandos Disponíveis

O pipeline é operado através do arquivo `src/main.py`.

### 0. interactive (Interface TUI)
Inicia um menu interativo para execução das etapas do pipeline sem necessidade de argumentos de linha de comando.

**Uso:**
```bash
uv run src/main.py interactive
```
Permite baixar ligantes do PubChem, preparar arquivos e executar dockings através de prompts no terminal.

---

### 1. validate (Controle Positivo / Redocking)
Valida a metodologia de docking comparando a pose gerada com a pose cristalográfica original.

**Etapas executadas:**
1. Download do arquivo PDB do RCSB.
2. Separação de receptor e ligante nativo.
3. Preparação dos arquivos em formato PDBQT.
4. Cálculo do Grid Box baseado no centroide do ligante original.
5. Execução do Vina e cálculo do RMSD final.

**Exemplo de uso:**
```bash
python src/main.py validate --pdb 4HG7 --ex 16
```
- `--pdb`: Código PDB de 4 caracteres.
- `--ex`: Exaustividade do Vina (padrão: 16).

---

### 2. screen (Triagem Virtual / Virtual Screening)
Executa o docking de novas moléculas contra um receptor preparado em coordenadas específicas.

**Etapas executadas:**
1. Criação de diretório de saída em `data/screening/{nome_do_ligante}`.
2. Execução do docking utilizando as coordenadas de centro (`cx, cy, cz`) fornecidas.
3. Extração da energia de afinidade (Score).

**Exemplo de uso:**
```bash
python src/main.py screen --receptor data/4HG7/processed/receptor.pdbqt --ligand novo_composto.pdbqt --cx -24.0 --cy 6.5 --cz -14.2 --size 22.0
```

**Parâmetros:**
- `--receptor`: Caminho para o receptor em formato `.pdbqt`.
- `--ligand`: Caminho para o ligante em formato `.pdbqt`.
- `--cx, --cy, --cz`: Coordenadas do centro do sítio ativo.
- `--size`: Tamanho da caixa cúbica em Å (padrão: 22.0).
- `--ex`: Exaustividade do Vina (padrão: 16).

---

## Estrutura de Arquivos

- `bin/`: Executável do AutoDock Vina.
- `data/`: Armazenamento de entradas e resultados.
  - `<PDB_ID>/`: Dados do fluxo `validate`.
  - `screening/<LIGAND_NAME>/`: Dados do fluxo `screen`.
- `src/`: Código-fonte.
  - `docking/preparation.py`: Conversão e preparação de estruturas.
  - `docking/vina_runner.py`: Interface de execução do Vina.
  - `docking/analysis.py`: Cálculo de RMSD e extração de scores.

---

## Instalação e Requisitos

1. **Python 3.13+**
2. **Dependências:** Gerenciadas via `uv`.
   ```bash
   uv sync
   ```
3. **Ambiente:** O executável `vina.exe` deve estar presente na pasta `bin/`.

---

## Interpretação de Resultados

- **Score (kcal/mol):** Indica a afinidade teórica (valores mais baixos indicam maior afinidade).
- **RMSD (Å):** Disponível apenas no comando `validate`. Valores ≤ 2.0 Å indicam precisão na predição da pose experimental.
