# mopac_batch_runner

Automated batch pipeline: **SMILES → 3D structure (Open Babel) → MOPAC input → MOPAC semi-empirical optimisation → structured results**.

The pipeline can be run as a **single combined call** or split into two independent phases:

| Phase | Function | What it does |
|---|---|---|
| 1 | `prepare_batch` | SMILES → 3D coordinates → `.mop` files |
| 2 | `run_mopac_batch` | Scans a directory for `.mop` files and submits each to MOPAC |
| Combined | `run_pipeline` | Equivalent to calling both phases in sequence |

---

## Requirements

### Python

- Python 3.9 or later
- `pandas`

```
pip install pandas
```

### External tools (must be on system PATH)

| Tool | Tested version | Download |
|---|---|---|
| [Open Babel](http://openbabel.org) | 3.1.1 | https://github.com/openbabel/openbabel/releases |
| [MOPAC](http://openmopac.net) | 23.2.1 | https://github.com/openmopac/mopac/releases |

Verify both are accessible before running:

```
obabel --version
mopac --version
```

---

## Installation

Copy `mopac_batch_runner.py` into your project directory. No package installation is required beyond `pandas`.

---

## Quick start

### Two-phase workflow (recommended)

```python
import pandas as pd
from mopac_batch_runner import prepare_batch, run_mopac_batch

df = pd.DataFrame({
    "compound_id": ["ethanol", "benzene", "aspirin"],
    "smiles":      ["CCO",     "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
})

# Phase 1: generate 3D structures and write .mop files.
prep = prepare_batch(
    df=df,
    smiles_col="smiles",
    id_col="compound_id",
    keywords="PM7 PRECISE GEO-OK",
    workdir="mopac_jobs",
    overwrite=True,
)

print(prep[["compound_id", "status", "duration_s"]])

# (Optional) inspect or edit .mop files here before running MOPAC.

# Phase 2: run MOPAC on every .mop file found in the workdir.
results = run_mopac_batch(
    workdir="mopac_jobs",
    overwrite=True,
    timeout=900,
)

print(results[["compound_id", "status", "duration_s"]])
```

### Combined single-call workflow

```python
from mopac_batch_runner import run_pipeline

results = run_pipeline(
    df=df,
    smiles_col="smiles",
    id_col="compound_id",
    keywords="PM7 PRECISE GEO-OK",
    workdir="mopac_jobs",
    overwrite=True,
    timeout=900,
)

print(results[["compound_id", "status", "duration_s"]])
```

Or run the built-in example directly:

```
py mopac_batch_runner.py
```

---

## Function reference

### `prepare_batch`

```python
prepare_batch(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: str | None = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: str | Path = "mopac_jobs",
    charge_col: str | None = None,
    multiplicity_col: str | None = None,
    overwrite: bool = False,
    obabel_options: list[str] | None = None,
    obabel_timeout: int = 120,
    max_workers: int = 1,
    auto_charge: bool = False,
) -> pd.DataFrame
```

Converts SMILES strings to 3D XYZ coordinates (via Open Babel) and writes a MOPAC `.mop` input file for each compound. No MOPAC calculations are performed.

**Key parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `charge_col` | `str \| None` | `None` | Column with explicit integer charges. `CHARGE=n` is appended to keywords for non-zero values. |
| `auto_charge` | `bool` | `False` | When `True` and `charge_col` is not set, the formal charge is computed from each SMILES using RDKit and added to the keywords as `CHARGE=n` for non-zero values. Requires RDKit (`pip install rdkit`). Ignored when `charge_col` is provided. |

**Returns**

A `DataFrame` with one row per input molecule, saved to `<workdir>/prepare_results.csv`:

| Column | Description |
|---|---|
| `compound_id` | Sanitised filesystem identifier. |
| `smiles` | Original SMILES string. |
| `status` | `ready`, `failed`, or `skipped`. |
| `mop_file` | Absolute path to the `.mop` file, or `None`. |
| `error_message` | Human-readable error description, or `None` on success. |
| `duration_s` | Wall-clock time for the Open Babel step. |

**Raises**

- `EnvironmentError` — if `obabel` is not found in PATH.
- `ImportError` — if `auto_charge=True` and RDKit is not installed.
- `KeyError` — if a referenced column is absent from `df`.

---

### `run_mopac_batch`

```python
run_mopac_batch(
    workdir: str | Path,
    overwrite: bool = False,
    timeout: int | None = 600,
    max_workers: int = 1,
) -> pd.DataFrame
```

Scans `<workdir>/*/*.mop` for prepared input files and submits each to MOPAC. The compound identifier and original SMILES are read directly from the `.mop` header — no DataFrame is required.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `workdir` | `str \| Path` | — | Root directory previously used for `prepare_batch`. All `<id>/<id>.mop` files beneath it are discovered automatically. |
| `overwrite` | `bool` | `False` | Re-run even when a `.out` file already exists. |
| `timeout` | `int \| None` | `600` | Per-compound MOPAC timeout in seconds. `None` = no limit. |
| `max_workers` | `int` | `1` | Number of parallel worker threads. |

**Returns**

A `DataFrame` with one row per `.mop` file found, saved to `<workdir>/batch_results.csv`:

| Column | Description |
|---|---|
| `compound_id` | Identifier read from the `.mop` file stem. |
| `smiles` | Original SMILES extracted from the `.mop` header. |
| `status` | `success`, `failed`, `skipped`, or `timeout`. |
| `mop_file` | Absolute path to the `.mop` input file. |
| `out_file` | Absolute path to the `.out` output file, or `None`. |
| `arc_file` | Absolute path to the `.arc` archive file, or `None`. |
| `error_message` | Human-readable error description, or `None` on success. |
| `return_code` | MOPAC process exit code, or `None` if MOPAC never ran. |
| `duration_s` | Wall-clock time for the MOPAC step. |

**Raises**

- `EnvironmentError` — if `mopac` is not found in PATH.
- `FileNotFoundError` — if `workdir` does not exist.

---

### `run_pipeline`

```python
run_pipeline(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: str | None = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: str | Path = "mopac_jobs",
    charge_col: str | None = None,
    multiplicity_col: str | None = None,
    overwrite: bool = False,
    timeout: int | None = 600,
    obabel_options: list[str] | None = None,
    obabel_timeout: int = 120,
    max_workers: int = 1,
    auto_charge: bool = False,
) -> pd.DataFrame
```

Runs both phases (Open Babel + MOPAC) in a single call. Equivalent to `prepare_batch` followed by `run_mopac_batch`. Use this when you do not need to inspect or edit `.mop` files between phases.

**Parameters** — same as `prepare_batch`, with the addition of `timeout` (MOPAC per-compound limit). `auto_charge` behaves identically to `prepare_batch`.

**Returns** — same schema as `run_mopac_batch`, saved to `<workdir>/batch_results.csv`.

**Raises**

- `EnvironmentError` — if `obabel` or `mopac` is not found in PATH.
- `KeyError` — if a referenced column is absent from `df`.

---

## Output directory structure

```
mopac_jobs/
  ethanol/
    input.smi          # SMILES written for Open Babel
    structure.xyz      # 3D coordinates generated by Open Babel
    ethanol.mop        # MOPAC input file
    ethanol.out        # MOPAC output (energies, geometry, etc.)
    ethanol.arc        # MOPAC archive / optimised geometry
    obabel.log         # Full obabel stdout + stderr (last attempt)
    mopac_process.log  # Full mopac stdout + stderr
  benzene/
    ...
  aspirin/
    ...
  prepare_results.csv  # Summary table from prepare_batch
  batch_results.csv    # Summary table from run_mopac_batch / run_pipeline
```

---

## MOPAC input file format

The `.mop` file written by this script follows the format required by MOPAC:

```
PM7 PRECISE GEO-OK
ethanol  SMILES: CCO

C      1.08535000  1     0.04989000  1    -0.04143000  1
O      3.08132000  1    -0.41155000  1     1.21850000  1
...
```

- **Line 1** — MOPAC keywords.
- **Line 2** — Compound title (identifier + SMILES for traceability).
- **Line 3** — Blank (required separator; a third header line here causes MOPAC to silently skip the calculation).
- **Lines 4+** — Atom coordinates with optimisation flag `1` per coordinate.

> The optimisation flag `1` means "optimise this degree of freedom". Change it to `0` in the source to freeze a coordinate.

---

## Common keywords

| Keyword | Effect |
|---|---|
| `PM7` | PM7 Hamiltonian (recommended general-purpose). |
| `PM6` | PM6 Hamiltonian (older, still widely used). |
| `AM1` | AM1 Hamiltonian (fast, less accurate). |
| `PRECISE` | Tightens convergence criteria (recommended for final geometries). |
| `GEO-OK` | Suppresses warnings about unusual bond lengths. Useful for Open Babel geometries before optimisation. |
| `1SCF` | Single-point energy only — no geometry optimisation. |
| `CHARGE=n` | Formal molecular charge (appended automatically when `charge_col` is set). |
| `DOUBLET` / `TRIPLET` | Spin state (appended automatically when `multiplicity_col` is set). |

---

## Handling charged and open-shell molecules

### Explicit charge column

Pass a column containing formal integer charges:

```python
df = pd.DataFrame({
    "id":           ["phenolate",  "oxygen"],
    "smiles":       ["c1ccc(cc1)[O-]", "[O][O]"],
    "charge":       [-1,           0],
    "multiplicity": [1,            3],   # triplet O2
})

prep = prepare_batch(
    df=df,
    smiles_col="smiles",
    id_col="id",
    keywords="PM7 PRECISE",
    charge_col="charge",
    multiplicity_col="multiplicity",
)

results = run_mopac_batch(workdir="mopac_jobs")
```

The generated keywords line for phenolate will be `PM7 PRECISE CHARGE=-1`.  
The generated keywords line for triplet O2 will be `PM7 PRECISE TRIPLET`.

### Automatic charge detection

When charge information is not available in the dataset, set `auto_charge=True` to have the formal charge computed directly from each SMILES using RDKit:

```python
prep = prepare_batch(
    df=df,
    smiles_col="smiles",
    id_col="id",
    keywords="PM7 PRECISE",
    auto_charge=True,   # detects [NH4+], [O-], etc. automatically
)
```

RDKit reads the charge tokens embedded in the SMILES (e.g. `[O-]`, `[NH4+]`, `[Fe+3]`) and appends `CHARGE=n` to the MOPAC keywords for any non-zero value. A log line is emitted for every compound where a non-zero charge is found:

```
Row 7: auto_charge detected charge=+1.
```

`auto_charge` is ignored when `charge_col` is provided — the explicit column always takes precedence. RDKit must be installed (`pip install rdkit`).

Supported multiplicities and their MOPAC keywords:

| Multiplicity (2S+1) | Keyword added |
|---|---|
| 1 | *(nothing)* |
| 2 | `DOUBLET` |
| 3 | `TRIPLET` |
| 4 | `QUARTET` |
| 5 | `QUINTET` |
| ≥ 6 | `UHF MS=<S>` where S = (mult−1)/2 |

---

## Skipping already-completed jobs

Each phase has its own skip logic, both controlled by `overwrite`.

| Phase | Skipped when | `overwrite=True` forces |
|---|---|---|
| `prepare_batch` | `.mop` file already exists | Re-run Open Babel and rewrite `.mop` |
| `run_mopac_batch` | `.out` file already exists | Re-run MOPAC |
| `run_pipeline` | `.out` file already exists | Re-run both Open Babel and MOPAC |

This makes it safe to resume an interrupted batch:

```python
# Resume preparation — skips compounds whose .mop already exists.
prep = prepare_batch(df=df, smiles_col="smiles", overwrite=False)

# Resume MOPAC — skips compounds whose .out already exists.
results = run_mopac_batch(workdir="mopac_jobs", overwrite=False)
```

To force full recomputation:

```python
prep = prepare_batch(df=df, smiles_col="smiles", overwrite=True)
results = run_mopac_batch(workdir="mopac_jobs", overwrite=True)
```

---

## 3D generation fallback

Open Babel 3.1.1 contains a bug where the default `--gen3d` algorithm (`best`) looks up fragment coordinates from `rigid-fragments.txt` and **segfaults** when a fragment is absent from the database. Common triggers include fluoroanilines, fluoropyrimidines, chloroanilines, and certain aliphatic fluorinated rings.

The script handles this automatically with a three-attempt retry chain:

| Attempt | Options used | Notes |
|---|---|---|
| 1 | `--gen3d` *(or caller-supplied `obabel_options`)* | Default; fastest and most accurate when it works. |
| 2 | `--gen3d dist` | Pure distance-geometry algorithm — bypasses the fragment database entirely. Used when attempt 1 crashes or produces extreme coordinates. |
| 3 | `--gen3d dist --minimize --ff UFF --steps 500` | Distance geometry followed by UFF force-field relaxation. Used when attempt 2 produces a valid but poorly optimised geometry. |

Each attempt also validates the resulting coordinates:

- Any atom with `|x|`, `|y|`, or `|z|` > 10⁶ Å (caused by XYZ column-width overflow in Open Babel) triggers the next fallback.
- Any two atoms closer than 0.1 Å (a duplicate-position artefact from the fragment database) also trigger the next fallback. This prevents MOPAC from aborting with `ATOMS N AND M ARE SEPARATED BY 0.0000 ANGSTROMS`.

The log reports which attempt succeeded:

```
[100280-10-2] 3D generation succeeded on attempt 2 (options: ['--gen3d', 'dist']).
```

To override the primary options while still keeping the automatic fallbacks:

```python
prep = prepare_batch(
    df=df,
    smiles_col="smiles",
    obabel_options=["--gen3d", "--minimize", "--ff", "MMFF94", "--steps", "300"],
)
```

To change the fallback chain globally, edit `OBABEL_FALLBACK_OPTIONS` at the top of the module.

---

## Parallel execution

Both `prepare_batch` and `run_mopac_batch` accept a `max_workers` parameter to run jobs simultaneously. Each job writes to its own subdirectory so there are no filesystem conflicts. The returned DataFrame always preserves the input row order (for `prepare_batch`) or file-discovery order (for `run_mopac_batch`) regardless of completion order.

```python
# 4 preparation workers, then 4 MOPAC workers.
prep = prepare_batch(
    df=df,
    smiles_col="smiles",
    id_col="compound_id",
    keywords="PM7 PRECISE GEO-OK",
    workdir="mopac_jobs",
    max_workers=4,
)

results = run_mopac_batch(
    workdir="mopac_jobs",
    max_workers=4,
    timeout=600,
)
```

**Choosing `max_workers`**

- A value of **4–8** is a practical starting point on a modern workstation.
- Keep it at or below the number of **physical CPU cores** — each MOPAC process is CPU-bound and spawning more workers than cores increases wall time rather than reducing it.
- Memory scales linearly: a typical PM7 optimisation uses 50–200 MB; multiply by `max_workers` to estimate peak RAM.
- `max_workers=1` (the default) is identical to sequential behaviour.

**Progress logging** differs slightly between modes.  Sequential mode logs before each job starts:

```
[1/5] ethanol
[ethanol] Prepared in 0.18s.
```

Parallel mode logs after each job completes (order reflects completion, not input order):

```
Parallel prepare: 5 jobs across 3 workers.
[1/5] benzene — ready
[2/5] ethanol — ready
[3/5] aspirin — ready
```

---

## Reading results after the batch

```python
import pandas as pd
import re
from pathlib import Path

# Preparation results (phase 1).
prep = pd.read_csv("mopac_jobs/prepare_results.csv")
print(prep[prep["status"] != "ready"])  # show any failures

# MOPAC results (phase 2).
results = pd.read_csv("mopac_jobs/batch_results.csv")

def heat_of_formation(out_path):
    """Extract final heat of formation in kcal/mol from a MOPAC .out file."""
    if not out_path or not Path(out_path).exists():
        return None
    text = Path(out_path).read_text(encoding="utf-8", errors="replace")
    m = re.search(r"HEAT OF FORMATION\s*=\s*([-+]?\d+\.\d+)", text)
    return float(m.group(1)) if m else None

results["hof_kcal_mol"] = results["out_file"].map(heat_of_formation)
print(results[["compound_id", "status", "hof_kcal_mol", "duration_s"]])
```

The `.arc` file contains the final optimised geometry in MOPAC format and can be converted back to XYZ with Open Babel:

```
obabel aspirin.arc -oxyz -O aspirin_opt.xyz
```

---

## Logging

The module uses Python's standard `logging` module under the logger name `mopac_batch_runner`. Enable it in your script:

```python
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
```

Set `level=logging.DEBUG` to see every subprocess command that is executed, including which obabel fallback attempt is being tried.

---

## Known limitations

- **Windows only** — tested on Windows 11 with MOPAC 23.2.1 and Open Babel 3.1.1. The code uses `pathlib` and list-based `subprocess` calls so it should be portable, but is untested on Linux/macOS.
- **MOPAC keyword line length** — MOPAC has an 80-character limit per keyword line. Long auto-generated keyword strings (charge + high spin multiplicity) may need manual splitting.
- **`obabel.log` shows only the last attempt** — when fallbacks are triggered, earlier attempt logs are overwritten. Inspect the `error_message` column in `prepare_results.csv` to see which attempt failed and why.
- **`run_mopac_batch` discovers `.mop` files alphabetically** — the row order in `batch_results.csv` reflects filesystem sort order, not the original DataFrame order. Join on `compound_id` if you need to align with the input.
