<p align="center">
  <img src="assets/logo2.png" alt="SMILE and optimize" width="640">
</p>

# SMILE and optimize

**Created by Maciej Gromelski** · v1.3.0 · Licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE) — free for academic, research, and personal use; **commercial use is not permitted**.

Automated batch pipeline: **SMILES → 3D structure (Open Babel) → MOPAC input → MOPAC semi-empirical optimisation → structured results**, available two ways:

- **Desktop GUI** (`gui.py`, "SMILE and optimize") — a point-and-click PySide6 app: load a spreadsheet of SMILES, configure the run, and step through Prepare → Run MOPAC → Export SDF while watching live results and logs. Long runs survive closing the app: see [Workspaces & resuming a session](#workspaces--resuming-a-session).
- **Python API** (`mopac_batch_runner.py`) — call the same pipeline programmatically from a script or Jupyter notebook.

The pipeline itself can be run as a **single combined call** or split into three independent phases:

| Phase | Function | What it does |
|---|---|---|
| 1 | `prepare_batch` | SMILES → 3D coordinates (MOL2 by default) → `.mop` files |
| 2 | `run_mopac_batch` | Scans a directory for `.mop` files and submits each to MOPAC |
| 3 | `export_sdf` | Converts optimised `.arc` files to `.sdf` via Open Babel |
| Combined | `run_pipeline` | Phases 1 + 2 in a single call (SDF export always runs separately) |
| Retry | `retry_failed_batch` | Re-runs only the compounds that failed/timed out in phases 1–2, with extended settings |

---

## Requirements

### Python

- Python 3.9 or later (the GUI has been tested on 3.14)
- Install everything via:

```
pip install -r requirements.txt
```

This installs `pandas` and `numpy` (data handling), `openpyxl` (Excel/CSV I/O), `PySide6` (the desktop GUI), and `rdkit` (optional automatic charge detection — used by `auto_charge=True` in the API and the GUI's "Auto-detect charge" checkbox).

If you only need the library and not the GUI or `auto_charge`, the minimal requirement is just `pandas`:

```
pip install pandas
```

### External tools (must be on system PATH)

| Tool | License | Tested version | Download (exact tested release) |
|---|---|---|---|
| [Open Babel](http://openbabel.org) | GPLv2 | 3.1.0 – 3.1.1 | [openbabel-3-1-1 release page](https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1) (Windows x64 `.exe` installer under Assets) |
| [MOPAC](http://openmopac.net) | Apache License 2.0 | 23.2.5 | [v23.2.5 release page](https://github.com/openmopac/mopac/releases/tag/v23.2.5) (Windows `.exe` installer under Assets) |

Neither tool is bundled with or distributed by this repository — download the installer from the linked release page and run it yourself. This keeps you on the exact version this pipeline was tested against and keeps each tool's own license (and any obligations that come with it, e.g. GPL source-availability for Open Babel) squarely with its own project. See [Third-party software](#third-party-software) below for what that means for this repo's own license.

> **Note:** Some Open Babel builds (including certain 3.1.0 installers) are missing the XYZ output format. The pipeline defaults to **MOL2** (`mol2`) to avoid this issue. See [3D coordinate format](#3d-coordinate-format) for details.

Verify both are accessible before running:

```
obabel --version
mopac --version
```

The GUI also checks this for you at launch (and via its *Check environment* button).

---

## Installation

Clone or copy this repository, then install dependencies:

```
pip install -r requirements.txt
```

If you only want the library (`mopac_batch_runner.py`) and not the GUI, you can instead just copy that single file into your project — no package installation is required beyond `pandas`. The GUI additionally needs `gui.py` and `workspace.py` (workspace save/restore) alongside it.

---

## Desktop GUI

Launch it with:

```
py gui.py
```

### Workflow

1. **Load input data** — under *1. Input Data*, click *Browse…* and pick an `.xlsx` or `.csv` file containing your compounds, then choose which column holds the SMILES strings and (optionally) which holds the compound IDs.
2. **Set the working directory** — where all generated files (`.mop`, `.out`, `.arc`, `.sdf`, per-phase result CSVs) get written. Defaults to `mopac_jobs`.
3. **Configure Settings** — MOPAC keywords, automatic RDKit-based charge detection, overwrite behaviour, per-phase worker counts, *MOPAC threads per job* (leave at *auto* — see [Choosing `max_workers`](#parallel-execution)), Prepare/MOPAC/Export timeouts, retry-specific timeouts and an optional retry MOPAC keyword override (see *Retrying failed jobs* below), and whether to also copy every successfully-exported `.sdf` into one flat "catalogue" folder for easy sharing.
4. **Click SAVE SETTINGS** to lock the Settings section. It prevents changing keywords/workers/timeouts by accident mid-run, and it's also what unlocks the *1. Prepare* button. Click *EDIT SETTINGS* to unlock and make changes again. Locking also writes the current settings into the workspace file (see [Workspaces](#workspaces--resuming-a-session)).
5. **Run the three phases in order:**
   - **1. Prepare (SMILES → .mop)**
   - **2. Run MOPAC**
   - **3. Export SDF** (also copies to the SDF catalogue folder, if enabled)
6. **Watch progress** on the right: a live results table (failed/timeout rows highlighted) refreshed after each phase, and a running log console.
7. **Retry Failed Jobs** — enabled whenever an input file is loaded and any compound is still `failed`/`timeout`/`cancelled`. It re-runs only those compounds (both Open-Babel-stage and MOPAC-stage failures) using the retry timeouts/keywords from Settings, then refreshes the table. Safe to click repeatedly — e.g. loosen the retry timeout further and click again if some jobs still don't converge.

*2. Run MOPAC* and *3. Export SDF* also unlock from what is already in the working directory, so you can point the app at a `workdir` from a previous session — or one produced by the Python API or the notebook — and carry on from there. Only *1. Prepare* and *Retry Failed Jobs* additionally require the input spreadsheet, because both rebuild compound IDs from it.

---

## Workspaces & resuming a session

Long runs outlive the app. A 20,000-compound preparation can take hours, and closing the window used to mean starting over: the input file, the Settings and the phase progress lived only in memory.

The GUI now keeps a **workspace** — a small `workspace.json` written into the working directory, next to the result CSVs and job folders it describes:

```
mopac_jobs/
  workspace.json       # ← input file + columns, Settings, per-phase progress
  prepare_results.csv
  batch_results.csv
  ...
```

### It saves itself

You do not have to remember to save. The workspace is written:

- **when a phase starts** — recorded as `running`, so a crash or a forced exit ninety minutes into a two-hour Prepare is still recoverable;
- **when a phase finishes** — with the per-status counts;
- **when you click SAVE SETTINGS**;
- **when you close the window.**

*File ▸ Save Workspace* (Ctrl+S) forces a save at any time. *File ▸ Save Workspace As…* exports a copy elsewhere — for sharing or archiving a configuration — without moving where autosave writes.

### Reopening one

Point *Working Directory* at a folder that contains a `workspace.json` and the app offers to restore it. *File ▸ Open Workspace* (Ctrl+O) and *File ▸ Recent Workspaces* do the same thing explicitly.

Restoring brings back the input file and its column choices, every Settings value, the lock state, and the results table (reloaded from the per-phase CSVs), then scans the working directory to report exactly what is on disk.

If the run was interrupted or left failures, you are asked how to restore it:

| Choice | Effect |
|---|---|
| **Resume — skip completed** | Turns **Overwrite off**, so compounds that already have their `.mop`/`.out`/`.sdf` are skipped instead of recomputed. This is what you want after an interrupted run. |
| **Restore settings as saved** | Restores the saved `overwrite` value verbatim. Use this when you deliberately want a full recomputation. |

> **Why this matters.** Overwrite defaults to *on* in the GUI. Restoring an interrupted 20,000-compound Prepare and pressing *1. Prepare* again with it on would redo every hour already spent. Resuming turns it off, and the run finishes in minutes.

### If the input file moved or changed

Compound IDs are rebuilt from the whole spreadsheet in row order, so the file matters. The workspace stores a SHA-256 of it and checks on restore:

- **Missing** — you are offered a *Locate file…* dialog. Declining is fine: *2. Run MOPAC* and *3. Export SDF* read the working directory alone and stay usable; only *1. Prepare* and *Retry Failed Jobs* need the spreadsheet.
- **Changed** — you get a warning. If rows were added, removed or reordered, *Retry Failed Jobs* can match compounds to the wrong job folders. Prepare into a fresh working directory instead.

### Menus

| Menu | Contents |
|---|---|
| **File** | Open / Recent / Save / Save As, Load Input File, Close Workspace, Exit |
| **Run** | The four phase actions (Ctrl+1…Ctrl+4) and Stop (Ctrl+.), mirroring the buttons |
| **Tools** | Check Environment, Scan Working Directory, open the working directory or SDF catalogue in Explorer, Clear Log |
| **Help** | Documentation, About |

*Tools ▸ Scan Working Directory* counts the `.mop`/`.out`/`.arc`/`.sdf` files actually present. This is more accurate than the result CSVs after an abrupt exit: those are only rewritten every 200 completed jobs, so they can lag by up to 199 compounds whose files are already on disk.

Deleting `workspace.json` costs you nothing but the GUI-side context — every `.mop`, `.out`, `.arc` and CSV stays where it is, and the later phases still unlock from the directory contents.

---

## Python API

The sections below describe `mopac_batch_runner.py`'s functions directly, for use in your own scripts or Jupyter notebooks. The GUI above wraps exactly these same functions.

### Quick start

#### Three-phase workflow (recommended)

```python
import pandas as pd
from mopac_batch_runner import prepare_batch, run_mopac_batch, export_sdf

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

# Phase 3: export optimised geometries to SDF.
sdf = export_sdf(workdir="mopac_jobs", overwrite=True)

print(sdf[["compound_id", "status", "sdf_file"]])
```

#### Combined single-call workflow

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

### Function reference

#### `prepare_batch`

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
    coord_format: str = "mol2",
    cancel_event: threading.Event | None = None,
) -> pd.DataFrame
```

Converts SMILES strings to 3D coordinates (via Open Babel) and writes a MOPAC `.mop` input file for each compound. No MOPAC calculations are performed.

**Key parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `charge_col` | `str \| None` | `None` | Column with explicit integer charges. `CHARGE=n` is appended to keywords for non-zero values. |
| `auto_charge` | `bool` | `False` | When `True` and `charge_col` is not set, the formal charge is computed from each SMILES using RDKit and added to the keywords as `CHARGE=n` for non-zero values. Requires RDKit (`pip install rdkit`). Ignored when `charge_col` is provided. |
| `coord_format` | `str` | `"mol2"` | Intermediate 3D coordinate format written by Open Babel. One of `"mol2"`, `"sdf"`, or `"xyz"`. Change this if your Open Babel build does not support the default. |
| `cancel_event` | `threading.Event \| None` | `None` | When set from another thread while this is running, any compound not yet started is skipped with `status="cancelled"` and the in-flight obabel process (if any) is killed immediately. See [Cancelling a running batch](#cancelling-a-running-batch). |

**Returns**

A `DataFrame` with one row per input molecule, saved to `<workdir>/prepare_results.csv`:

| Column | Description |
|---|---|
| `compound_id` | Sanitised filesystem identifier. |
| `smiles` | Original SMILES string. |
| `status` | `ready`, `failed`, `skipped`, or `cancelled`. |
| `mop_file` | Absolute path to the `.mop` file, or `None`. |
| `error_message` | Human-readable error description, or `None` on success. |
| `duration_s` | Wall-clock time for the Open Babel step. |

**Raises**

- `EnvironmentError` — if `obabel` is not found in PATH.
- `ImportError` — if `auto_charge=True` and RDKit is not installed.
- `KeyError` — if a referenced column is absent from `df`.

---

#### `run_mopac_batch`

```python
run_mopac_batch(
    workdir: str | Path,
    overwrite: bool = False,
    timeout: int | None = 600,
    max_workers: int = 1,
    cancel_event: threading.Event | None = None,
    mopac_threads: int | None = None,
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
| `cancel_event` | `threading.Event \| None` | `None` | When set from another thread while this is running, any compound not yet started is skipped with `status="cancelled"` and the in-flight MOPAC process(es) (if any) are killed immediately. See [Cancelling a running batch](#cancelling-a-running-batch). |
| `mopac_threads` | `int \| None` | `None` | CPU threads each MOPAC process may use. `None`/`0` = automatic: CPU cores ÷ `max_workers` (at least 1). See [Choosing `max_workers`](#parallel-execution). |

**Returns**

A `DataFrame` with one row per `.mop` file found, saved to `<workdir>/batch_results.csv`:

| Column | Description |
|---|---|
| `compound_id` | Identifier read from the `.mop` file stem. |
| `smiles` | Original SMILES extracted from the `.mop` header. |
| `status` | `success`, `failed`, `skipped`, `timeout`, or `cancelled`. |
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

#### `run_pipeline`

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
    coord_format: str = "mol2",
    cancel_event: threading.Event | None = None,
) -> pd.DataFrame
```

Runs both phases (Open Babel + MOPAC) in a single call. Equivalent to `prepare_batch` followed by `run_mopac_batch`. Use this when you do not need to inspect or edit `.mop` files between phases.

**Parameters** — same as `prepare_batch`, with the addition of `timeout` (MOPAC per-compound limit). `auto_charge`, `coord_format`, and `cancel_event` behave identically to `prepare_batch`.

**Returns** — same schema as `run_mopac_batch`, saved to `<workdir>/batch_results.csv`.

**Raises**

- `EnvironmentError` — if `obabel` or `mopac` is not found in PATH.
- `KeyError` — if a referenced column is absent from `df`.

---

#### `export_sdf`

```python
export_sdf(
    workdir: str | Path,
    overwrite: bool = False,
    timeout: int = 60,
    max_workers: int = 1,
    cancel_event: threading.Event | None = None,
) -> pd.DataFrame
```

Phase 3: converts optimised MOPAC `.arc` files to `.sdf` format using Open Babel. Scans `<workdir>/*/*.arc` — only compounds that completed a successful MOPAC optimisation (i.e. have an `.arc` file) are processed. Compounds whose MOPAC run failed are automatically excluded.

The geometry used is extracted from the `FINAL GEOMETRY OBTAINED` block at the end of each `.arc` file, which contains the fully converged Cartesian coordinates. Open Babel reads this block using the MOPAC Cartesian input format (`-imop`) rather than the `.arc` extension, which produces incorrect output in Open Babel 3.1.1.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `workdir` | `str \| Path` | — | Root directory used for the previous pipeline phases. |
| `overwrite` | `bool` | `False` | Re-convert even when a `.sdf` file already exists. |
| `timeout` | `int` | `60` | Per-compound obabel timeout in seconds. |
| `max_workers` | `int` | `1` | Number of parallel worker threads. |
| `cancel_event` | `threading.Event \| None` | `None` | When set from another thread while this is running, any compound not yet started is skipped with `status="cancelled"` and the in-flight obabel process (if any) is killed immediately. See [Cancelling a running batch](#cancelling-a-running-batch). |

**Returns**

A `DataFrame` with one row per `.arc` file found, saved to `<workdir>/sdf_results.csv`:

| Column | Description |
|---|---|
| `compound_id` | Identifier from the `.arc` file stem. |
| `arc_file` | Absolute path to the source `.arc` file. |
| `sdf_file` | Absolute path to the produced `.sdf` file, or `None`. |
| `status` | `success`, `failed`, `skipped`, or `cancelled`. |
| `error_message` | Human-readable error description, or `None` on success. |
| `duration_s` | Wall-clock time for the Open Babel conversion. |

**Raises**

- `EnvironmentError` — if `obabel` is not found in PATH.
- `FileNotFoundError` — if `workdir` does not exist.

---

#### `retry_failed_batch`

```python
retry_failed_batch(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: str | None = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: str | Path = "mopac_jobs",
    charge_col: str | None = None,
    multiplicity_col: str | None = None,
    auto_charge: bool = False,
    coord_format: str = "mol2",
    obabel_options: list[str] | None = None,
    obabel_timeout: int = 240,
    mopac_timeout: int | None = 1800,
    retry_keywords: str | None = None,
    max_workers: int = 1,
    cancel_event: threading.Event | None = None,
    mopac_threads: int | None = None,
) -> pd.DataFrame
```

Re-runs only the compounds that failed, timed out, or were cancelled in a previous batch, using extended settings, and merges the outcomes back into `prepare_results.csv` / `batch_results.csv`. Compounds that already succeeded are left completely untouched.

It reads `<workdir>/prepare_results.csv` to find compounds where Open Babel never produced a `.mop` file (`status` in `failed`/`cancelled`), and `<workdir>/batch_results.csv` to find compounds where MOPAC failed, timed out, or was cancelled (`status` in `failed`/`timeout`/`cancelled`). Each group is retried independently — Open Babel failures with *obabel_timeout*, MOPAC failures (plus any compound that just got a `.mop` file from the first retry) with *mopac_timeout*. This means a compound stopped via [`cancel_event`](#cancelling-a-running-batch) on an earlier run is retried exactly like a genuine failure or timeout.

**Key parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df`, `smiles_col`, `id_col`, `keywords`, `charge_col`, `multiplicity_col`, `auto_charge`, `coord_format` | — | — | Must match the values used for the **original** `prepare_batch` call. See warning below. |
| `obabel_timeout` | `int` | `240` | Per-compound Open Babel timeout for the retry pass. |
| `mopac_timeout` | `int \| None` | `1800` | Per-compound MOPAC timeout for the retry pass. |
| `retry_keywords` | `str \| None` | `None` | When given, used instead of `keywords` for every retried compound (e.g. to relax convergence criteria). `None` reuses `keywords` unchanged. |
| `cancel_event` | `threading.Event \| None` | `None` | Same as on `prepare_batch`/`run_mopac_batch` — cancels the retry pass itself. See [Cancelling a running batch](#cancelling-a-running-batch). |

> **Warning:** `df` must be the exact same DataFrame (same rows, same order) passed to the original `prepare_batch` call. Compound-ID assignment — auto-generated sequential IDs (`mol_000007`, …) and duplicate-ID suffixing — depends on enumerating the *full* DataFrame in its original order. Passing a filtered or reordered DataFrame will not line up with the existing job directories and their prior results.

**Returns**

The updated `batch_results.csv` DataFrame (same schema as `run_mopac_batch`), covering **all** compounds, not just the retried ones. Also saved back to `<workdir>/batch_results.csv` (and `<workdir>/prepare_results.csv` is updated in place too).

**Raises**

- `EnvironmentError` — if `obabel` or `mopac` is not found in PATH.
- `FileNotFoundError` — if `workdir` does not exist.
- `KeyError` — if a referenced column is absent from `df`.

**Example**

```python
from mopac_batch_runner import prepare_batch, run_mopac_batch, retry_failed_batch

prep = prepare_batch(df=df, smiles_col="smiles", id_col="compound_id", workdir="mopac_jobs")
results = run_mopac_batch(workdir="mopac_jobs", timeout=600)

# Some compounds timed out at 600s -- retry just those with a much longer
# timeout and a slightly relaxed keyword line.
retried = retry_failed_batch(
    df=df,
    smiles_col="smiles",
    id_col="compound_id",
    workdir="mopac_jobs",
    mopac_timeout=3600,
    retry_keywords="PM7 PRECISE GEO-OK GNORM=1.0",
)
print(retried[["compound_id", "status", "duration_s"]])
```

---

### Cancelling a running batch

Every batch-level function (`prepare_batch`, `run_mopac_batch`, `run_pipeline`, `export_sdf`, `retry_failed_batch`) accepts an optional `cancel_event: threading.Event | None`. Setting it **from another thread** while a batch is running kills whatever obabel/MOPAC process is currently in flight within a fraction of a second, and skips every compound that hasn't started yet — instead of waiting for the current job to finish naturally or hit its own `timeout`.

This only matters when you're driving the pipeline from a background thread (e.g. a GUI, as `gui.py` does for its STOP button) so that setting the event from the main thread can actually interrupt a call that's blocking elsewhere:

```python
import threading
from mopac_batch_runner import run_mopac_batch

cancel_event = threading.Event()

def run():
    run_mopac_batch(workdir="mopac_jobs", timeout=1800, cancel_event=cancel_event)

t = threading.Thread(target=run)
t.start()

# ... later, from the main thread (or another thread) ...
cancel_event.set()   # kills the in-flight MOPAC process(es) immediately
t.join()
```

Compounds that were still queued or running when `cancel_event` was set come back with `status="cancelled"` in the returned DataFrame/CSV, distinct from `"failed"` and `"timeout"`. `retry_failed_batch` treats `"cancelled"` exactly like a failure or timeout, so stopping a batch partway through and retrying afterward compose naturally — see [`retry_failed_batch`](#retry_failed_batch) above.

If you call these functions directly on the main thread with no `cancel_event`, nothing changes — this is purely opt-in.

---

### Output directory structure

```
mopac_jobs/
  ethanol/
    input.smi          # SMILES written for Open Babel
    structure.mol2     # 3D coordinates generated by Open Babel (format matches coord_format)
    ethanol.mop        # MOPAC input file
    ethanol.out        # MOPAC output (energies, geometry, etc.)
    ethanol.arc        # MOPAC archive / optimised geometry
    ethanol.sdf        # Optimised structure in SDF format (from export_sdf)
    obabel.log         # obabel stdout + stderr for 3D generation (last attempt)
    obabel_sdf.log     # obabel stdout + stderr for SDF conversion
    mopac_process.log  # Full mopac stdout + stderr
  benzene/
    ...
  aspirin/
    ...
  prepare_results.csv  # Summary table from prepare_batch
  batch_results.csv    # Summary table from run_mopac_batch / run_pipeline
  sdf_results.csv      # Summary table from export_sdf
  workspace.json       # GUI session state (see Workspaces) — not written by the library
```

---

### MOPAC input file format

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

### Common keywords

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

### Handling charged and open-shell molecules

#### Explicit charge column

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

#### Automatic charge detection

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

> **Note:** `auto_charge` (and RDKit in general) requires the ionic center in the SMILES to carry an explicit formal-charge sign (e.g. `[N+]`, `[O-]`) — a SMILES missing that sign will either fail to parse ("Explicit valence ... greater than permitted") or silently parse as a neutral radical instead of the intended ion. If you're working with a library of pre-formed ions (e.g. ionic liquid cations/anions) and see either symptom, check that every charged atom's bracket actually includes its `+`/`-`.

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

### Skipping already-completed jobs

Each phase has its own skip logic, both controlled by `overwrite`.

| Phase | Skipped when | `overwrite=True` forces |
|---|---|---|
| `prepare_batch` | `.mop` file already exists | Re-run Open Babel and rewrite `.mop` |
| `run_mopac_batch` | `.out` file already exists | Re-run MOPAC |
| `run_pipeline` | `.out` file already exists | Re-run both Open Babel and MOPAC |

This makes it safe to resume an interrupted batch. (In the GUI, restoring a workspace in *Resume* mode is exactly this — it sets `overwrite` to `False` for you; see [Workspaces](#workspaces--resuming-a-session).)

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

### 3D coordinate format

Open Babel is used to convert SMILES to a 3D coordinate file before writing the MOPAC input. The format of this intermediate file is controlled by the `coord_format` parameter (and the `OBABEL_3D_FORMAT` constant at the top of the module).

| Format | `coord_format` value | Notes |
|---|---|---|
| Tripos MOL2 | `"mol2"` | **Default.** Widely supported across all Open Babel builds. |
| MDL SDF | `"sdf"` | Also broadly supported; functionally equivalent to MOL2 for this pipeline. |
| XYZ | `"xyz"` | Not available in some Open Babel installers (e.g. certain 3.1.0 builds on Windows). Use only if you have confirmed it works. |

Check which formats your Open Babel build supports:

```
obabel -L formats
```

To use SDF instead of the default MOL2:

```python
prep = prepare_batch(df, smiles_col="smiles", coord_format="sdf")
```

To change the default globally for all calls in a session, edit `OBABEL_3D_FORMAT` at the top of `mopac_batch_runner.py`:

```python
OBABEL_3D_FORMAT: str = "sdf"   # or "xyz" if your build supports it
```

---

### 3D generation fallback

Open Babel 3.1.1 contains a bug where the default `--gen3d` algorithm (`best`) looks up fragment coordinates from `rigid-fragments.txt` and **segfaults** when a fragment is absent from the database. Common triggers include fluoroanilines, fluoropyrimidines, chloroanilines, and certain aliphatic fluorinated rings.

The script handles this automatically with a three-attempt retry chain:

| Attempt | Options used | Notes |
|---|---|---|
| 1 | `--gen3d` *(or caller-supplied `obabel_options`)* | Default; fastest and most accurate when it works. |
| 2 | `--gen3d dist` | Pure distance-geometry algorithm — bypasses the fragment database entirely. Used when attempt 1 crashes or produces extreme coordinates. |
| 3 | `--gen3d dist --minimize --ff UFF --steps 500` | Distance geometry followed by UFF force-field relaxation. Used when attempt 2 produces a valid but poorly optimised geometry. |

Each attempt also validates the resulting coordinates:

- Any atom with `|x|`, `|y|`, or `|z|` > 10⁶ Å (a 3D-generation artefact from Open Babel) triggers the next fallback.
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

### Parallel execution

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

**MOPAC threads vs. workers**

MOPAC is built with OpenMP and, by default, opens **one thread per CPU core in every process it runs**. With `max_workers=4` on a 32-core machine that is 128 busy threads competing for 32 cores; they spend most of their time waiting on each other and a 20-atom PM7 optimisation that should take seconds takes ten minutes.

`run_mopac_batch`, `run_pipeline` and `retry_failed_batch` therefore cap each MOPAC process through its environment (`OMP_NUM_THREADS`, plus the MKL/OpenBLAS equivalents). The cap is controlled by `mopac_threads` (the *MOPAC threads per job* field in the GUI):

| `mopac_threads` | Effect |
|---|---|
| `None` / `0` (default, GUI *auto*) | CPU cores ÷ `max_workers`, at least 1. The workers together use the whole machine without oversubscribing it. |
| `n > 0` | Exactly `n` threads per MOPAC process. Keep `max_workers × n` at or below the core count. |

For typical drug-like molecules (up to a few hundred atoms) MOPAC gains little from multithreading, so **many workers × 1 thread** gives the best throughput — on 32 cores, `max_workers=28` finishes a 35-compound PFAS batch in half a minute. Reserve several threads per job (few workers) only for very large systems where a single optimisation is itself the bottleneck.

Calling `run_mopac` directly (one file, no batch) leaves MOPAC's own default in place unless you pass `threads=`.

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

### Reading results after the batch

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

### Logging

The module uses Python's standard `logging` module under the logger name `mopac_batch_runner`. Enable it in your script:

```python
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
```

Set `level=logging.DEBUG` to see every subprocess command that is executed, including which obabel fallback attempt is being tried.

The GUI attaches its own log handler to this same logger and streams everything into its on-screen log console, so this applies whether you're driving the pipeline from a script or from `gui.py`.

---

### Known limitations

- **Windows only** — tested on Windows 11 with MOPAC 23.2.1 and Open Babel 3.1.0–3.1.1. The code uses `pathlib` and list-based `subprocess` calls so it should be portable, but is untested on Linux/macOS.
- **XYZ format not available in all Open Babel builds** — some installers (including certain 3.1.0 Windows builds) omit the XYZ output format. The pipeline defaults to MOL2 to avoid this. If you see `obabel: cannot write output format!` in `obabel.log`, confirm your format with `obabel -L formats` and pass `coord_format="mol2"` (or `"sdf"`) explicitly.
- **MOPAC keyword line length** — MOPAC has an 80-character limit per keyword line. Long auto-generated keyword strings (charge + high spin multiplicity) may need manual splitting.
- **`obabel.log` shows only the last attempt** — when fallbacks are triggered, earlier attempt logs are overwritten. Inspect the `error_message` column in `prepare_results.csv` to see which attempt failed and why.
- **`run_mopac_batch` discovers `.mop` files alphabetically** — the row order in `batch_results.csv` reflects filesystem sort order, not the original DataFrame order. Join on `compound_id` if you need to align with the input.
- **`retry_failed_batch` requires the original DataFrame, unmodified** — it reconstructs compound IDs by re-running the same ID-assignment logic as `prepare_batch` over the full `df`. If you used auto-generated sequential IDs (no `id_col`) and pass a different or reordered DataFrame on retry, compounds will be matched to the wrong directories. Always keep the original `df` around for the lifetime of a `workdir`. The GUI guards against this by fingerprinting the input file in the workspace and warning on restore if it changed, but the library cannot check for you.
- **No numeric progress** — the GUI progress bar is indeterminate. Per-compound progress is visible only as `[i/n]` lines in the log console, and as counts in the results table once a phase finishes.
- **`cancel_event` kills the whole obabel/MOPAC process, not a sub-step within it** — cancellation terminates the child process tree, so a partially-completed optimisation cannot be resumed from where it stopped. A killed MOPAC job leaves behind whatever partial `.out` file it had written; the next `overwrite=True` run (including `retry_failed_batch`) simply overwrites it, same as an ordinary timeout.
- **A wedged external process leaks two reader threads** — if obabel or MOPAC survives a process-tree kill *and* something still holds its output pipes open, the pipes are abandoned so the batch can continue. This leaks the pair of daemon threads `subprocess` uses to read them, for the remaining life of the process. It is bounded and rare (once per wedged job), and preferable to the alternative, which was hanging the batch.

---

## Changelog

### v1.3.0

- **Workspaces — the GUI can now pick up where it left off.** A long run used to end at the window's close button: the input file, the Settings and the phase progress lived only in memory, so a 20,000-compound preparation that took two hours had to be started over. The GUI now writes a `workspace.json` into the working directory and restores from it. See [Workspaces & resuming a session](#workspaces--resuming-a-session).
  - **Saves itself** at the start of each phase (recorded as `running`, which is what makes a crash mid-run recoverable), at the end of each phase, on SAVE SETTINGS, and on close — plus *File ▸ Save Workspace*. Written atomically, so an interrupted write cannot leave a corrupt file.
  - **Restores** the input file and its column choices, every Settings value, the lock state, and the results table (reloaded from the per-phase CSVs), then scans the working directory to report exactly what is on disk.
  - **Resume mode** turns Overwrite *off* when a restored run was interrupted or left failures, so compounds that already have their `.mop`/`.out`/`.sdf` are skipped rather than recomputed. Overwrite defaults to on, so without this a resumed run would silently redo everything.
  - **Input-file fingerprinting** — a SHA-256 is stored with the workspace. A missing file offers a *Locate file…* dialog; a changed one warns, because compound IDs are rebuilt from the spreadsheet in row order.
- **Application menu bar** — File / Run / Tools / Help, with keyboard shortcuts for every phase (Ctrl+1…Ctrl+4, Ctrl+. to stop). Adds a recent-workspaces list, *Tools ▸ Scan Working Directory*, shortcuts to open the working directory and SDF catalogue in Explorer, and About. The status bar now shows which workspace is loaded and when it was last written.
- **Known limitation retired: GUI phases are no longer strictly sequential.** *2. Run MOPAC* and *3. Export SDF* now unlock from what is in the working directory (`*/*.mop` and `*/*.arc` respectively), matching what the library always allowed. A `workdir` produced by the Python API, by `runner.ipynb`, or by an earlier session is picked up as-is — no workspace file needed. *1. Prepare* and *Retry Failed Jobs* still require the input spreadsheet, since both rebuild compound IDs from it.
- **Fixed: `retry_failed_batch` wrote its result CSVs non-atomically** — unlike every other writer in the pipeline. A crash partway through truncated exactly the files a resumed session reads back to work out what already succeeded. Both now go through the same temp-file-and-rename path as the phase functions.
- **Fixed: a dead `QTimer` leaked per phase run** — the worker poll timer was constructed inside `_start_phase`, leaving an orphaned timer parented to the window after every phase. It is now created once.

### v1.2.0

- **Fixed: one bad molecule could hang an entire batch and discard every result.** A compound whose obabel process survived being killed — because its output pipe handles had been inherited by another process — left the worker thread blocked forever in an unbounded `communicate()`. Since results were only written after the *last* job returned, a 5,590-compound run produced 5,582 `.mop` files and no `prepare_results.csv` at all. Four independent changes now make that impossible:
  - **Bounded process teardown** — killing a child is capped at ~5s and kills the whole process tree first (killing the direct child alone orphans anything it spawned, putting it out of reach). If the pipes are still held after that they are abandoned rather than waited on, so `run_command` always returns.
  - **Per-compound time budget** — new `molecule_timeout` parameter on `prepare_batch` and `prepare_one_molecule` caps the *sum* of all obabel fallback attempts, not just each one individually. Defaults to `obabel_timeout * attempts + 30`.
  - **Incremental result CSVs** — every phase now rewrites its results CSV every 200 completed jobs (atomically, via a temp file), so an interrupted or wedged run keeps everything finished up to that point instead of losing all of it.
  - **Stall detection** — batches now watch the gap *between* completions rather than waiting indefinitely for stragglers. If no job finishes for the per-job timeout plus a 5-minute grace period, the outstanding jobs are recorded as `failed` and the batch returns. Every submitted compound still gets exactly one row in the results.
- **Fixed: closing the GUI window did not end the process.** The batch worker ran inside a `ThreadPoolExecutor`, whose threads are joined by an `atexit` hook — so a single stuck worker kept the app alive after the window was gone and it had to be killed from the terminal. Closing now cancels the running phase and force-exits after a 10-second grace period.
- **Known limitation retired** — the `cancel_event` note about child processes not spawning their own children no longer applies; cancellation and timeouts now kill the full process tree.
- **Stopped bundling the MOPAC/Open Babel installers** — they were briefly committed under `assets/`; removed so the repo doesn't carry (and isn't responsible for satisfying) Open Babel's GPLv2 source-availability obligation for redistributed binaries. The Requirements table now links directly to the exact tested release pages instead. See [Third-party software](#third-party-software).

### v1.1.0

- **Retry failed/timed-out jobs** — new `retry_failed_batch()` function (library) and **Retry Failed Jobs** button (GUI). Re-runs only the compounds that failed, timed out, or were cancelled in a previous batch — with extended timeouts and an optional MOPAC keyword override — and merges the outcomes back into the existing result CSVs. See [`retry_failed_batch`](#retry_failed_batch).
- **Stop a running batch** — new **STOP** button in the GUI, backed by a `cancel_event: threading.Event` parameter on every batch function (`prepare_batch`, `run_mopac_batch`, `run_pipeline`, `export_sdf`, `retry_failed_batch`). Kills the in-flight obabel/MOPAC process immediately (within a fraction of a second) instead of waiting for it to finish or time out. See [Cancelling a running batch](#cancelling-a-running-batch).
- **New `"cancelled"` status** — compounds interrupted by STOP get this status (distinct from `"failed"`/`"timeout"`) and are automatically picked up by a subsequent Retry Failed Jobs run.
- **Prepare (obabel) timeout is now configurable in the GUI** — previously silently fixed at 120s with no UI control.
- **Clearer batch-completion log messages** — the "finished" log line for each phase now reports a per-status breakdown and explicitly notes `(stopped by user)` when STOP was used, instead of a bare compound count that looked identical whether the run completed cleanly or was interrupted.
- **GUI layout** — the three workflow buttons (`1. Prepare`, `2. Run MOPAC`, `3. Export SDF`) plus STOP and the progress bar moved to the top of the results side (right column), above the results table and log console. **Retry Failed Jobs** moved to the left column, under Settings, next to the retry timeout/keyword fields it reads from.
- **Footer provenance** — author, version, copyright year, and contact email are now loaded from [`provenance.json`](provenance.json) at startup instead of being hardcoded in `gui.py`.

### v1.0.0

Initial release: three-phase pipeline (`prepare_batch` → `run_mopac_batch` → `export_sdf`), combined `run_pipeline`, desktop GUI, RDKit-based `auto_charge`, obabel 3D-generation fallback chain, and parallel execution via `max_workers`.

---

## Third-party software

This repository does not include or redistribute Open Babel or MOPAC. Both are external tools you download and install yourself from their official release pages (see [External tools](#external-tools-must-be-on-system-path)), and this pipeline calls them via `subprocess` as separate command-line programs — they are never linked into or bundled with this code.

Each keeps its own, independent license:

| Tool | License | License text |
|---|---|---|
| Open Babel | GNU General Public License v2 (GPLv2) | https://github.com/openbabel/openbabel/blob/master/COPYING |
| MOPAC | Apache License 2.0 | https://github.com/openmopac/mopac/blob/main/LICENSE |

The **[PolyForm Noncommercial License](LICENSE)** below covers only the code in *this* repository (`mopac_batch_runner.py`, `gui.py`, and everything else here). It does not apply to, and does not change the terms of, Open Babel or MOPAC — installing and using those tools is governed entirely by their own licenses linked above, independent of anything in this repo.

---

## License

This project is licensed under the **[PolyForm Noncommercial License 1.0.0](LICENSE)**.

In short: you're free to use, modify, and redistribute this software for **any noncommercial purpose** — academic research, education, personal projects, nonprofit and government use — at no cost. **Commercial use of any kind is not permitted** without a separate license from the author. See the [LICENSE](LICENSE) file for the full legal text.
