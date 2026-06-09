"""
mopac_batch_runner.py

Production-quality batch runner for the SMILES → Open Babel → MOPAC workflow.

Dependencies (Python):
    pandas

External tools (must be available in PATH):
    obabel   (Open Babel CLI)
    mopac    (MOPAC semi-empirical quantum chemistry package)

Typical usage::

    import pandas as pd
    from mopac_batch_runner import run_mopac_batch

    df = pd.DataFrame({
        "compound_id": ["ethanol", "benzene"],
        "smiles":      ["CCO",     "c1ccccc1"],
    })

    results = run_mopac_batch(
        df=df,
        smiles_col="smiles",
        id_col="compound_id",
        keywords="PM7 PRECISE GEO-OK",
        workdir="mopac_jobs",
        overwrite=True,
        timeout=900,
    )
    print(results)
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Module logger
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default Open Babel CLI options.
# Change this list to adjust how 3D structures are generated globally,
# or pass obabel_options= per call to run_mopac_batch / generate_3d_with_obabel.
# ---------------------------------------------------------------------------

OBABEL_DEFAULT_OPTIONS: list[str] = ["--gen3d"]

# Fallback option sets tried in order when the primary obabel call produces a
# crash, NaN coordinates, or an XYZ file with overflowed column values.
#
# Rationale for each fallback:
#   1. "--gen3d dist" — pure distance-geometry algorithm; does NOT look up
#      fragment coordinates from rigid-fragments.txt, which is the root cause
#      of the "all zero coordinates" segfault in Open Babel 3.1.1.
#      ("fast", "rough", and "best" all share the same fragment-lookup code
#      path and crash identically on missing fragments.)
#   2. "--gen3d dist" + UFF minimisation — distance geometry followed by a
#      force-field relaxation to clean up NaN/extreme values.
OBABEL_FALLBACK_OPTIONS: list[list[str]] = [
    ["--gen3d", "dist"],
    ["--gen3d", "dist", "--minimize", "--ff", "UFF", "--steps", "500"],
]

# Coordinates whose absolute value exceeds this threshold (Å) are physically
# impossible and indicate a 3D-generation artefact (XYZ column overflow, etc.).
_MAX_COORD_ABS: float = 1e6

# Minimum allowed distance between any two atoms (Å).  Open Babel occasionally
# places two atoms at identical positions (a 3D-generation bug triggered by
# certain ring fluorine patterns).  Atoms closer than this threshold are
# flagged so the fallback obabel options are tried instead of passing a broken
# geometry to MOPAC.
_MIN_ATOM_DIST: float = 0.1

# ---------------------------------------------------------------------------
# Spin-multiplicity → MOPAC keyword mapping.
# Extend for higher multiplicities; _build_keywords has a UHF/MS fallback.
# Verify these keywords against your specific MOPAC version's manual.
# ---------------------------------------------------------------------------

_SPIN_KEYWORDS: dict[int, str] = {
    2: "DOUBLET",
    3: "TRIPLET",
    4: "QUARTET",
    5: "QUINTET",
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class MoleculeJob:
    """All inputs required to process a single molecule."""

    compound_id: str
    smiles: str
    charge: int
    multiplicity: int
    job_dir: Path
    keywords: str


@dataclass
class JobResult:
    """Outcome of processing one molecule through the full pipeline."""

    compound_id: str
    smiles: str
    status: str = "failed"  # success | failed | skipped | timeout
    mop_file: Optional[Path] = None
    out_file: Optional[Path] = None
    arc_file: Optional[Path] = None
    error_message: Optional[str] = None
    return_code: Optional[int] = None
    duration: Optional[float] = None

    def to_dict(self) -> dict:
        """Serialise to a flat dict suitable for a pandas DataFrame row."""
        return {
            "compound_id": self.compound_id,
            "smiles": self.smiles,
            "status": self.status,
            "mop_file": str(self.mop_file) if self.mop_file else None,
            "out_file": str(self.out_file) if self.out_file else None,
            "arc_file": str(self.arc_file) if self.arc_file else None,
            "error_message": self.error_message,
            "return_code": self.return_code,
            "duration_s": round(self.duration, 3) if self.duration is not None else None,
        }


@dataclass
class CommandResult:
    """Result of a subprocess call."""

    returncode: int
    stdout: str
    stderr: str
    timed_out: bool = False


# ---------------------------------------------------------------------------
# Environment validation
# ---------------------------------------------------------------------------


def validate_environment() -> None:
    """Verify that obabel and mopac executables are accessible from PATH.

    Raises:
        EnvironmentError: with a clear message listing every missing executable.
    """
    missing = [tool for tool in ("obabel", "mopac") if shutil.which(tool) is None]
    if missing:
        raise EnvironmentError(
            f"Required executable(s) not found in PATH: {', '.join(missing)}. "
            "Install Open Babel and/or MOPAC and ensure they appear in your "
            "system PATH, then retry."
        )


# ---------------------------------------------------------------------------
# Identifier sanitisation
# ---------------------------------------------------------------------------


def sanitize_id(raw_id: str) -> str:
    """Return a filesystem-safe version of an arbitrary compound identifier.

    Any character that is not alphanumeric, a hyphen, or an underscore is
    replaced by ``_``.  Consecutive underscores are collapsed and
    leading/trailing underscores are stripped.

    Args:
        raw_id: Original identifier (may contain spaces, dots, brackets, …).

    Returns:
        Sanitised string safe for use as a directory or file stem.
        Falls back to ``"mol"`` if the result would otherwise be empty.
    """
    sanitized = re.sub(r"[^\w\-]", "_", str(raw_id))
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    return sanitized or "mol"


# ---------------------------------------------------------------------------
# Generic subprocess helper
# ---------------------------------------------------------------------------


def _decode_stream(data: "str | bytes | None") -> str:
    """Safely decode subprocess stdout/stderr to a plain string."""
    if data is None:
        return ""
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return data


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    timeout: Optional[int] = None,
    log_file: Optional[Path] = None,
) -> CommandResult:
    """Execute an external command and capture its output.

    Does **not** use ``shell=True``.  Pass arguments as a list so that paths
    containing spaces are handled correctly on Windows.

    Args:
        cmd:      Command and arguments as a list.
        cwd:      Working directory for the child process.
        timeout:  Seconds before the process is killed; ``None`` = no limit.
        log_file: When given, a summary + full stdout/stderr is written here.

    Returns:
        :class:`CommandResult` with ``returncode``, ``stdout``, ``stderr``,
        and ``timed_out``.
    """
    logger.debug("run_command: %s  (cwd=%s)", " ".join(cmd), cwd)

    timed_out = False
    stdout = ""
    stderr = ""
    returncode = -1

    try:
        proc = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        returncode = proc.returncode
        stdout = proc.stdout or ""
        stderr = proc.stderr or ""

    except subprocess.TimeoutExpired as exc:
        timed_out = True
        returncode = -1
        # With text=True the partial output may be str or None; bytes if text
        # mode was bypassed internally on some platforms — handle all cases.
        stdout = _decode_stream(exc.stdout)
        stderr = _decode_stream(exc.stderr)

    if log_file is not None:
        try:
            with log_file.open("w", encoding="utf-8") as fh:
                fh.write(f"command    : {' '.join(cmd)}\n")
                fh.write(f"returncode : {returncode}\n")
                fh.write(f"timed_out  : {timed_out}\n")
                fh.write("\n--- stdout ---\n")
                fh.write(stdout)
                fh.write("\n--- stderr ---\n")
                fh.write(stderr)
        except OSError as write_err:
            logger.warning("Could not write log %s: %s", log_file, write_err)

    return CommandResult(
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        timed_out=timed_out,
    )


# ---------------------------------------------------------------------------
# SMILES file writer
# ---------------------------------------------------------------------------


def write_smiles_file(smiles: str, compound_id: str, path: Path) -> None:
    """Write a SMILES string to a ``.smi`` file for Open Babel.

    Uses the tab-separated ``<SMILES>\\t<title>`` convention so that Open Babel
    can read the molecule name from the second field.

    Args:
        smiles:      SMILES string.
        compound_id: Molecule title used by Open Babel as the molecule name.
        path:        Destination file path (created or overwritten).
    """
    path.write_text(f"{smiles}\t{compound_id}\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# 3D structure generation via Open Babel
# ---------------------------------------------------------------------------


def generate_3d_with_obabel(
    smi_file: Path,
    xyz_file: Path,
    obabel_options: Optional[list[str]] = None,
    timeout: int = 120,
) -> CommandResult:
    """Convert a ``.smi`` file to a 3D ``.xyz`` file using the obabel CLI.

    The ``-ismi`` input format flag and ``-oxyz -O <xyz_file>`` output flags
    are always included.  Additional options (e.g. force-field minimisation)
    can be injected via *obabel_options*; by default only ``--gen3d`` is used.

    An ``obabel.log`` file is written next to *smi_file*.

    Example equivalent shell command::

        obabel -ismi input.smi -oxyz -O structure.xyz --gen3d

    Args:
        smi_file:       Input ``.smi`` file path.
        xyz_file:       Desired output ``.xyz`` file path.
        obabel_options: Extra obabel flags appended after the core flags.
                        Defaults to ``OBABEL_DEFAULT_OPTIONS`` (``["--gen3d"]``).
        timeout:        Seconds before the obabel process is killed.

    Returns:
        :class:`CommandResult` from the obabel call.
    """
    if obabel_options is None:
        obabel_options = list(OBABEL_DEFAULT_OPTIONS)

    cmd: list[str] = [
        "obabel",
        "-ismi",
        str(smi_file),
        "-oxyz",
        "-O",
        str(xyz_file),
        *obabel_options,
    ]

    return run_command(
        cmd,
        cwd=smi_file.parent,
        timeout=timeout,
        log_file=smi_file.parent / "obabel.log",
    )


# ---------------------------------------------------------------------------
# XYZ parser
# ---------------------------------------------------------------------------


@dataclass
class XYZAtom:
    """Cartesian coordinates for a single atom."""

    element: str
    x: float
    y: float
    z: float


def parse_xyz(xyz_file: Path) -> list[XYZAtom]:
    """Parse a standard XYZ file and return a list of atoms.

    Expected format::

        <n_atoms>
        <comment line (ignored)>
        <Element>  <x>  <y>  <z>
        ...

    Args:
        xyz_file: Path to the ``.xyz`` file to read.

    Returns:
        List of :class:`XYZAtom` objects in file order.

    Raises:
        ValueError: If the file is too short, the first line is not an integer,
                    the declared atom count is zero, fewer coordinate lines than
                    declared atoms are present, or any coordinate cannot be
                    parsed as a float.
    """
    text = xyz_file.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    if len(lines) < 2:
        raise ValueError(
            f"{xyz_file.name} is too short to be a valid XYZ file "
            f"({len(lines)} line(s) found)."
        )

    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(
            f"First line of {xyz_file.name} must be an integer atom count; "
            f"got: {lines[0]!r}"
        )

    if n_atoms <= 0:
        raise ValueError(
            f"{xyz_file.name} declares {n_atoms} atoms — nothing to process."
        )

    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) < n_atoms:
        raise ValueError(
            f"{xyz_file.name} declares {n_atoms} atoms but only "
            f"{len(atom_lines)} coordinate line(s) were found."
        )

    atoms: list[XYZAtom] = []
    for file_line_no, line in enumerate(atom_lines, start=3):  # 1-based
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(
                f"{xyz_file.name} line {file_line_no}: expected "
                f"'<Element> <x> <y> <z>', got: {line!r}"
            )
        try:
            atoms.append(
                XYZAtom(
                    element=parts[0],
                    x=float(parts[1]),
                    y=float(parts[2]),
                    z=float(parts[3]),
                )
            )
        except ValueError:
            raise ValueError(
                f"{xyz_file.name} line {file_line_no}: cannot parse "
                f"coordinates: {line!r}"
            )

    return atoms


def _validate_atom_coords(atoms: list[XYZAtom]) -> Optional[str]:
    """Return an error string if any coordinate is physically unreasonable.

    Detects three classes of Open Babel 3D-generation artefacts:

    * **Column overflow** — when a float is too wide for the XYZ fixed-width
      field, adjacent numbers merge into one token.  After splitting on
      whitespace the merged token parses to an astronomically large float.
    * **NaN / Inf** — non-finite values produced during force-field building.
    * **Duplicate positions** — two atoms placed at the same point in space
      (distance < ``_MIN_ATOM_DIST``).  Occurs with certain fluorinated-ring
      SMILES when the fragment-database lookup produces overlapping coordinates.
      MOPAC will report "ATOMS N AND M ARE SEPARATED BY 0.0000 ANGSTROMS" and
      abort without producing a ``.arc`` file, even when ``GEO-OK`` is set.

    Args:
        atoms: Parsed atom list from :func:`parse_xyz`.

    Returns:
        A human-readable error string, or ``None`` if all coordinates look
        reasonable.
    """
    import math

    for atom in atoms:
        for label, v in (("x", atom.x), ("y", atom.y), ("z", atom.z)):
            if not math.isfinite(v):
                return (
                    f"Atom {atom.element} has non-finite {label}-coordinate: {v}"
                )
            if abs(v) > _MAX_COORD_ABS:
                return (
                    f"Atom {atom.element} has extreme {label}-coordinate "
                    f"({v:.6g} Å) — likely an XYZ column-overflow artefact."
                )

    min_dist_sq = _MIN_ATOM_DIST ** 2
    for i, a in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            b = atoms[j]
            d_sq = (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            if d_sq < min_dist_sq:
                return (
                    f"Atoms {i + 1} ({a.element}) and {j + 1} ({b.element}) are "
                    f"only {math.sqrt(d_sq):.4f} Å apart — likely a 3D-generation artefact."
                )

    return None


# ---------------------------------------------------------------------------
# MOPAC keyword builder (charge + spin multiplicity)
# ---------------------------------------------------------------------------


def _build_keywords(base_keywords: str, charge: int, multiplicity: int) -> str:
    """Append charge and spin-multiplicity tokens to a MOPAC keyword string.

    MOPAC conventions used here (verify against your MOPAC version's docs):

    * **Charge**:       ``CHARGE=n``  (e.g. ``CHARGE=-1`` for a monoanion).
    * **Multiplicity**: ``DOUBLET`` / ``TRIPLET`` / ``QUARTET`` / ``QUINTET``
      for 2S+1 = 2–5.  For higher multiplicities ``UHF MS=<S>`` is appended,
      where S = (mult − 1) / 2 is the total spin quantum number.

    A neutral singlet (charge=0, multiplicity=1) returns *base_keywords*
    unchanged.

    Args:
        base_keywords: Keyword string provided by the caller.
        charge:        Formal molecular charge (integer).
        multiplicity:  Spin multiplicity 2S+1 ≥ 1.

    Returns:
        Possibly extended keyword string.
    """
    extra: list[str] = []

    if charge != 0:
        extra.append(f"CHARGE={charge}")

    if multiplicity > 1:
        spin_kw = _SPIN_KEYWORDS.get(multiplicity)
        if spin_kw:
            extra.append(spin_kw)
        else:
            # For mult > 5: UHF with explicit total spin S = (mult-1)/2.
            s = (multiplicity - 1) / 2
            ms_str = f"{s:.1f}" if s != int(s) else str(int(s))
            extra.append("UHF")
            extra.append(f"MS={ms_str}")

    if not extra:
        return base_keywords

    return base_keywords.rstrip() + " " + " ".join(extra)


# ---------------------------------------------------------------------------
# MOPAC input file writer
# ---------------------------------------------------------------------------


def write_mop_file(
    mop_file: Path,
    compound_id: str,
    smiles: str,
    atoms: list[XYZAtom],
    keywords: str,
    charge: int = 0,
    multiplicity: int = 1,
) -> None:
    """Write a MOPAC Cartesian-coordinate input (``.mop``) file.

    File structure::

        <keywords line>
        <compound_id>  SMILES: <smiles>
        <blank line>
        <Element>  <x>  1  <y>  1  <z>  1
        ...

    MOPAC reads exactly two header lines (keywords + title) then expects a
    blank separator before the coordinate block.  The SMILES is embedded in
    the title line to preserve traceability without adding a third header line.

    The optimisation flag ``1`` after each coordinate tells MOPAC to optimise
    that degree of freedom.  Change to ``0`` to fix a coordinate.

    Args:
        mop_file:     Destination path (created or overwritten).
        compound_id:  Title placed on the second header line.
        smiles:       Original SMILES string embedded in the description line.
        atoms:        Atom list from :func:`parse_xyz`.
        keywords:     Base MOPAC keyword string (first line of the file).
        charge:       Formal molecular charge (0 = neutral).
        multiplicity: Spin multiplicity (1 = singlet).
    """
    effective_keywords = _build_keywords(keywords, charge, multiplicity)

    lines: list[str] = [
        effective_keywords,
        f"{compound_id}  SMILES: {smiles}",
        "",
    ]

    for atom in atoms:
        # Format: Element(left-padded to 3)  x 1  y 1  z 1
        # Flag value 1 = optimise this coordinate.
        lines.append(
            f"{atom.element:<3s}  {atom.x:>14.8f}  1  {atom.y:>14.8f}  1  {atom.z:>14.8f}  1"
        )

    # MOPAC requires a trailing newline after the last coordinate line.
    mop_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# MOPAC runner
# ---------------------------------------------------------------------------


def run_mopac(
    mop_file: Path,
    timeout: Optional[int] = 600,
) -> CommandResult:
    """Run MOPAC on a single ``.mop`` input file.

    MOPAC is executed in the directory that contains *mop_file* so that all
    output files (``.out``, ``.arc``, ``.den``, …) are written there.

    A ``mopac_process.log`` file is created in the same directory.

    Args:
        mop_file: Path to the ``.mop`` input file.
        timeout:  Seconds before MOPAC is killed; ``None`` = no limit.

    Returns:
        :class:`CommandResult` with ``returncode``, ``stdout``, ``stderr``,
        and ``timed_out``.
    """
    return run_command(
        cmd=["mopac", mop_file.name],
        cwd=mop_file.parent,
        timeout=timeout,
        log_file=mop_file.parent / "mopac_process.log",
    )


# ---------------------------------------------------------------------------
# Single-molecule pipeline
# ---------------------------------------------------------------------------


def process_one_molecule(
    job: MoleculeJob,
    overwrite: bool = False,
    obabel_timeout: int = 120,
    mopac_timeout: Optional[int] = 600,
    obabel_options: Optional[list[str]] = None,
) -> JobResult:
    """Execute the complete SMILES → 3D → MOPAC pipeline for one compound.

    Steps:

    1. Guard against blank SMILES.
    2. Create the job directory.
    3. Skip early if a ``.out`` file already exists and *overwrite* is False.
    4. Write the ``.smi`` input file for Open Babel.
    5. Run Open Babel to produce a 3D ``.xyz`` file.
    6. Parse the ``.xyz`` file.
    7. Write the MOPAC ``.mop`` input file.
    8. Run MOPAC.
    9. Collect output paths and return a :class:`JobResult`.

    Any per-step error is stored in :attr:`JobResult.error_message`; the
    overall batch loop is never interrupted by a single compound's failure.

    Args:
        job:             Molecule inputs encapsulated in :class:`MoleculeJob`.
        overwrite:       When False, skip if a ``.out`` file already exists.
        obabel_timeout:  Timeout (s) for the Open Babel subprocess.
        mopac_timeout:   Timeout (s) for the MOPAC subprocess; ``None`` = no
                         limit.
        obabel_options:  Extra CLI flags forwarded to obabel.

    Returns:
        :class:`JobResult` summarising the outcome of this compound.
    """
    result = JobResult(compound_id=job.compound_id, smiles=job.smiles)

    # -- Guard: blank SMILES ------------------------------------------------
    if not (job.smiles or "").strip():
        result.error_message = "Empty or blank SMILES string."
        return result  # status remains "failed"

    # -- Create job directory -----------------------------------------------
    job.job_dir.mkdir(parents=True, exist_ok=True)

    mop_file = job.job_dir / f"{job.compound_id}.mop"
    out_file = job.job_dir / f"{job.compound_id}.out"
    arc_file = job.job_dir / f"{job.compound_id}.arc"

    # -- Skip existing results when overwrite=False -------------------------
    if not overwrite and out_file.exists():
        result.status = "skipped"
        result.mop_file = mop_file if mop_file.exists() else None
        result.out_file = out_file
        result.arc_file = arc_file if arc_file.exists() else None
        logger.info("[%s] Skipped — .out file already exists.", job.compound_id)
        return result

    t0 = time.monotonic()

    # -- Step 1: write SMILES file ------------------------------------------
    smi_file = job.job_dir / "input.smi"
    try:
        write_smiles_file(job.smiles, job.compound_id, smi_file)
    except OSError as exc:
        result.error_message = f"Cannot write SMILES file: {exc}"
        return result

    # -- Steps 2–3: generate 3D coordinates + parse XYZ (with fallbacks) -----
    # Build the ordered list of option sets to try: caller-supplied first,
    # then each entry from OBABEL_FALLBACK_OPTIONS.
    primary = obabel_options if obabel_options is not None else list(OBABEL_DEFAULT_OPTIONS)
    option_sets = [primary] + OBABEL_FALLBACK_OPTIONS

    xyz_file = job.job_dir / "structure.xyz"
    atoms: list[XYZAtom] = []
    attempt_error = "No obabel attempt was made."

    for attempt, opts in enumerate(option_sets, start=1):
        logger.debug(
            "[%s] obabel attempt %d/%d: %s",
            job.compound_id, attempt, len(option_sets), opts,
        )
        obabel_res = generate_3d_with_obabel(
            smi_file=smi_file,
            xyz_file=xyz_file,
            obabel_options=opts,
            timeout=obabel_timeout,
        )

        if obabel_res.timed_out:
            attempt_error = f"Open Babel timed out after {obabel_timeout}s."
            break  # timeout is not recoverable by retrying

        if obabel_res.returncode != 0:
            attempt_error = (
                f"Open Babel failed (rc={obabel_res.returncode}). "
                f"stderr: {obabel_res.stderr.strip()[:300]}"
            )
            logger.debug("[%s] attempt %d failed: %s", job.compound_id, attempt, attempt_error)
            continue

        if not xyz_file.exists() or xyz_file.stat().st_size == 0:
            attempt_error = (
                f"Open Babel returned rc=0 but {xyz_file.name} is missing or empty. "
                f"stderr: {obabel_res.stderr.strip()[:300]}"
            )
            logger.debug("[%s] attempt %d: empty XYZ.", job.compound_id, attempt)
            continue

        try:
            parsed = parse_xyz(xyz_file)
        except ValueError as exc:
            attempt_error = f"XYZ parse error: {exc}"
            logger.debug("[%s] attempt %d: %s", job.compound_id, attempt, attempt_error)
            continue

        coord_err = _validate_atom_coords(parsed)
        if coord_err:
            attempt_error = f"Coordinate validation failed: {coord_err}"
            logger.debug("[%s] attempt %d: %s", job.compound_id, attempt, attempt_error)
            continue

        # This attempt succeeded.
        atoms = parsed
        if attempt > 1:
            logger.info(
                "[%s] 3D generation succeeded on attempt %d (options: %s).",
                job.compound_id, attempt, opts,
            )
        attempt_error = ""
        break

    if not atoms:
        result.error_message = attempt_error
        return result

    # -- Step 4: write .mop -------------------------------------------------
    try:
        write_mop_file(
            mop_file=mop_file,
            compound_id=job.compound_id,
            smiles=job.smiles,
            atoms=atoms,
            keywords=job.keywords,
            charge=job.charge,
            multiplicity=job.multiplicity,
        )
    except OSError as exc:
        result.error_message = f"Cannot write .mop file: {exc}"
        return result

    result.mop_file = mop_file

    # -- Step 5: run MOPAC --------------------------------------------------
    logger.info("[%s] Running MOPAC ...", job.compound_id)
    mopac_res = run_mopac(mop_file, timeout=mopac_timeout)

    result.return_code = mopac_res.returncode
    result.duration = time.monotonic() - t0

    if mopac_res.timed_out:
        result.status = "timeout"
        result.error_message = f"MOPAC timed out after {mopac_timeout}s."
        logger.warning("[%s] MOPAC timed out.", job.compound_id)
        return result

    # -- Collect output files -----------------------------------------------
    result.out_file = out_file if out_file.exists() else None
    result.arc_file = arc_file if arc_file.exists() else None

    if result.out_file is None:
        result.status = "failed"
        result.error_message = (
            f"MOPAC finished (rc={mopac_res.returncode}) but no .out file was "
            f"produced.  stderr: {mopac_res.stderr.strip()[:400]}"
        )
        return result

    if mopac_res.returncode != 0:
        # .out exists but MOPAC signalled failure — record as failed but keep
        # both output paths so the user can inspect them.
        result.status = "failed"
        result.error_message = (
            f"MOPAC exited with rc={mopac_res.returncode}. "
            f"Inspect {out_file.name} for details."
        )
        logger.warning(
            "[%s] MOPAC rc=%d.", job.compound_id, mopac_res.returncode
        )
        return result

    result.status = "success"
    logger.info("[%s] Done in %.2fs.", job.compound_id, result.duration)
    return result


# ---------------------------------------------------------------------------
# Job-list builder (sequential — must not be parallelised)
# ---------------------------------------------------------------------------


def _build_job_list(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: Optional[str],
    keywords: str,
    workdir: Path,
    charge_col: Optional[str],
    multiplicity_col: Optional[str],
) -> list[MoleculeJob]:
    """Resolve every DataFrame row into a :class:`MoleculeJob`.

    This step must run sequentially so that duplicate sanitised IDs are
    detected and given unique suffixes before any worker touches the
    filesystem.

    Args:
        df:               Input DataFrame (already column-validated).
        smiles_col:       Column holding SMILES strings.
        id_col:           Optional column for compound identifiers.
        keywords:         MOPAC keyword string.
        workdir:          Resolved absolute root directory.
        charge_col:       Optional column for formal charges.
        multiplicity_col: Optional column for spin multiplicities.

    Returns:
        Ordered list of :class:`MoleculeJob` objects, one per DataFrame row.
    """
    jobs: list[MoleculeJob] = []
    seen_ids: dict[str, int] = {}

    for enum_idx, (_, row) in enumerate(df.iterrows(), start=1):

        # -- SMILES ---------------------------------------------------------
        raw_smiles = row[smiles_col]
        smiles = str(raw_smiles) if pd.notna(raw_smiles) else ""

        # -- Charge ---------------------------------------------------------
        charge = 0
        if charge_col:
            raw_charge = row[charge_col]
            if pd.notna(raw_charge):
                try:
                    charge = int(raw_charge)
                except (ValueError, TypeError):
                    logger.warning(
                        "Row %d: cannot convert charge %r to int; using 0.",
                        enum_idx, raw_charge,
                    )

        # -- Multiplicity ---------------------------------------------------
        multiplicity = 1
        if multiplicity_col:
            raw_mult = row[multiplicity_col]
            if pd.notna(raw_mult):
                try:
                    multiplicity = int(raw_mult)
                except (ValueError, TypeError):
                    logger.warning(
                        "Row %d: cannot convert multiplicity %r to int; using 1.",
                        enum_idx, raw_mult,
                    )
            if multiplicity < 1:
                logger.warning(
                    "Row %d: multiplicity %d < 1; clamping to 1.",
                    enum_idx, multiplicity,
                )
                multiplicity = 1

        # -- Compound ID ----------------------------------------------------
        if id_col and pd.notna(row[id_col]):
            raw_id = str(row[id_col])
        else:
            raw_id = f"mol_{enum_idx:06d}"

        base_id = sanitize_id(raw_id)

        if base_id in seen_ids:
            seen_ids[base_id] += 1
            compound_id = f"{base_id}_{seen_ids[base_id]}"
        else:
            seen_ids[base_id] = 0
            compound_id = base_id

        jobs.append(MoleculeJob(
            compound_id=compound_id,
            smiles=smiles,
            charge=charge,
            multiplicity=multiplicity,
            job_dir=workdir / compound_id,
            keywords=keywords,
        ))

    return jobs


# ---------------------------------------------------------------------------
# Public batch entry point
# ---------------------------------------------------------------------------


def run_mopac_batch(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: Optional[str] = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: "str | Path" = "mopac_jobs",
    charge_col: Optional[str] = None,
    multiplicity_col: Optional[str] = None,
    overwrite: bool = False,
    timeout: Optional[int] = 600,
    obabel_options: Optional[list[str]] = None,
    obabel_timeout: int = 120,
    max_workers: int = 1,
) -> pd.DataFrame:
    """Run SMILES → Open Babel → MOPAC for every row of a DataFrame.

    For each row the function:

    1. Resolves a compound identifier from *id_col* or generates a sequential
       one (``mol_000001``, ``mol_000002``, …).
    2. Creates ``<workdir>/<compound_id>/``.
    3. Generates a 3D structure with Open Babel (``--gen3d`` by default,
       with automatic fallback to ``--gen3d dist`` on crashes).
    4. Writes a MOPAC Cartesian-coordinate input file.
    5. Executes MOPAC and waits up to *timeout* seconds.
    6. Records status, output paths, return code, and elapsed time.

    When *max_workers* > 1 the molecule jobs run in a
    :class:`~concurrent.futures.ThreadPoolExecutor`.  Each job writes to its
    own subdirectory so there are no filesystem conflicts.  ID deduplication
    always runs sequentially before the pool is started.

    Individual failures never abort the batch.  At the end the full results
    table is saved to ``<workdir>/batch_results.csv``.

    Args:
        df:               Input DataFrame.  Must contain *smiles_col*.
        smiles_col:       Column name holding SMILES strings.
        id_col:           Optional column to use as compound identifiers.
                          When omitted, sequential IDs ``mol_000001`` …
                          are generated.
        keywords:         MOPAC keyword string written as the first line of
                          every ``.mop`` file.
        workdir:          Root directory for all job subdirectories and for
                          ``batch_results.csv``.
        charge_col:       Optional column containing formal integer charges.
                          When provided, ``CHARGE=n`` is appended to *keywords*
                          for non-zero values.
        multiplicity_col: Optional column containing spin multiplicities (2S+1).
                          When provided, the appropriate MOPAC spin keyword is
                          appended (see :func:`_build_keywords`).
        overwrite:        Re-run even when a ``.out`` file already exists.
        timeout:          Per-compound MOPAC timeout in seconds.
                          ``None`` = no limit.
        obabel_options:   Extra CLI flags forwarded to obabel.
                          Default: ``["--gen3d"]``.
        obabel_timeout:   Per-compound timeout for the Open Babel subprocess.
        max_workers:      Number of parallel worker threads.  ``1`` (default)
                          runs sequentially.  A value of ``4``–``8`` is a
                          reasonable starting point; keep it below the number
                          of physical CPU cores to avoid memory pressure from
                          concurrent MOPAC processes.

    Returns:
        DataFrame with columns: ``compound_id``, ``smiles``, ``status``,
        ``mop_file``, ``out_file``, ``arc_file``, ``error_message``,
        ``return_code``, ``duration_s``.  Row order matches the input
        DataFrame regardless of completion order.

    Raises:
        EnvironmentError: If obabel or mopac is not found in PATH.
        KeyError:         If a referenced column is absent from *df*.
    """
    validate_environment()

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    # -- Column validation --------------------------------------------------
    optional_cols = [c for c in (id_col, charge_col, multiplicity_col) if c]
    missing_cols = [c for c in [smiles_col, *optional_cols] if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Column(s) not found in DataFrame: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    # -- Build ordered job list (sequential — deduplicates IDs) ------------
    jobs = _build_job_list(
        df, smiles_col, id_col, keywords, workdir, charge_col, multiplicity_col
    )
    n = len(jobs)

    # Shared kwargs forwarded to every process_one_molecule call.
    run_kwargs = dict(
        overwrite=overwrite,
        obabel_timeout=obabel_timeout,
        mopac_timeout=timeout,
        obabel_options=obabel_options,
    )

    # -- Execute jobs -------------------------------------------------------
    results: list[Optional[JobResult]] = [None] * n

    if max_workers == 1:
        for idx, job in enumerate(jobs):
            logger.info("[%d/%d] %s", idx + 1, n, job.compound_id)
            try:
                results[idx] = process_one_molecule(job=job, **run_kwargs)
            except Exception as exc:  # pylint: disable=broad-except
                logger.exception("[%s] Unexpected exception: %s", job.compound_id, exc)
                results[idx] = JobResult(
                    compound_id=job.compound_id, smiles=job.smiles,
                    status="failed", error_message=f"Unexpected exception: {exc}",
                )
    else:
        logger.info("Parallel batch: %d jobs across %d workers.", n, max_workers)
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            future_to_idx = {
                pool.submit(process_one_molecule, job=job, **run_kwargs): idx
                for idx, job in enumerate(jobs)
            }
            completed = 0
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                completed += 1
                try:
                    results[idx] = future.result()
                except Exception as exc:  # pylint: disable=broad-except
                    logger.exception(
                        "[%s] Unexpected exception: %s",
                        jobs[idx].compound_id, exc,
                    )
                    results[idx] = JobResult(
                        compound_id=jobs[idx].compound_id,
                        smiles=jobs[idx].smiles,
                        status="failed",
                        error_message=f"Unexpected exception: {exc}",
                    )
                logger.info(
                    "[%d/%d] %s — %s",
                    completed, n, results[idx].compound_id, results[idx].status,
                )

    # -- Build results DataFrame --------------------------------------------
    if results:
        results_df = pd.DataFrame([r.to_dict() for r in results])
    else:
        results_df = pd.DataFrame(
            columns=[
                "compound_id", "smiles", "status",
                "mop_file", "out_file", "arc_file",
                "error_message", "return_code", "duration_s",
            ]
        )

    csv_path = workdir / "batch_results.csv"
    results_df.to_csv(csv_path, index=False)
    logger.info("Batch finished — %d compound(s).  Results → %s", n, csv_path)

    return results_df


# ---------------------------------------------------------------------------
# Script entry point / example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%H:%M:%S",
    )

    sample_df = pd.DataFrame(
        {
            "compound_id": ["ethanol", "benzene", "aspirin"],
            "smiles": [
                "CCO",
                "c1ccccc1",
                "CC(=O)Oc1ccccc1C(=O)O",
            ],
        }
    )

    results = run_mopac_batch(
        df=sample_df,
        smiles_col="smiles",
        id_col="compound_id",
        keywords="PM7 PRECISE GEO-OK",
        workdir="mopac_jobs",
        overwrite=True,
        timeout=900,
    )

    print("\n=== Batch Results ===")
    with pd.option_context("display.max_colwidth", 60):
        print(results.to_string(index=False))
