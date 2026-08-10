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
import threading
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

# Output format used by Open Babel when generating 3D coordinates.
# Supported values: "xyz" (default), "mol2" (Tripos MOL2), "sdf" (MDL SDF/MOL).
# Change this if your Open Babel installation does not support XYZ output.
OBABEL_3D_FORMAT: str = "mol2"

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
# Formal-charge detection via RDKit (used by auto_charge)
# ---------------------------------------------------------------------------


def _compute_formal_charge(smiles: str) -> int:
    """Return the formal charge of *smiles* using RDKit.

    Args:
        smiles: A valid SMILES string.

    Returns:
        Integer formal charge, or 0 if the SMILES cannot be parsed.

    Raises:
        ImportError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem  # lazy import — only required when auto_charge=True
    except ImportError:
        raise ImportError(
            "RDKit is required for auto_charge=True. "
            "Install it with:  pip install rdkit"
        )
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(
            "auto_charge: could not parse SMILES %r — defaulting to charge=0.", smiles
        )
        return 0
    return Chem.GetFormalCharge(mol)


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
    status: str = "failed"  # success | failed | skipped | timeout | cancelled
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
class PrepResult:
    """Outcome of phase 1: SMILES → 3D coordinates → .mop file."""

    compound_id: str
    smiles: str
    status: str = "failed"  # ready | failed | skipped | cancelled
    mop_file: Optional[Path] = None
    error_message: Optional[str] = None
    duration: Optional[float] = None

    def to_dict(self) -> dict:
        """Serialise to a flat dict suitable for a pandas DataFrame row."""
        return {
            "compound_id": self.compound_id,
            "smiles": self.smiles,
            "status": self.status,
            "mop_file": str(self.mop_file) if self.mop_file else None,
            "error_message": self.error_message,
            "duration_s": round(self.duration, 3) if self.duration is not None else None,
        }


@dataclass
class SdfResult:
    """Outcome of phase 3: .arc → .sdf conversion via Open Babel."""

    compound_id: str
    arc_file: Optional[Path] = None
    sdf_file: Optional[Path] = None
    status: str = "failed"  # success | failed | skipped | cancelled
    error_message: Optional[str] = None
    duration: Optional[float] = None

    def to_dict(self) -> dict:
        """Serialise to a flat dict suitable for a pandas DataFrame row."""
        return {
            "compound_id": self.compound_id,
            "arc_file": str(self.arc_file) if self.arc_file else None,
            "sdf_file": str(self.sdf_file) if self.sdf_file else None,
            "status": self.status,
            "error_message": self.error_message,
            "duration_s": round(self.duration, 3) if self.duration is not None else None,
        }


@dataclass
class CommandResult:
    """Result of a subprocess call."""

    returncode: int
    stdout: str
    stderr: str
    timed_out: bool = False
    cancelled: bool = False


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

# How often the wait loop below wakes up to check *cancel_event* and the
# wall-clock *timeout*. Short enough that a STOP click kills the child
# process within a fraction of a second; long enough not to busy-loop.
_POLL_INTERVAL: float = 0.2


def run_command(
    cmd: list[str],
    cwd: Optional[Path] = None,
    timeout: Optional[int] = None,
    log_file: Optional[Path] = None,
    cancel_event: Optional[threading.Event] = None,
) -> CommandResult:
    """Execute an external command and capture its output.

    Does **not** use ``shell=True``.  Pass arguments as a list so that paths
    containing spaces are handled correctly on Windows.

    Runs the child via :class:`subprocess.Popen` and waits on it in short
    ``_POLL_INTERVAL`` slices (via repeated ``Popen.communicate(timeout=...)``
    calls — the pattern the :mod:`subprocess` docs themselves recommend for a
    cancellable wait) instead of a single blocking ``subprocess.run`` call.
    This lets *cancel_event* be honoured promptly: setting it from another
    thread kills the child process within about one poll interval, rather
    than waiting for it to finish on its own or hit *timeout*.

    Args:
        cmd:          Command and arguments as a list.
        cwd:          Working directory for the child process.
        timeout:      Seconds before the process is killed; ``None`` = no limit.
        log_file:     When given, a summary + full stdout/stderr is written here.
        cancel_event: When set (from another thread) while the process is
                      running, it is killed immediately instead of being
                      allowed to run to completion or to *timeout*.

    Returns:
        :class:`CommandResult` with ``returncode``, ``stdout``, ``stderr``,
        ``timed_out``, and ``cancelled``.
    """
    logger.debug("run_command: %s  (cwd=%s)", " ".join(cmd), cwd)

    timed_out = False
    cancelled = False
    stdout = ""
    stderr = ""
    returncode = -1

    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    start = time.monotonic()
    while True:
        try:
            stdout, stderr = proc.communicate(timeout=_POLL_INTERVAL)
            returncode = proc.returncode
            break
        except subprocess.TimeoutExpired:
            if cancel_event is not None and cancel_event.is_set():
                cancelled = True
                proc.kill()
                stdout, stderr = proc.communicate()
                returncode = proc.returncode
                break
            if timeout is not None and (time.monotonic() - start) >= timeout:
                timed_out = True
                proc.kill()
                stdout, stderr = proc.communicate()
                returncode = proc.returncode
                break

    stdout = stdout or ""
    stderr = stderr or ""

    if log_file is not None:
        try:
            with log_file.open("w", encoding="utf-8") as fh:
                fh.write(f"command    : {' '.join(cmd)}\n")
                fh.write(f"returncode : {returncode}\n")
                fh.write(f"timed_out  : {timed_out}\n")
                fh.write(f"cancelled  : {cancelled}\n")
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
        cancelled=cancelled,
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
    coord_file: Path,
    obabel_options: Optional[list[str]] = None,
    timeout: int = 120,
    coord_format: str = "xyz",
    cancel_event: Optional[threading.Event] = None,
) -> CommandResult:
    """Convert a ``.smi`` file to a 3D coordinate file using the obabel CLI.

    The ``-ismi`` input format flag and ``-o<fmt> -O <coord_file>`` output
    flags are always included.  Additional options (e.g. force-field
    minimisation) can be injected via *obabel_options*.

    Supported *coord_format* values: ``"xyz"`` (default), ``"mol2"``,
    ``"sdf"``.  Use ``"mol2"`` or ``"sdf"`` if your Open Babel build lacks
    XYZ output support.

    An ``obabel.log`` file is written next to *smi_file*.

    Example equivalent shell command (xyz)::

        obabel -ismi input.smi -oxyz -O structure.xyz --gen3d

    Args:
        smi_file:       Input ``.smi`` file path.
        coord_file:     Desired output file path (extension should match *coord_format*).
        obabel_options: Extra obabel flags appended after the core flags.
                        Defaults to ``OBABEL_DEFAULT_OPTIONS`` (``["--gen3d"]``).
        timeout:        Seconds before the obabel process is killed.
        coord_format:   Output format passed to obabel (``-o<fmt>``).
                        Must be one of ``"xyz"``, ``"mol2"``, or ``"sdf"``.
        cancel_event:   When set from another thread, the obabel process is
                        killed immediately. See :func:`run_command`.

    Returns:
        :class:`CommandResult` from the obabel call.
    """
    if obabel_options is None:
        obabel_options = list(OBABEL_DEFAULT_OPTIONS)

    cmd: list[str] = [
        "obabel",
        "-ismi",
        str(smi_file),
        f"-o{coord_format}",
        "-O",
        str(coord_file),
        *obabel_options,
    ]

    return run_command(
        cmd,
        cwd=smi_file.parent,
        timeout=timeout,
        log_file=smi_file.parent / "obabel.log",
        cancel_event=cancel_event,
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


def parse_mol2(mol2_file: Path) -> list[XYZAtom]:
    """Parse a Tripos MOL2 file and return a list of atoms.

    Reads the ``@<TRIPOS>ATOM`` section.  The atom type (e.g. ``C.ar``) is
    split on ``.`` to extract the element symbol.

    Args:
        mol2_file: Path to the ``.mol2`` file to read.

    Returns:
        List of :class:`XYZAtom` objects in file order.

    Raises:
        ValueError: If no atoms are found in the ``@<TRIPOS>ATOM`` section.
    """
    text = mol2_file.read_text(encoding="utf-8", errors="replace")
    in_atom = False
    atoms: list[XYZAtom] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "@<TRIPOS>ATOM":
            in_atom = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atom:
            break
        if not in_atom or not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        # cols: id  name  x  y  z  type  ...
        element = parts[5].split(".")[0]
        try:
            atoms.append(XYZAtom(element=element, x=float(parts[2]), y=float(parts[3]), z=float(parts[4])))
        except ValueError as exc:
            raise ValueError(f"{mol2_file.name}: cannot parse atom line: {line!r}") from exc
    if not atoms:
        raise ValueError(f"{mol2_file.name}: no atoms found in @<TRIPOS>ATOM section.")
    return atoms


def parse_sdf(sdf_file: Path) -> list[XYZAtom]:
    """Parse the first molecule in an MDL SDF/MOL (V2000) file.

    The counts line (line 4, 0-indexed line 3) encodes the atom count in the
    first three characters.  The following atom lines each encode x, y, z in
    fixed-width fields (columns 0-9, 10-19, 20-29) and the element symbol in
    columns 31-33.

    Args:
        sdf_file: Path to the ``.sdf`` or ``.mol`` file to read.

    Returns:
        List of :class:`XYZAtom` objects in file order.

    Raises:
        ValueError: If the file is too short, the atom count cannot be parsed,
                    or any atom line cannot be parsed.
    """
    text = sdf_file.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError(f"{sdf_file.name}: too short to be a valid SDF file ({len(lines)} lines).")
    counts_line = lines[3]
    try:
        n_atoms = int(counts_line[:3])
    except ValueError:
        raise ValueError(
            f"{sdf_file.name}: cannot parse atom count from counts line: {counts_line!r}"
        )
    if n_atoms <= 0:
        raise ValueError(f"{sdf_file.name}: declares {n_atoms} atoms — nothing to process.")
    atom_lines = lines[4 : 4 + n_atoms]
    if len(atom_lines) < n_atoms:
        raise ValueError(
            f"{sdf_file.name}: expected {n_atoms} atom lines but found {len(atom_lines)}."
        )
    atoms: list[XYZAtom] = []
    for i, line in enumerate(atom_lines):
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            element = line[31:34].strip()
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"{sdf_file.name} atom line {i + 1}: cannot parse coordinates: {line!r}"
            ) from exc
        if not element:
            raise ValueError(f"{sdf_file.name} atom line {i + 1}: empty element symbol: {line!r}")
        atoms.append(XYZAtom(element=element, x=x, y=y, z=z))
    return atoms


_SUPPORTED_3D_FORMATS: set[str] = {"xyz", "mol2", "sdf"}


def parse_3d_coords(coord_file: Path, fmt: str) -> list[XYZAtom]:
    """Dispatch to the appropriate parser based on *fmt*.

    Args:
        coord_file: Path to the 3D coordinate file.
        fmt:        Format string — one of ``"xyz"``, ``"mol2"``, ``"sdf"``.

    Returns:
        List of :class:`XYZAtom` objects.

    Raises:
        ValueError: If *fmt* is not a supported format.
    """
    if fmt == "xyz":
        return parse_xyz(coord_file)
    if fmt == "mol2":
        return parse_mol2(coord_file)
    if fmt == "sdf":
        return parse_sdf(coord_file)
    raise ValueError(
        f"Unsupported 3D coordinate format: {fmt!r}. "
        f"Choose from: {sorted(_SUPPORTED_3D_FORMATS)}"
    )


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
    cancel_event: Optional[threading.Event] = None,
) -> CommandResult:
    """Run MOPAC on a single ``.mop`` input file.

    MOPAC is executed in the directory that contains *mop_file* so that all
    output files (``.out``, ``.arc``, ``.den``, …) are written there.

    A ``mopac_process.log`` file is created in the same directory.

    Args:
        mop_file:     Path to the ``.mop`` input file.
        timeout:      Seconds before MOPAC is killed; ``None`` = no limit.
        cancel_event: When set from another thread, MOPAC is killed
                      immediately. See :func:`run_command`.

    Returns:
        :class:`CommandResult` with ``returncode``, ``stdout``, ``stderr``,
        ``timed_out``, and ``cancelled``.
    """
    return run_command(
        cmd=["mopac", mop_file.name],
        cwd=mop_file.parent,
        timeout=timeout,
        log_file=mop_file.parent / "mopac_process.log",
        cancel_event=cancel_event,
    )


# ---------------------------------------------------------------------------
# Phase 3: convert optimised .arc files to .sdf via Open Babel
# ---------------------------------------------------------------------------


def _extract_final_geometry(arc_file: Path) -> Optional[str]:
    """Return the MOPAC Cartesian block after 'FINAL GEOMETRY OBTAINED'.

    The tail of a ``.arc`` file is valid MOPAC Cartesian input (same format
    as ``.mop``).  Open Babel reads it correctly when the input format is
    specified as ``-imop``.

    Returns:
        The geometry block as a string, or ``None`` if the marker is absent
        (e.g. the optimisation did not converge).
    """
    text = arc_file.read_text(encoding="utf-8", errors="replace")
    marker = "FINAL GEOMETRY OBTAINED"
    idx = text.find(marker)
    if idx == -1:
        return None
    return text[idx + len(marker):].lstrip("\n")


def _convert_arc_to_sdf(
    arc_file: Path,
    sdf_file: Path,
    timeout: int = 60,
    overwrite: bool = False,
    cancel_event: Optional[threading.Event] = None,
) -> SdfResult:
    """Convert a single MOPAC ``.arc`` file to ``.sdf`` using obabel.

    The optimised geometry is extracted from the ``FINAL GEOMETRY OBTAINED``
    block of the ``.arc`` file and written to a temporary ``.mop`` file.
    Open Babel converts that with ``-imop`` (which reads MOPAC Cartesian
    format correctly) rather than relying on ``.arc`` extension detection,
    which produces incorrect output in Open Babel 3.1.1.

    Args:
        arc_file:     Path to the ``.arc`` file produced by MOPAC.
        sdf_file:     Desired output ``.sdf`` path (same directory).
        timeout:      Seconds before obabel is killed.
        overwrite:    When False, skip if a ``.sdf`` file already exists.
        cancel_event: When already set, this compound is skipped immediately
                      with ``status="cancelled"`` without spawning obabel.
                      When set while obabel is running, the process is
                      killed immediately. See :func:`run_command`.

    Returns:
        :class:`SdfResult` summarising the conversion.
    """
    compound_id = arc_file.stem
    result = SdfResult(compound_id=compound_id, arc_file=arc_file)

    if cancel_event is not None and cancel_event.is_set():
        result.status = "cancelled"
        result.error_message = "Cancelled by user before this compound started."
        return result

    if not overwrite and sdf_file.exists():
        result.status = "skipped"
        result.sdf_file = sdf_file
        logger.info("[%s] Skipped — .sdf file already exists.", compound_id)
        return result

    geom_block = _extract_final_geometry(arc_file)
    if geom_block is None:
        result.error_message = (
            "Could not find 'FINAL GEOMETRY OBTAINED' in .arc file — "
            "the geometry optimisation may not have completed."
        )
        return result

    tmp_mop = arc_file.parent / "_opt_tmp.mop"
    t0 = time.monotonic()
    try:
        tmp_mop.write_text(geom_block, encoding="utf-8")
        res = run_command(
            cmd=["obabel", "-imop", tmp_mop.name, "-osdf", "-O", sdf_file.name],
            cwd=arc_file.parent,
            timeout=timeout,
            log_file=arc_file.parent / "obabel_sdf.log",
            cancel_event=cancel_event,
        )
    finally:
        if tmp_mop.exists():
            tmp_mop.unlink()

    result.duration = time.monotonic() - t0

    if res.cancelled:
        result.status = "cancelled"
        result.error_message = "Cancelled by user."
        return result

    if res.timed_out:
        result.error_message = f"Open Babel timed out after {timeout}s."
        return result

    if not sdf_file.exists() or sdf_file.stat().st_size == 0:
        result.error_message = (
            f"Open Babel finished (rc={res.returncode}) but {sdf_file.name} is "
            f"missing or empty. stderr: {res.stderr.strip()[:300]}"
        )
        return result

    # Fix the molecule title: obabel inherits the temp filename; replace with
    # the compound_id so the SDF record is correctly labelled.
    try:
        content = sdf_file.read_text(encoding="utf-8", errors="replace")
        first_newline = content.index("\n")
        sdf_file.write_text(compound_id + content[first_newline:], encoding="utf-8")
    except (ValueError, OSError):
        pass  # title fix is cosmetic; don't fail the whole result

    result.status = "success"
    result.sdf_file = sdf_file
    logger.info("[%s] SDF exported in %.2fs.", compound_id, result.duration)
    return result


# ---------------------------------------------------------------------------
# Phase 1: SMILES → 3D coordinates → .mop file
# ---------------------------------------------------------------------------


def prepare_one_molecule(
    job: MoleculeJob,
    overwrite: bool = False,
    obabel_timeout: int = 120,
    obabel_options: Optional[list[str]] = None,
    coord_format: str = OBABEL_3D_FORMAT,
    cancel_event: Optional[threading.Event] = None,
) -> PrepResult:
    """Run phase 1 of the pipeline for one compound: SMILES → .mop file.

    Steps:

    1. Guard against blank SMILES.
    2. Create the job directory.
    3. Skip early if a ``.mop`` file already exists and *overwrite* is False.
    4. Write the ``.smi`` input file for Open Babel.
    5. Run Open Babel to produce a 3D coordinate file (with automatic fallbacks).
    6. Parse and validate the coordinate file.
    7. Write the MOPAC ``.mop`` input file.

    Args:
        job:             Molecule inputs encapsulated in :class:`MoleculeJob`.
        overwrite:       When False, skip if a ``.mop`` file already exists.
        obabel_timeout:  Timeout (s) for the Open Babel subprocess.
        obabel_options:  Extra CLI flags forwarded to obabel.
        coord_format:    3D coordinate format produced by Open Babel.
                         One of ``"xyz"`` (default), ``"mol2"``, ``"sdf"``.
                         Override when your Open Babel build lacks XYZ support.
        cancel_event:    When already set, this compound is skipped
                         immediately with ``status="cancelled"`` without
                         spawning obabel. When set while obabel is running,
                         the process is killed immediately and no further
                         fallback attempts are tried.

    Returns:
        :class:`PrepResult` with ``status`` ``"ready"``, ``"failed"``,
        ``"skipped"``, or ``"cancelled"``.
    """
    result = PrepResult(compound_id=job.compound_id, smiles=job.smiles)

    if cancel_event is not None and cancel_event.is_set():
        result.status = "cancelled"
        result.error_message = "Cancelled by user before this compound started."
        return result

    if not (job.smiles or "").strip():
        result.error_message = "Empty or blank SMILES string."
        return result

    job.job_dir.mkdir(parents=True, exist_ok=True)
    mop_file = job.job_dir / f"{job.compound_id}.mop"

    if not overwrite and mop_file.exists():
        result.status = "skipped"
        result.mop_file = mop_file
        logger.info("[%s] Skipped — .mop file already exists.", job.compound_id)
        return result

    t0 = time.monotonic()

    smi_file = job.job_dir / "input.smi"
    try:
        write_smiles_file(job.smiles, job.compound_id, smi_file)
    except OSError as exc:
        result.error_message = f"Cannot write SMILES file: {exc}"
        result.duration = time.monotonic() - t0
        return result

    primary = obabel_options if obabel_options is not None else list(OBABEL_DEFAULT_OPTIONS)
    option_sets = [primary] + OBABEL_FALLBACK_OPTIONS

    coord_file = job.job_dir / f"structure.{coord_format}"
    atoms: list[XYZAtom] = []
    attempt_error = "No obabel attempt was made."

    for attempt, opts in enumerate(option_sets, start=1):
        logger.debug(
            "[%s] obabel attempt %d/%d: %s",
            job.compound_id, attempt, len(option_sets), opts,
        )
        obabel_res = generate_3d_with_obabel(
            smi_file=smi_file,
            coord_file=coord_file,
            obabel_options=opts,
            timeout=obabel_timeout,
            coord_format=coord_format,
            cancel_event=cancel_event,
        )

        if obabel_res.cancelled:
            result.status = "cancelled"
            result.error_message = "Cancelled by user."
            result.duration = time.monotonic() - t0
            return result

        if obabel_res.timed_out:
            attempt_error = f"Open Babel timed out after {obabel_timeout}s."
            break

        if obabel_res.returncode != 0:
            attempt_error = (
                f"Open Babel failed (rc={obabel_res.returncode}). "
                f"stderr: {obabel_res.stderr.strip()[:300]}"
            )
            logger.debug("[%s] attempt %d failed: %s", job.compound_id, attempt, attempt_error)
            continue

        if not coord_file.exists() or coord_file.stat().st_size == 0:
            attempt_error = (
                f"Open Babel returned rc=0 but {coord_file.name} is missing or empty. "
                f"stderr: {obabel_res.stderr.strip()[:300]}"
            )
            logger.debug("[%s] attempt %d: empty coord file.", job.compound_id, attempt)
            continue

        try:
            parsed = parse_3d_coords(coord_file, coord_format)
        except ValueError as exc:
            attempt_error = f"3D coordinate parse error: {exc}"
            logger.debug("[%s] attempt %d: %s", job.compound_id, attempt, attempt_error)
            continue

        coord_err = _validate_atom_coords(parsed)
        if coord_err:
            attempt_error = f"Coordinate validation failed: {coord_err}"
            logger.debug("[%s] attempt %d: %s", job.compound_id, attempt, attempt_error)
            continue

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
        result.duration = time.monotonic() - t0
        return result

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
        result.duration = time.monotonic() - t0
        return result

    result.status = "ready"
    result.mop_file = mop_file
    result.duration = time.monotonic() - t0
    logger.info("[%s] Prepared in %.2fs.", job.compound_id, result.duration)
    return result


# ---------------------------------------------------------------------------
# Phase 2: run MOPAC on a single .mop file
# ---------------------------------------------------------------------------


def _run_mopac_one(
    mop_file: Path,
    timeout: Optional[int] = 600,
    overwrite: bool = False,
    cancel_event: Optional[threading.Event] = None,
) -> JobResult:
    """Run MOPAC on a single prepared ``.mop`` file.

    Reads the compound identifier and original SMILES from the ``.mop`` header
    (line 2 has the form ``<id>  SMILES: <smiles>``).

    Args:
        mop_file:     Absolute path to the ``.mop`` input file.
        timeout:      Seconds before MOPAC is killed; ``None`` = no limit.
        overwrite:    When False, skip if a ``.out`` file already exists.
        cancel_event: When already set, this compound is skipped immediately
                      with ``status="cancelled"`` without spawning MOPAC.
                      When set while MOPAC is running, it is killed
                      immediately.

    Returns:
        :class:`JobResult` summarising the MOPAC run.
    """
    compound_id = mop_file.stem
    smiles = ""
    try:
        lines = mop_file.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(lines) >= 2 and "SMILES:" in lines[1]:
            smiles = lines[1].split("SMILES:", 1)[1].strip()
    except OSError:
        pass

    out_file = mop_file.parent / f"{compound_id}.out"
    arc_file = mop_file.parent / f"{compound_id}.arc"

    result = JobResult(compound_id=compound_id, smiles=smiles, mop_file=mop_file)

    if cancel_event is not None and cancel_event.is_set():
        result.status = "cancelled"
        result.error_message = "Cancelled by user before this compound started."
        return result

    if not overwrite and out_file.exists():
        result.status = "skipped"
        result.out_file = out_file
        result.arc_file = arc_file if arc_file.exists() else None
        logger.info("[%s] Skipped — .out file already exists.", compound_id)
        return result

    t0 = time.monotonic()
    logger.info("[%s] Running MOPAC ...", compound_id)
    mopac_res = run_mopac(mop_file, timeout=timeout, cancel_event=cancel_event)

    result.return_code = mopac_res.returncode
    result.duration = time.monotonic() - t0

    if mopac_res.cancelled:
        result.status = "cancelled"
        result.error_message = "Cancelled by user."
        logger.info("[%s] MOPAC cancelled by user.", compound_id)
        return result

    if mopac_res.timed_out:
        result.status = "timeout"
        result.error_message = f"MOPAC timed out after {timeout}s."
        logger.warning("[%s] MOPAC timed out.", compound_id)
        return result

    result.out_file = out_file if out_file.exists() else None
    result.arc_file = arc_file if arc_file.exists() else None

    if result.out_file is None:
        result.status = "failed"
        result.error_message = (
            f"MOPAC finished (rc={mopac_res.returncode}) but no .out file was produced. "
            f"stderr: {mopac_res.stderr.strip()[:400]}"
        )
        return result

    if mopac_res.returncode != 0:
        result.status = "failed"
        result.error_message = (
            f"MOPAC exited with rc={mopac_res.returncode}. "
            f"Inspect {out_file.name} for details."
        )
        logger.warning("[%s] MOPAC rc=%d.", compound_id, mopac_res.returncode)
        return result

    result.status = "success"
    logger.info("[%s] Done in %.2fs.", compound_id, result.duration)
    return result


# ---------------------------------------------------------------------------
# Combined single-molecule pipeline (used by run_pipeline)
# ---------------------------------------------------------------------------


def process_one_molecule(
    job: MoleculeJob,
    overwrite: bool = False,
    obabel_timeout: int = 120,
    mopac_timeout: Optional[int] = 600,
    obabel_options: Optional[list[str]] = None,
    coord_format: str = OBABEL_3D_FORMAT,
    cancel_event: Optional[threading.Event] = None,
) -> JobResult:
    """Execute the complete SMILES → 3D → MOPAC pipeline for one compound.

    Delegates to :func:`prepare_one_molecule` (phase 1) and
    :func:`_run_mopac_one` (phase 2).  Skip logic for the combined pipeline
    is based on whether a ``.out`` file already exists (i.e. MOPAC has
    already run), not whether a ``.mop`` file exists.

    Args:
        job:             Molecule inputs encapsulated in :class:`MoleculeJob`.
        overwrite:       When False, skip if a ``.out`` file already exists.
        obabel_timeout:  Timeout (s) for the Open Babel subprocess.
        mopac_timeout:   Timeout (s) for the MOPAC subprocess; ``None`` = no
                         limit.
        obabel_options:  Extra CLI flags forwarded to obabel.
        coord_format:    3D coordinate format (``"xyz"``, ``"mol2"``, ``"sdf"``).
        cancel_event:    When already set, this compound is skipped
                         immediately with ``status="cancelled"``. When set
                         while obabel/MOPAC is running, it is killed
                         immediately. See :func:`run_command`.

    Returns:
        :class:`JobResult` summarising the outcome of this compound.
    """
    out_file = job.job_dir / f"{job.compound_id}.out"
    arc_file = job.job_dir / f"{job.compound_id}.arc"
    mop_file = job.job_dir / f"{job.compound_id}.mop"

    if cancel_event is not None and cancel_event.is_set():
        result = JobResult(compound_id=job.compound_id, smiles=job.smiles)
        result.status = "cancelled"
        result.error_message = "Cancelled by user before this compound started."
        return result

    # Combined-pipeline skip: based on .out, not .mop.
    if not overwrite and out_file.exists():
        result = JobResult(compound_id=job.compound_id, smiles=job.smiles)
        result.status = "skipped"
        result.mop_file = mop_file if mop_file.exists() else None
        result.out_file = out_file
        result.arc_file = arc_file if arc_file.exists() else None
        logger.info("[%s] Skipped — .out file already exists.", job.compound_id)
        return result

    # Phase 1: always regenerate the .mop when we are not skipping.
    prep = prepare_one_molecule(
        job, overwrite=True, obabel_timeout=obabel_timeout,
        obabel_options=obabel_options, coord_format=coord_format,
        cancel_event=cancel_event,
    )

    if prep.status in ("failed", "cancelled"):
        result = JobResult(compound_id=job.compound_id, smiles=job.smiles)
        result.status = "failed" if prep.status == "failed" else "cancelled"
        result.error_message = prep.error_message
        result.duration = prep.duration
        return result

    # Phase 2: run MOPAC (always, since we just wrote a fresh .mop).
    result = _run_mopac_one(
        prep.mop_file, timeout=mopac_timeout, overwrite=True, cancel_event=cancel_event,
    )
    # Total duration = obabel time + MOPAC time.
    result.duration = (prep.duration or 0.0) + (result.duration or 0.0)
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
    auto_charge: bool = False,
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
        elif auto_charge and smiles:
            charge = _compute_formal_charge(smiles)
            if charge != 0:
                logger.info("Row %d: auto_charge detected charge=%+d.", enum_idx, charge)

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
# Public batch entry points
# ---------------------------------------------------------------------------


def _log_batch_finished(
    label: str,
    results_df: pd.DataFrame,
    csv_path: Path,
    cancel_event: Optional[threading.Event] = None,
) -> None:
    """Log a completion summary that reports a per-status breakdown and
    explicitly flags a user-initiated stop.

    A plain "finished — N compound(s)" message (the previous behaviour)
    reads as if everything completed normally even when the run was
    interrupted via *cancel_event* partway through — this makes both the
    outcome mix and the interruption visible in the same line.
    """
    n = len(results_df)
    if n and "status" in results_df.columns:
        counts_text = ", ".join(
            f"{status}: {count}" for status, count in results_df["status"].value_counts().items()
        )
    else:
        counts_text = "no rows"
    stopped_note = " (stopped by user)" if cancel_event is not None and cancel_event.is_set() else ""
    logger.info(
        "%s finished%s — %d compound(s): %s.  Results → %s",
        label, stopped_note, n, counts_text, csv_path,
    )


def prepare_batch(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: Optional[str] = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: "str | Path" = "mopac_jobs",
    charge_col: Optional[str] = None,
    multiplicity_col: Optional[str] = None,
    overwrite: bool = False,
    obabel_options: Optional[list[str]] = None,
    obabel_timeout: int = 120,
    max_workers: int = 1,
    auto_charge: bool = False,
    coord_format: str = OBABEL_3D_FORMAT,
    cancel_event: Optional[threading.Event] = None,
) -> pd.DataFrame:
    """Phase 1 batch: SMILES → 3D coordinates → ``.mop`` files.

    For each row the function generates a job directory, converts the SMILES
    to a 3D coordinate file with Open Babel (with automatic fallbacks), and
    writes a MOPAC ``.mop`` input file.  No MOPAC calculations are performed.

    Run :func:`run_mopac_batch` afterwards to submit the prepared jobs to
    MOPAC.

    Args:
        df:               Input DataFrame.  Must contain *smiles_col*.
        smiles_col:       Column name holding SMILES strings.
        id_col:           Optional column to use as compound identifiers.
                          When omitted, sequential IDs ``mol_000001`` …
                          are generated.
        keywords:         MOPAC keyword string written as the first line of
                          every ``.mop`` file.
        workdir:          Root directory for all job subdirectories and for
                          ``prepare_results.csv``.
        charge_col:       Optional column containing formal integer charges.
        multiplicity_col: Optional column containing spin multiplicities (2S+1).
        overwrite:        Re-generate even when a ``.mop`` file already exists.
        obabel_options:   Extra CLI flags forwarded to obabel.
                          Default: ``["--gen3d"]``.
        obabel_timeout:   Per-compound timeout for the Open Babel subprocess.
        max_workers:      Number of parallel worker threads.
        auto_charge:      When True and *charge_col* is not set, the formal
                          charge is computed from each SMILES string using
                          RDKit and appended to the MOPAC keywords as
                          ``CHARGE=n`` for non-zero values.  Requires RDKit
                          (``pip install rdkit``).  Ignored when *charge_col*
                          is provided.
        coord_format:     3D coordinate format produced by Open Babel.
                          One of ``"xyz"`` (default), ``"mol2"``, ``"sdf"``.
                          Use ``"mol2"`` or ``"sdf"`` if your Open Babel build
                          lacks XYZ output support.
        cancel_event:     When set (from another thread) while this is
                          running, any compound not yet started is skipped
                          with ``status="cancelled"`` and the in-flight
                          obabel process (if any) is killed immediately.

    Returns:
        DataFrame with columns: ``compound_id``, ``smiles``, ``status``
        (``ready`` / ``failed`` / ``skipped`` / ``cancelled``), ``mop_file``,
        ``error_message``, ``duration_s``.  Also saved to
        ``<workdir>/prepare_results.csv``.

    Raises:
        EnvironmentError: If ``obabel`` is not found in PATH.
        ImportError:      If ``auto_charge=True`` and RDKit is not installed.
        KeyError:         If a referenced column is absent from *df*.
    """
    if shutil.which("obabel") is None:
        raise EnvironmentError(
            "obabel not found in PATH. Install Open Babel and ensure it is on "
            "your system PATH, then retry."
        )

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    optional_cols = [c for c in (id_col, charge_col, multiplicity_col) if c]
    missing_cols = [c for c in [smiles_col, *optional_cols] if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Column(s) not found in DataFrame: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    jobs = _build_job_list(
        df, smiles_col, id_col, keywords, workdir, charge_col, multiplicity_col,
        auto_charge=auto_charge,
    )
    n = len(jobs)

    prep_kwargs = dict(
        overwrite=overwrite,
        obabel_timeout=obabel_timeout,
        obabel_options=obabel_options,
        coord_format=coord_format,
        cancel_event=cancel_event,
    )

    results: list[Optional[PrepResult]] = [None] * n

    if max_workers == 1:
        for idx, job in enumerate(jobs):
            logger.info("[%d/%d] %s", idx + 1, n, job.compound_id)
            try:
                results[idx] = prepare_one_molecule(job=job, **prep_kwargs)
            except Exception as exc:  # pylint: disable=broad-except
                logger.exception("[%s] Unexpected exception: %s", job.compound_id, exc)
                results[idx] = PrepResult(
                    compound_id=job.compound_id, smiles=job.smiles,
                    status="failed", error_message=f"Unexpected exception: {exc}",
                )
    else:
        logger.info("Parallel prepare: %d jobs across %d workers.", n, max_workers)
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            future_to_idx = {
                pool.submit(prepare_one_molecule, job=job, **prep_kwargs): idx
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
                        "[%s] Unexpected exception: %s", jobs[idx].compound_id, exc,
                    )
                    results[idx] = PrepResult(
                        compound_id=jobs[idx].compound_id,
                        smiles=jobs[idx].smiles,
                        status="failed",
                        error_message=f"Unexpected exception: {exc}",
                    )
                logger.info(
                    "[%d/%d] %s — %s",
                    completed, n, results[idx].compound_id, results[idx].status,
                )

    results_df = pd.DataFrame([r.to_dict() for r in results]) if results else pd.DataFrame(
        columns=["compound_id", "smiles", "status", "mop_file", "error_message", "duration_s"]
    )

    csv_path = workdir / "prepare_results.csv"
    results_df.to_csv(csv_path, index=False)
    _log_batch_finished("Prepare", results_df, csv_path, cancel_event)
    return results_df


def run_mopac_batch(
    workdir: "str | Path",
    overwrite: bool = False,
    timeout: Optional[int] = 600,
    max_workers: int = 1,
    cancel_event: Optional[threading.Event] = None,
) -> pd.DataFrame:
    """Phase 2 batch: run MOPAC on all ``.mop`` files found in *workdir*.

    Scans ``<workdir>/*/*.mop`` for prepared input files (as written by
    :func:`prepare_batch`) and submits each to MOPAC.  The compound
    identifier and original SMILES are read from the ``.mop`` header so no
    DataFrame is required.

    Args:
        workdir:      Root directory that was used for :func:`prepare_batch`.
                      All ``<id>/<id>.mop`` files beneath it are discovered
                      automatically.
        overwrite:    Re-run even when a ``.out`` file already exists.
        timeout:      Per-compound MOPAC timeout in seconds.  ``None`` = no
                      limit.
        max_workers:  Number of parallel worker threads.
        cancel_event: When set (from another thread) while this is running,
                      any compound not yet started is skipped with
                      ``status="cancelled"`` and the in-flight MOPAC
                      process(es) (if any) are killed immediately.

    Returns:
        DataFrame with columns: ``compound_id``, ``smiles``, ``status``
        (``success`` / ``failed`` / ``skipped`` / ``timeout`` /
        ``cancelled``), ``mop_file``, ``out_file``, ``arc_file``,
        ``error_message``, ``return_code``, ``duration_s``.  Also saved to
        ``<workdir>/batch_results.csv``.

    Raises:
        EnvironmentError: If ``mopac`` is not found in PATH.
        FileNotFoundError: If *workdir* does not exist.
    """
    if shutil.which("mopac") is None:
        raise EnvironmentError(
            "mopac not found in PATH. Install MOPAC and ensure it is on your "
            "system PATH, then retry."
        )

    workdir = Path(workdir).resolve()
    if not workdir.exists():
        raise FileNotFoundError(f"workdir does not exist: {workdir}")

    mop_files = sorted(workdir.glob("*/*.mop"))
    n = len(mop_files)
    if n == 0:
        logger.warning("No .mop files found under %s — nothing to run.", workdir)
        return pd.DataFrame(
            columns=[
                "compound_id", "smiles", "status",
                "mop_file", "out_file", "arc_file",
                "error_message", "return_code", "duration_s",
            ]
        )

    logger.info("Found %d .mop file(s) in %s.", n, workdir)

    run_kwargs = dict(timeout=timeout, overwrite=overwrite, cancel_event=cancel_event)
    results: list[Optional[JobResult]] = [None] * n

    if max_workers == 1:
        for idx, mop_file in enumerate(mop_files):
            logger.info("[%d/%d] %s", idx + 1, n, mop_file.stem)
            try:
                results[idx] = _run_mopac_one(mop_file, **run_kwargs)
            except Exception as exc:  # pylint: disable=broad-except
                logger.exception("[%s] Unexpected exception: %s", mop_file.stem, exc)
                results[idx] = JobResult(
                    compound_id=mop_file.stem, smiles="",
                    status="failed", error_message=f"Unexpected exception: {exc}",
                )
    else:
        logger.info("Parallel MOPAC: %d jobs across %d workers.", n, max_workers)
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            future_to_idx = {
                pool.submit(_run_mopac_one, mop_file, **run_kwargs): idx
                for idx, mop_file in enumerate(mop_files)
            }
            completed = 0
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                completed += 1
                try:
                    results[idx] = future.result()
                except Exception as exc:  # pylint: disable=broad-except
                    logger.exception(
                        "[%s] Unexpected exception: %s", mop_files[idx].stem, exc,
                    )
                    results[idx] = JobResult(
                        compound_id=mop_files[idx].stem, smiles="",
                        status="failed", error_message=f"Unexpected exception: {exc}",
                    )
                logger.info(
                    "[%d/%d] %s — %s",
                    completed, n, results[idx].compound_id, results[idx].status,
                )

    results_df = pd.DataFrame([r.to_dict() for r in results]) if results else pd.DataFrame(
        columns=[
            "compound_id", "smiles", "status",
            "mop_file", "out_file", "arc_file",
            "error_message", "return_code", "duration_s",
        ]
    )

    csv_path = workdir / "batch_results.csv"
    results_df.to_csv(csv_path, index=False)
    _log_batch_finished("MOPAC batch", results_df, csv_path, cancel_event)
    return results_df


def run_pipeline(
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
    auto_charge: bool = False,
    coord_format: str = OBABEL_3D_FORMAT,
    cancel_event: Optional[threading.Event] = None,
) -> pd.DataFrame:
    """Combined pipeline: SMILES → Open Babel → MOPAC in a single call.

    Equivalent to calling :func:`prepare_batch` followed by
    :func:`run_mopac_batch`.  Use the two-step approach when you want to
    inspect or edit ``.mop`` files between the geometry generation and the
    MOPAC run.

    Args:
        df:               Input DataFrame.  Must contain *smiles_col*.
        smiles_col:       Column name holding SMILES strings.
        id_col:           Optional column to use as compound identifiers.
        keywords:         MOPAC keyword string (first line of every ``.mop``).
        workdir:          Root directory for all job subdirectories.
        charge_col:       Optional column containing formal integer charges.
        multiplicity_col: Optional column containing spin multiplicities (2S+1).
        overwrite:        Re-run even when a ``.out`` file already exists.
        timeout:          Per-compound MOPAC timeout in seconds.
        obabel_options:   Extra CLI flags forwarded to obabel.
        obabel_timeout:   Per-compound timeout for the Open Babel subprocess.
        max_workers:      Number of parallel worker threads.
        coord_format:     3D coordinate format produced by Open Babel.
                          One of ``"xyz"`` (default), ``"mol2"``, ``"sdf"``.
        cancel_event:     When set (from another thread) while this is
                          running, any compound not yet started is skipped
                          with ``status="cancelled"`` and the in-flight
                          obabel/MOPAC process(es) (if any) are killed
                          immediately.

    Returns:
        DataFrame with columns: ``compound_id``, ``smiles``, ``status``,
        ``mop_file``, ``out_file``, ``arc_file``, ``error_message``,
        ``return_code``, ``duration_s``.  Also saved to
        ``<workdir>/batch_results.csv``.

    Raises:
        EnvironmentError: If obabel or mopac is not found in PATH.
        KeyError:         If a referenced column is absent from *df*.
    """
    validate_environment()

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    optional_cols = [c for c in (id_col, charge_col, multiplicity_col) if c]
    missing_cols = [c for c in [smiles_col, *optional_cols] if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Column(s) not found in DataFrame: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    jobs = _build_job_list(
        df, smiles_col, id_col, keywords, workdir, charge_col, multiplicity_col,
        auto_charge=auto_charge,
    )
    n = len(jobs)

    run_kwargs = dict(
        overwrite=overwrite,
        obabel_timeout=obabel_timeout,
        mopac_timeout=timeout,
        obabel_options=obabel_options,
        coord_format=coord_format,
        cancel_event=cancel_event,
    )

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
                        "[%s] Unexpected exception: %s", jobs[idx].compound_id, exc,
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

    results_df = pd.DataFrame([r.to_dict() for r in results]) if results else pd.DataFrame(
        columns=[
            "compound_id", "smiles", "status",
            "mop_file", "out_file", "arc_file",
            "error_message", "return_code", "duration_s",
        ]
    )

    csv_path = workdir / "batch_results.csv"
    results_df.to_csv(csv_path, index=False)
    _log_batch_finished("Pipeline", results_df, csv_path, cancel_event)
    return results_df


def export_sdf(
    workdir: "str | Path",
    overwrite: bool = False,
    timeout: int = 60,
    max_workers: int = 1,
    cancel_event: Optional[threading.Event] = None,
) -> pd.DataFrame:
    """Phase 3 batch: convert optimised MOPAC ``.arc`` files to ``.sdf``.

    Scans ``<workdir>/*/*.arc`` for files produced by :func:`run_mopac_batch`
    and converts each to an SDF file using Open Babel.  The SDF is written
    next to the ``.arc`` file as ``<compound_id>.sdf``.

    Only ``.arc`` files are processed; compounds whose MOPAC run failed (no
    ``.arc`` produced) are automatically excluded.

    Args:
        workdir:      Root directory used for the previous pipeline phases.
        overwrite:    Re-convert even when a ``.sdf`` file already exists.
        timeout:      Per-compound obabel timeout in seconds.
        max_workers:  Number of parallel worker threads.
        cancel_event: When set (from another thread) while this is running,
                      any compound not yet started is skipped with
                      ``status="cancelled"`` and the in-flight obabel
                      process (if any) is killed immediately.

    Returns:
        DataFrame with columns: ``compound_id``, ``arc_file``, ``sdf_file``,
        ``status`` (``success`` / ``failed`` / ``skipped`` / ``cancelled``),
        ``error_message``, ``duration_s``.  Also saved to
        ``<workdir>/sdf_results.csv``.

    Raises:
        EnvironmentError:  If ``obabel`` is not found in PATH.
        FileNotFoundError: If *workdir* does not exist.
    """
    if shutil.which("obabel") is None:
        raise EnvironmentError(
            "obabel not found in PATH. Install Open Babel and ensure it is on "
            "your system PATH, then retry."
        )

    workdir = Path(workdir).resolve()
    if not workdir.exists():
        raise FileNotFoundError(f"workdir does not exist: {workdir}")

    arc_files = sorted(workdir.glob("*/*.arc"))
    n = len(arc_files)
    if n == 0:
        logger.warning("No .arc files found under %s — nothing to convert.", workdir)
        return pd.DataFrame(
            columns=["compound_id", "arc_file", "sdf_file", "status", "error_message", "duration_s"]
        )

    logger.info("Found %d .arc file(s) to convert.", n)

    run_kwargs = dict(timeout=timeout, overwrite=overwrite, cancel_event=cancel_event)
    results: list[Optional[SdfResult]] = [None] * n

    if max_workers == 1:
        for idx, arc_file in enumerate(arc_files):
            sdf_file = arc_file.parent / f"{arc_file.stem}.sdf"
            logger.info("[%d/%d] %s", idx + 1, n, arc_file.stem)
            try:
                results[idx] = _convert_arc_to_sdf(arc_file, sdf_file, **run_kwargs)
            except Exception as exc:  # pylint: disable=broad-except
                logger.exception("[%s] Unexpected exception: %s", arc_file.stem, exc)
                results[idx] = SdfResult(
                    compound_id=arc_file.stem, arc_file=arc_file,
                    status="failed", error_message=f"Unexpected exception: {exc}",
                )
    else:
        logger.info("Parallel SDF export: %d files across %d workers.", n, max_workers)
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            future_to_idx = {
                pool.submit(
                    _convert_arc_to_sdf,
                    arc_file,
                    arc_file.parent / f"{arc_file.stem}.sdf",
                    **run_kwargs,
                ): idx
                for idx, arc_file in enumerate(arc_files)
            }
            completed = 0
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                completed += 1
                try:
                    results[idx] = future.result()
                except Exception as exc:  # pylint: disable=broad-except
                    logger.exception(
                        "[%s] Unexpected exception: %s", arc_files[idx].stem, exc,
                    )
                    results[idx] = SdfResult(
                        compound_id=arc_files[idx].stem, arc_file=arc_files[idx],
                        status="failed", error_message=f"Unexpected exception: {exc}",
                    )
                logger.info(
                    "[%d/%d] %s — %s",
                    completed, n, results[idx].compound_id, results[idx].status,
                )

    results_df = pd.DataFrame([r.to_dict() for r in results]) if results else pd.DataFrame(
        columns=["compound_id", "arc_file", "sdf_file", "status", "error_message", "duration_s"]
    )

    csv_path = workdir / "sdf_results.csv"
    results_df.to_csv(csv_path, index=False)
    _log_batch_finished("SDF export", results_df, csv_path, cancel_event)
    return results_df


# ---------------------------------------------------------------------------
# Retry: re-run failed/timed-out jobs with extended settings
# ---------------------------------------------------------------------------


def _merge_retry_results(
    existing_df: pd.DataFrame,
    new_df: pd.DataFrame,
    key: str = "compound_id",
) -> pd.DataFrame:
    """Merge *new_df* rows into *existing_df* by *key*, fully overwriting matches.

    Unlike ``DataFrame.update()``, this also overwrites columns that go back
    to ``None``/NaN in *new_df* (e.g. ``error_message`` clearing to ``None``
    on a successful retry) instead of silently keeping the stale value —
    ``update()`` treats NaN in the incoming frame as "leave the old value
    alone," which would leave a recovered compound's old error message in
    place.

    Args:
        existing_df: Previously saved results (e.g. loaded from
                     ``batch_results.csv``). May be empty.
        new_df:      Fresh results for a subset of the same *key* values.
        key:         Column to align rows on.

    Returns:
        *existing_df* with matching rows overwritten and any brand-new keys
        appended.
    """
    if new_df.empty:
        return existing_df
    if existing_df.empty:
        return new_df.reset_index(drop=True)

    existing = existing_df.set_index(key)
    new = new_df.set_index(key)

    for col in new.columns:
        if col not in existing.columns:
            existing[col] = None
        else:
            # Columns loaded from CSV that were all-NaN in `existing` (e.g.
            # error_message on a prior timeout) get an incompatible numeric
            # dtype (float64). Widen to object so string values from `new`
            # can be assigned without pandas raising a dtype-coercion error.
            existing[col] = existing[col].astype(object)

    matched = new.index.intersection(existing.index)
    new_rows = new.index.difference(existing.index)

    if len(matched) > 0:
        existing.loc[matched, new.columns] = new.loc[matched, new.columns]
    if len(new_rows) > 0:
        existing = pd.concat([existing, new.loc[new_rows]])

    return existing.reset_index()


def retry_failed_batch(
    df: pd.DataFrame,
    smiles_col: str,
    id_col: Optional[str] = None,
    keywords: str = "PM7 PRECISE GEO-OK",
    workdir: "str | Path" = "mopac_jobs",
    charge_col: Optional[str] = None,
    multiplicity_col: Optional[str] = None,
    auto_charge: bool = False,
    coord_format: str = OBABEL_3D_FORMAT,
    obabel_options: Optional[list[str]] = None,
    obabel_timeout: int = 240,
    mopac_timeout: Optional[int] = 1800,
    retry_keywords: Optional[str] = None,
    max_workers: int = 1,
    cancel_event: Optional[threading.Event] = None,
) -> pd.DataFrame:
    """Re-run compounds that failed, timed out, or were cancelled in a
    previous batch.

    Reads ``<workdir>/prepare_results.csv`` and ``<workdir>/batch_results.csv``
    (as written by :func:`prepare_batch` / :func:`run_mopac_batch`) to find
    which compounds still need work, retries only those — using extended
    timeouts and, optionally, a different MOPAC keyword string — and merges
    the outcomes back into both CSVs. Compounds that already succeeded are
    left untouched.

    Two independent retry passes happen:

    1. Compounds with ``status in {"failed", "cancelled"}`` in
       ``prepare_results.csv`` (Open Babel never produced a ``.mop`` file)
       are retried with *obabel_timeout*.
    2. Compounds with ``status in {"failed", "timeout", "cancelled"}`` in
       ``batch_results.csv``, plus any compound that just became ``ready``
       in step 1, are retried with *mopac_timeout*.

    A compound that was previously stopped via *cancel_event* (i.e. ended up
    with ``status="cancelled"``) is therefore retried exactly like a
    genuine failure or timeout — stopping a run and retrying afterward
    compose naturally with no extra bookkeeping needed.

    Args:
        df:               The **same** DataFrame (same rows, same order)
                          passed to the original :func:`prepare_batch` call.
                          Compound-ID assignment (auto sequential IDs,
                          duplicate-ID suffixing) depends on enumerating the
                          full DataFrame in its original order — a filtered
                          or reordered DataFrame will not line up with the
                          existing job directories.
        smiles_col:       Column name holding SMILES strings.
        id_col:           Same *id_col* used for the original run.
        keywords:         Same *keywords* used for the original run. Used to
                          rebuild charge/multiplicity/ID assignment for every
                          job; superseded by *retry_keywords* (if given) when
                          actually writing retried ``.mop`` files.
        workdir:          Root directory used for the previous pipeline runs.
        charge_col:       Same *charge_col* used for the original run.
        multiplicity_col: Same *multiplicity_col* used for the original run.
        auto_charge:      Same *auto_charge* used for the original run.
        coord_format:     Same *coord_format* used for the original run.
        obabel_options:   Extra obabel CLI flags for the retry's Open Babel
                          attempts. Defaults to ``OBABEL_DEFAULT_OPTIONS``.
        obabel_timeout:   Per-compound Open Babel timeout for the retry pass.
        mopac_timeout:    Per-compound MOPAC timeout for the retry pass.
        retry_keywords:   When given, used instead of *keywords* for every
                          retried compound (e.g. to relax convergence
                          criteria). ``None`` reuses *keywords* unchanged.
        max_workers:      Number of parallel worker threads.
        cancel_event:     When set (from another thread) while this is
                          running, any compound not yet started is skipped
                          with ``status="cancelled"`` and the in-flight
                          obabel/MOPAC process(es) (if any) are killed
                          immediately.

    Returns:
        The updated ``batch_results.csv`` DataFrame (same schema as
        :func:`run_mopac_batch`), covering all compounds — not just the
        retried ones. Also saved back to ``<workdir>/batch_results.csv``
        (and ``<workdir>/prepare_results.csv`` is updated in place too).

    Raises:
        EnvironmentError:  If ``obabel`` or ``mopac`` is not found in PATH.
        FileNotFoundError: If *workdir* does not exist.
        KeyError:          If a referenced column is absent from *df*.
    """
    validate_environment()

    workdir = Path(workdir).resolve()
    if not workdir.exists():
        raise FileNotFoundError(f"workdir does not exist: {workdir}")

    optional_cols = [c for c in (id_col, charge_col, multiplicity_col) if c]
    missing_cols = [c for c in [smiles_col, *optional_cols] if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Column(s) not found in DataFrame: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    prepare_csv = workdir / "prepare_results.csv"
    batch_csv = workdir / "batch_results.csv"

    prepare_cols = ["compound_id", "smiles", "status", "mop_file", "error_message", "duration_s"]
    batch_cols = [
        "compound_id", "smiles", "status",
        "mop_file", "out_file", "arc_file",
        "error_message", "return_code", "duration_s",
    ]

    prepare_df = pd.read_csv(prepare_csv) if prepare_csv.exists() else pd.DataFrame(columns=prepare_cols)
    batch_df = pd.read_csv(batch_csv) if batch_csv.exists() else pd.DataFrame(columns=batch_cols)

    obabel_retry_ids = (
        set(prepare_df.loc[prepare_df["status"].isin(["failed", "cancelled"]), "compound_id"])
        if not prepare_df.empty else set()
    )
    mopac_retry_ids = (
        set(batch_df.loc[batch_df["status"].isin(["failed", "timeout", "cancelled"]), "compound_id"])
        if not batch_df.empty else set()
    )

    if not obabel_retry_ids and not mopac_retry_ids:
        logger.info("Nothing to retry — no failed/timeout compounds found in %s.", workdir)
        return batch_df

    effective_keywords = retry_keywords if retry_keywords else keywords

    jobs = _build_job_list(
        df, smiles_col, id_col, effective_keywords, workdir,
        charge_col, multiplicity_col, auto_charge=auto_charge,
    )
    jobs_by_id = {j.compound_id: j for j in jobs}

    missing = (obabel_retry_ids | mopac_retry_ids) - jobs_by_id.keys()
    if missing:
        logger.warning(
            "Retry: %d compound id(s) from previous results not found in the "
            "current DataFrame — skipping: %s", len(missing), sorted(missing),
        )

    # -- Phase-1 retry: obabel-stage failures --------------------------------
    obabel_jobs = [jobs_by_id[i] for i in obabel_retry_ids if i in jobs_by_id]
    prep_retry_results: list[PrepResult] = []

    if obabel_jobs:
        logger.info(
            "Retrying %d obabel-stage failure(s) with timeout=%ds.",
            len(obabel_jobs), obabel_timeout,
        )
        prep_kwargs = dict(
            overwrite=True, obabel_timeout=obabel_timeout,
            obabel_options=obabel_options, coord_format=coord_format,
            cancel_event=cancel_event,
        )
        if max_workers == 1:
            for job in obabel_jobs:
                try:
                    prep_retry_results.append(prepare_one_molecule(job=job, **prep_kwargs))
                except Exception as exc:  # pylint: disable=broad-except
                    logger.exception("[%s] Unexpected exception during retry: %s", job.compound_id, exc)
                    prep_retry_results.append(PrepResult(
                        compound_id=job.compound_id, smiles=job.smiles,
                        status="failed", error_message=f"Unexpected exception: {exc}",
                    ))
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                future_to_job = {
                    pool.submit(prepare_one_molecule, job=job, **prep_kwargs): job
                    for job in obabel_jobs
                }
                for future in as_completed(future_to_job):
                    job = future_to_job[future]
                    try:
                        prep_retry_results.append(future.result())
                    except Exception as exc:  # pylint: disable=broad-except
                        logger.exception("[%s] Unexpected exception during retry: %s", job.compound_id, exc)
                        prep_retry_results.append(PrepResult(
                            compound_id=job.compound_id, smiles=job.smiles,
                            status="failed", error_message=f"Unexpected exception: {exc}",
                        ))

        prep_retry_df = pd.DataFrame([r.to_dict() for r in prep_retry_results])
        prepare_df = _merge_retry_results(prepare_df, prep_retry_df)
        prepare_df.to_csv(prepare_csv, index=False)

    newly_ready_ids = {r.compound_id for r in prep_retry_results if r.status == "ready"}

    # -- Phase-2 retry: mopac-stage failures + freshly-prepared retries -----
    mopac_targets = mopac_retry_ids | newly_ready_ids
    mopac_target_jobs = [
        jobs_by_id[i] for i in mopac_targets
        if i in jobs_by_id and (workdir / i / f"{i}.mop").exists()
    ]

    mopac_retry_results: list[JobResult] = []

    if mopac_target_jobs:
        logger.info(
            "Retrying %d MOPAC-stage job(s) with timeout=%s.",
            len(mopac_target_jobs), mopac_timeout,
        )
        run_kwargs = dict(timeout=mopac_timeout, overwrite=True, cancel_event=cancel_event)

        if max_workers == 1:
            for job in mopac_target_jobs:
                mop_file = job.job_dir / f"{job.compound_id}.mop"
                try:
                    mopac_retry_results.append(_run_mopac_one(mop_file, **run_kwargs))
                except Exception as exc:  # pylint: disable=broad-except
                    logger.exception("[%s] Unexpected exception during retry: %s", job.compound_id, exc)
                    mopac_retry_results.append(JobResult(
                        compound_id=job.compound_id, smiles=job.smiles,
                        status="failed", error_message=f"Unexpected exception: {exc}",
                    ))
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                future_to_job = {
                    pool.submit(_run_mopac_one, job.job_dir / f"{job.compound_id}.mop", **run_kwargs): job
                    for job in mopac_target_jobs
                }
                for future in as_completed(future_to_job):
                    job = future_to_job[future]
                    try:
                        mopac_retry_results.append(future.result())
                    except Exception as exc:  # pylint: disable=broad-except
                        logger.exception("[%s] Unexpected exception during retry: %s", job.compound_id, exc)
                        mopac_retry_results.append(JobResult(
                            compound_id=job.compound_id, smiles=job.smiles,
                            status="failed", error_message=f"Unexpected exception: {exc}",
                        ))

        mopac_retry_df = pd.DataFrame([r.to_dict() for r in mopac_retry_results])
        batch_df = _merge_retry_results(batch_df, mopac_retry_df)
        batch_df.to_csv(batch_csv, index=False)

    recovered = sum(1 for r in mopac_retry_results if r.status == "success")
    still_failed_obabel = len(obabel_jobs) - len(newly_ready_ids)
    stopped_note = " (stopped by user)" if cancel_event is not None and cancel_event.is_set() else ""
    logger.info(
        "Retry finished%s — MOPAC: %d/%d recovered. Obabel: %d/%d still failing. "
        "Results → %s",
        stopped_note, recovered, len(mopac_target_jobs), still_failed_obabel, len(obabel_jobs), batch_csv,
    )

    return batch_df


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

    # Step 1: generate 3D structures and write .mop files.
    prep = prepare_batch(
        df=sample_df,
        smiles_col="smiles",
        id_col="compound_id",
        keywords="PM7 PRECISE GEO-OK",
        workdir="mopac_jobs",
        overwrite=True,
    )

    print("\n=== Prepare Results ===")
    print(prep[["compound_id", "status", "duration_s"]].to_string(index=False))

    # Step 2: run MOPAC on all prepared .mop files.
    results = run_mopac_batch(
        workdir="mopac_jobs",
        overwrite=True,
        timeout=900,
    )

    print("\n=== MOPAC Results ===")
    with pd.option_context("display.max_colwidth", 60):
        print(results[["compound_id", "status", "duration_s"]].to_string(index=False))
