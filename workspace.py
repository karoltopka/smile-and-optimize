"""Workspace persistence for the SMILE and optimize GUI.

A *workspace* is a small JSON sidecar written into the working directory --
``<workdir>/workspace.json`` -- recording everything the GUI holds in memory
that the working directory itself cannot tell you: which input file the
compounds came from, which columns were used, the Settings values, and how
far each phase got.

Why a sidecar rather than a database or a pickle: the pipeline in
mopac_batch_runner.py is *already* resumable from the filesystem alone --
``run_mopac_batch`` globs ``<workdir>/*/*.mop`` and reads each compound's ID
and SMILES back out of the .mop header, and every phase skips finished
compounds when ``overwrite=False``. The only thing genuinely lost when the
app closes is the GUI-side context, and that is a few kilobytes of plain
values. Keeping it as readable JSON next to the CSVs means a workspace can
be inspected, hand-edited, or thrown away without consequence, and moving or
zipping the working directory keeps it consistent.

Deliberately free of any Qt import: this module is exercised by the GUI but
does not depend on it, so the format can be tested from a plain REPL.
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

WORKSPACE_FILENAME = "workspace.json"

# Bumped only for changes a previous version could not read correctly.
# load_workspace() refuses lower-and-incompatible, and warns-but-continues on
# higher, since a newer file is usually just additive.
SCHEMA_VERSION = 1

PHASE_KEYS = ("prepare", "mopac", "export")

# Which results CSV mopac_batch_runner writes for each phase. Used on restore
# to repopulate the results table without re-running anything.
PHASE_RESULT_CSV = {
    "prepare": "prepare_results.csv",
    "mopac": "batch_results.csv",
    "export": "sdf_results.csv",
}

# Phase states. "running" is persisted at phase *start*, which is what makes a
# process killed mid-run recoverable: on restore it reads back as interrupted
# rather than as never-started.
STATE_NOT_STARTED = "not_started"
STATE_RUNNING = "running"
STATE_COMPLETED = "completed"
STATE_STOPPED = "stopped"

# Read size for input-file hashing. A 20k-row .xlsx is a few MB, so this is
# only ever a handful of iterations, but streaming costs nothing and keeps a
# pathological input from being slurped whole.
_HASH_CHUNK = 1 << 20


class WorkspaceError(Exception):
    """Raised when a workspace file is missing, unreadable, or not a workspace."""


def now_iso() -> str:
    """Timestamp for ``saved_at`` / ``started_at`` / ``finished_at`` fields."""
    return datetime.now().isoformat(timespec="seconds")


# ----------------------------------------------------------------------
# Input-file fingerprinting
# ----------------------------------------------------------------------
def fingerprint_input(path: "str | Path") -> dict:
    """Hash and stat an input file so a later restore can detect edits.

    This matters more than it looks. Compound IDs are not stored anywhere
    durable -- ``_build_job_list`` in mopac_batch_runner.py derives them by
    enumerating the whole DataFrame in order (sequential ``mol_000001`` names
    for missing IDs, ``_1``/``_2`` suffixes to break sanitised duplicates),
    and ``retry_failed_batch`` re-derives them the same way. So a restored
    session pointed at a reordered or edited spreadsheet would silently map
    compounds onto the wrong job directories.

    Returns an empty dict if the file cannot be read, so a fingerprint can
    never be the thing that stops a workspace from being saved.
    """
    path = Path(path)
    try:
        stat = path.stat()
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(_HASH_CHUNK), b""):
                digest.update(chunk)
    except OSError:
        return {}
    return {
        "sha256": digest.hexdigest(),
        "size_bytes": stat.st_size,
        "mtime": stat.st_mtime,
    }


def input_matches(saved: dict, current: dict) -> bool:
    """True when two fingerprints describe the same file contents.

    Compares the hash alone -- size and mtime are recorded for diagnostics
    only. A file copied or re-saved without content changes must still count
    as a match, otherwise every restore from a synced folder would warn.
    """
    saved_hash = saved.get("sha256")
    current_hash = current.get("sha256")
    if not saved_hash or not current_hash:
        return True  # nothing to compare -- don't cry wolf
    return saved_hash == current_hash


# ----------------------------------------------------------------------
# Save / load
# ----------------------------------------------------------------------
def workspace_path_for(workdir: "str | Path") -> Path:
    """Canonical workspace location for a working directory."""
    return Path(workdir) / WORKSPACE_FILENAME


def new_workspace(app_version: str = "") -> dict:
    """An empty workspace dict with every expected key present.

    Callers fill in what they know; absent values stay at these defaults so
    restore never has to guard against missing keys.
    """
    return {
        "schema_version": SCHEMA_VERSION,
        "app_version": app_version,
        "saved_at": now_iso(),
        "workdir": "",
        "input": {},
        "settings": {},
        "settings_locked": False,
        "phases": {key: new_phase(key) for key in PHASE_KEYS},
        "workdir_scan": {},
    }


def new_phase(key: str) -> dict:
    return {
        "state": STATE_NOT_STARTED,
        "started_at": None,
        "finished_at": None,
        "counts": {},
        "results_csv": PHASE_RESULT_CSV.get(key, ""),
    }


def save_workspace(path: "str | Path", data: dict) -> None:
    """Write *data* to *path* atomically.

    Same temp-file-then-``os.replace`` approach as ``_write_partial_csv`` in
    mopac_batch_runner.py, and for the same reason: the GUI force-exits via
    ``os._exit`` when a phase is still winding down at close, so a write must
    never be interruptible into a half-written file. A truncated workspace
    would be worse than none at all -- it would fail to parse on the next
    launch, right when the user is trying to recover a long run.

    Raises OSError; callers in the GUI treat a failed save as a warning in
    the log rather than an error dialog.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {**data, "saved_at": now_iso(), "schema_version": SCHEMA_VERSION}
    tmp_path = path.with_name(path.name + ".tmp")
    tmp_path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    os.replace(tmp_path, path)


def load_workspace(path: "str | Path") -> dict:
    """Read a workspace file, normalising it to the current schema shape.

    Raises:
        WorkspaceError: if the file is missing, unreadable, not valid JSON,
            not a JSON object, or was written by an incompatible schema.
    """
    path = Path(path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise WorkspaceError(f"No workspace file at {path}") from exc
    except OSError as exc:
        raise WorkspaceError(f"Could not read {path}: {exc}") from exc
    except ValueError as exc:
        raise WorkspaceError(f"{path.name} is not valid JSON: {exc}") from exc

    if not isinstance(raw, dict) or "schema_version" not in raw:
        raise WorkspaceError(f"{path.name} does not look like a workspace file.")

    version = raw.get("schema_version")
    if not isinstance(version, int) or version < 1:
        raise WorkspaceError(f"{path.name} has an unrecognised schema_version: {version!r}")

    # Fill in anything this version expects but the file predates, so the
    # rest of the GUI can read keys unconditionally.
    merged = new_workspace()
    merged.update(raw)
    phases = dict(merged.get("phases") or {})
    for key in PHASE_KEYS:
        phase = dict(new_phase(key))
        phase.update(phases.get(key) or {})
        phases[key] = phase
    merged["phases"] = phases
    for key in ("input", "settings", "workdir_scan"):
        if not isinstance(merged.get(key), dict):
            merged[key] = {}
    return merged


def is_future_version(data: dict) -> bool:
    """True when the file was written by a newer build than this one."""
    try:
        return int(data.get("schema_version", 0)) > SCHEMA_VERSION
    except (TypeError, ValueError):
        return False


# ----------------------------------------------------------------------
# Filesystem scan
# ----------------------------------------------------------------------
def scan_workdir(workdir: "str | Path", cancel_event: Optional[Any] = None) -> dict:
    """Count the pipeline artefacts actually present under *workdir*.

    The results CSVs are only rewritten every ``_CHECKPOINT_EVERY`` (200)
    completed jobs, so after a hard exit they can be missing up to 199
    compounds whose .mop/.out files are sitting on disk. Restoring from the
    CSVs alone therefore under-reports progress. This walks the job
    directories instead and gives the true picture.

    One ``os.scandir`` per job directory, extensions counted by suffix -- at
    20,000 folders that is seconds, not minutes, but it is still slow enough
    that the GUI runs it on a background thread (see MainWindow's "Scan"
    pseudo-phase) rather than blocking the event loop.

    Args:
        workdir:      Working directory to scan.
        cancel_event: Optional ``threading.Event``; polled between job
                      directories so a scan over a slow network share can be
                      stopped with the STOP button. A cancelled scan returns
                      the counts gathered so far with ``cancelled: True``.

    Returns:
        Dict with ``job_dirs``, ``mop``, ``out``, ``arc``, ``sdf``,
        ``scanned_at``, and ``cancelled``.
    """
    workdir = Path(workdir)
    counts = {"job_dirs": 0, "mop": 0, "out": 0, "arc": 0, "sdf": 0}
    cancelled = False

    if workdir.is_dir():
        tracked = {".mop": "mop", ".out": "out", ".arc": "arc", ".sdf": "sdf"}
        try:
            with os.scandir(workdir) as top_level:
                job_dirs = [entry for entry in top_level if entry.is_dir()]
        except OSError:
            job_dirs = []

        for job_dir in job_dirs:
            if cancel_event is not None and cancel_event.is_set():
                cancelled = True
                break
            counts["job_dirs"] += 1
            try:
                with os.scandir(job_dir.path) as files:
                    for entry in files:
                        key = tracked.get(Path(entry.name).suffix.lower())
                        if key is not None and entry.is_file():
                            counts[key] += 1
            except OSError:
                # An unreadable job directory is worth skipping, not worth
                # failing the whole scan over.
                continue

    counts["scanned_at"] = now_iso()
    counts["cancelled"] = cancelled
    return counts


def describe_scan(scan: dict) -> str:
    """One-line human summary of a scan result, for the summary label."""
    if not scan:
        return "Working directory not scanned yet."
    text = (
        f"Working directory: {scan.get('job_dirs', 0):,} job folder(s) — "
        f"{scan.get('mop', 0):,} .mop, {scan.get('out', 0):,} .out, "
        f"{scan.get('arc', 0):,} .arc, {scan.get('sdf', 0):,} .sdf"
    )
    if scan.get("cancelled"):
        text += "  (scan stopped early)"
    return text


def describe_workspace(data: dict) -> str:
    """Multi-line summary shown in the "restore this workspace?" dialog."""
    lines = []
    saved_at = data.get("saved_at")
    if saved_at:
        lines.append(f"Saved: {saved_at}")

    source = (data.get("input") or {}).get("path")
    if source:
        rows = (data.get("input") or {}).get("n_rows")
        lines.append(f"Input: {Path(source).name}" + (f" ({rows:,} rows)" if rows else ""))

    labels = {"prepare": "Prepare", "mopac": "Run MOPAC", "export": "Export SDF"}
    for key in PHASE_KEYS:
        phase = (data.get("phases") or {}).get(key) or {}
        state = phase.get("state", STATE_NOT_STARTED)
        if state == STATE_NOT_STARTED:
            continue
        counts = phase.get("counts") or {}
        detail = ", ".join(f"{name}: {value:,}" for name, value in counts.items())
        readable = {
            STATE_RUNNING: "interrupted while running",
            STATE_COMPLETED: "completed",
            STATE_STOPPED: "stopped by user",
        }.get(state, state)
        lines.append(f"{labels[key]}: {readable}" + (f" — {detail}" if detail else ""))

    scan = data.get("workdir_scan") or {}
    if scan:
        lines.append(describe_scan(scan))

    return "\n".join(lines) if lines else "No phases have been run yet."


def has_unfinished_work(data: dict) -> bool:
    """True when at least one phase was interrupted or reported failures.

    Drives the restore dialog's default action: an interrupted run should be
    resumed with ``overwrite`` off so finished compounds are skipped, rather
    than redone from scratch.
    """
    for key in PHASE_KEYS:
        phase = (data.get("phases") or {}).get(key) or {}
        if phase.get("state") in (STATE_RUNNING, STATE_STOPPED):
            return True
        counts = phase.get("counts") or {}
        if any(counts.get(status) for status in ("failed", "timeout", "cancelled")):
            return True
    return False
