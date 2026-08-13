"""SMILE and optimize -- desktop GUI for mopac_batch_runner.py.

Wraps the three-phase pipeline (prepare -> run MOPAC -> export SDF) in a
PySide6 window so it can be driven without a Jupyter notebook open. Run
directly with:

    py gui.py
"""

from __future__ import annotations

import functools
import json
import logging
import os
import queue
import shutil
import sys
import threading
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QSettings, Qt, QTimer, QUrl
from PySide6.QtGui import (
    QAction,
    QColor,
    QDesktopServices,
    QFontDatabase,
    QIcon,
    QKeySequence,
    QPixmap,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import workspace as ws
from mopac_batch_runner import (
    export_sdf,
    prepare_batch,
    retry_failed_batch,
    run_mopac_batch,
    validate_environment,
)

DEFAULT_WORKERS = min(4, os.cpu_count() or 4)
ASSETS_DIR = Path(__file__).resolve().parent / "assets"
ICON_PATH = ASSETS_DIR / "icon.png"
README_PATH = Path(__file__).resolve().parent / "README.md"

# QSettings scope for things that belong to the *application* rather than to
# any one workspace: the recent-workspaces list. Everything that describes a
# run lives in <workdir>/workspace.json instead, so a working directory stays
# self-describing when it is copied to another machine.
_ORG_NAME = "SMILE and optimize"
_APP_NAME = "SMILE and optimize"
_RECENT_KEY = "recentWorkspaces"
_MAX_RECENT = 8

# Single source of truth for the Settings group: (workspace key, MainWindow
# attribute, kind). Drives workspace save/restore AND the SAVE SETTINGS lock,
# so a new setting only has to be added in one place. Previously the lock kept
# its own hand-maintained copy of this list.
_SETTINGS_FIELDS: tuple[tuple[str, str, str], ...] = (
    ("keywords", "keywords_edit", "text"),
    ("auto_charge", "auto_charge_check", "bool"),
    ("overwrite", "overwrite_check", "bool"),
    ("prepare_workers", "prepare_workers_spin", "int"),
    ("mopac_workers", "mopac_workers_spin", "int"),
    ("export_workers", "export_workers_spin", "int"),
    ("obabel_timeout", "obabel_timeout_spin", "int"),
    ("mopac_timeout", "mopac_timeout_spin", "int"),
    ("export_timeout", "export_timeout_spin", "int"),
    ("retry_obabel_timeout", "retry_obabel_timeout_spin", "int"),
    ("retry_mopac_timeout", "retry_mopac_timeout_spin", "int"),
    ("retry_keywords", "retry_keywords_edit", "text"),
    ("copy_sdf", "copy_sdf_check", "bool"),
    ("sdf_catalogue", "sdf_catalogue_edit", "text"),
)

# sdf_catalogue_edit is excluded: its enabled state is owned by
# copy_sdf_check (see _build_settings_group), so the lock must not set it
# blindly -- see _on_toggle_settings_lock_clicked.
_LOCKABLE_SETTINGS: tuple[str, ...] = tuple(
    attr for _, attr, _ in _SETTINGS_FIELDS if attr != "sdf_catalogue_edit"
)

# Phase name -> workspace phase key. The phase names are the same magic
# strings _poll_worker already dispatches on. "Retry Failed" maps to the
# MOPAC phase because that is the CSV it rewrites, and "Scan" is a
# pseudo-phase with no results of its own.
_PHASE_TO_WORKSPACE_KEY = {
    "Prepare": "prepare",
    "Run MOPAC": "mopac",
    "Export SDF": "export",
    "Retry Failed": "mopac",
}
_SCAN_PHASE_NAME = "Scan"

# Seconds a running phase gets to wind down after the window is closed before
# the process exits regardless. See MainWindow.closeEvent.
_FORCED_EXIT_GRACE_S: float = 10.0

PROVENANCE_PATH = Path(__file__).resolve().parent / "provenance.json"
# Kept in sync with provenance.json; used only if that file is missing or
# malformed so a broken/deleted provenance file can't stop the app launching.
_DEFAULT_PROVENANCE = {
    "author": "Maciej Gromelski",
    "version": "1.3.0",
    "year": 2026,
    "contact_email": "maciek.gromelski@gmail.com",
}


def _load_provenance() -> dict:
    """Load author/version/year/contact info for the footer from provenance.json."""
    try:
        data = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
        return {**_DEFAULT_PROVENANCE, **data}
    except (OSError, ValueError) as exc:
        print(
            f"Warning: could not load {PROVENANCE_PATH.name} ({exc}); using built-in defaults.",
            file=sys.stderr,
        )
        return dict(_DEFAULT_PROVENANCE)


PROVENANCE = _load_provenance()

APP_STYLESHEET = """
/* ============================================================
   Base / page
   ============================================================ */
QMainWindow {
    background-color: #F5F5F5;
}
QWidget {
    background-color: #F5F5F5;
    color: #303841;
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 13px;
}

/* ============================================================
   Menu bar (see MainWindow._build_menu_bar()). Without these rules the
   menus fall through to the generic QWidget rule above -- a flat bar and
   a borderless popup with no hover feedback.

   The bottom border is deliberately a quiet 1px grey, NOT another orange
   rule: headerBar's 2px #FF5722 border-bottom sits directly underneath
   and has to stay the page's dominant divider.
   ============================================================ */
QMenuBar {
    background-color: #F5F5F5;
    color: #303841;
    border-bottom: 1px solid #DADADA;
    padding: 2px 8px;
}
QMenuBar::item {
    background: transparent;
    padding: 5px 12px;
    border-radius: 4px;
}
QMenuBar::item:selected {
    background-color: #D9EBEC;
    color: #185861;
}
QMenuBar::item:pressed {
    background-color: #76ABAE;
    color: #F5F5F5;
}
QMenu {
    background-color: #FFFFFF;
    color: #303841;
    border: 1px solid #C7D3D4;
    border-radius: 6px;
    padding: 4px;
}
QMenu::item {
    padding: 6px 28px 6px 20px;
    border-radius: 4px;
}
QMenu::item:selected {
    background-color: #76ABAE;
    color: #F5F5F5;
}
QMenu::item:disabled {
    color: #9AA0A2;
}
QMenu::separator {
    height: 1px;
    background-color: #DDE4E4;
    margin: 4px 10px;
}

/* ============================================================
   Header bar (see MainWindow._build_header()): icon on the left, a
   vertical accent divider, then a title/subtitle stack. Title/divider
   color (#185861) and subtitle color (#FF5722) were sampled directly
   from assets/logo2.png and assets/icon.png rather than picked freehand.
   ============================================================ */
QWidget#headerBar {
    background-color: #F5F5F5;
    border-bottom: 2px solid #FF5722;
}
QLabel#headerLogo {
    background: transparent;
}
QFrame#headerDivider {
    background-color: #185861;
    border: none;
}
QLabel#headerTitle {
    background: transparent;
    color: #185861;
    font-size: 22px;
    font-weight: 700;
    padding: 2px 0px;
}
QLabel#headerSubtitle {
    background: transparent;
    color: #FF5722;
    font-size: 12px;
    font-weight: 600;
}

/* ============================================================
   Group boxes ("cards")
   ============================================================ */
QGroupBox {
    background-color: #F5F5F5;
    border: 1px solid #C7D3D4;
    border-left: 3px solid #76ABAE;
    border-radius: 6px;
    margin-top: 14px;
    padding: 16px 12px 12px 12px;
    font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 6px;
    color: #303841;
    background-color: #F5F5F5;
}

/* ============================================================
   Buttons -- secondary/utility (default QPushButton)
   ============================================================ */
QPushButton {
    background-color: #FFFFFF;
    color: #303841;
    border: 1px solid #76ABAE;
    border-radius: 6px;
    padding: 6px 14px;
    font-weight: 500;
}
QPushButton:hover {
    background-color: #EAF3F3;
    border: 1px solid #5C8A8D;
}
QPushButton:pressed {
    background-color: #D9EBEC;
    border: 1px solid #303841;
}
QPushButton:disabled {
    background-color: #D9D9D9;
    color: #9AA0A2;
    border: 1px solid #C7C7C7;
}

/* ============================================================
   Buttons -- primary phase actions (objectName="primaryButton")
   NOTE: :disabled kept LAST in this block -- see plan risk notes.
   ============================================================ */
QPushButton#primaryButton {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                 stop:0 #8FC1C4, stop:1 #76ABAE);
    color: #303841;
    border: 1px solid #5C8A8D;
    border-radius: 6px;
    padding: 8px 18px;
    font-weight: 700;
}
QPushButton#primaryButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                 stop:0 #FF8A65, stop:1 #76ABAE);
    border: 1px solid #FF5722;
    color: #303841;
}
QPushButton#primaryButton:pressed {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                 stop:0 #5C8A8D, stop:1 #4A7679);
    border: 1px solid #303841;
    color: #F5F5F5;
    padding-top: 9px;
    padding-bottom: 7px;
}
QPushButton#primaryButton:disabled {
    background: #D9D9D9;
    color: #9AA0A2;
    border: 1px solid #C7C7C7;
}

/* ============================================================
   STOP button (objectName="stopButton") -- kills the running phase's
   in-flight obabel/MOPAC process(es) immediately. Disabled (grey, same
   convention as the other buttons' :disabled state) whenever no phase is
   running; red outline/fill otherwise so it reads as a danger action.
   ============================================================ */
QPushButton#stopButton {
    background-color: #FFFFFF;
    color: #B00020;
    border: 2px solid #B00020;
    border-radius: 6px;
    padding: 8px 18px;
    font-weight: 700;
}
QPushButton#stopButton:hover {
    background-color: #FDECEA;
    border: 2px solid #8E0018;
}
QPushButton#stopButton:pressed {
    background-color: #B00020;
    color: #FFFFFF;
}
QPushButton#stopButton:disabled {
    background-color: #D9D9D9;
    color: #9AA0A2;
    border: 2px solid #C7C7C7;
}

/* ============================================================
   Settings lock/unlock toggle (objectName="lockToggleButton")
   This is the "unlocked / editable" baseline look; the "locked" look
   is applied as a per-instance override in
   MainWindow._on_toggle_settings_lock_clicked (same pattern as
   env_status_label's dynamic color) since it's just two fixed states.
   ============================================================ */
QPushButton#lockToggleButton {
    background-color: #FFFFFF;
    color: #303841;
    border: 2px solid #76ABAE;
    border-radius: 6px;
    padding: 8px 28px;
    font-weight: 700;
}
QPushButton#lockToggleButton:hover {
    background-color: #EAF3F3;
    border: 2px solid #5C8A8D;
}
QPushButton#lockToggleButton:pressed {
    background-color: #D9EBEC;
}

/* ============================================================
   Progress bar
   ============================================================ */
QProgressBar {
    background-color: #E4E4E4;
    border: 1px solid #C7D3D4;
    border-radius: 6px;
    text-align: center;
    color: #303841;
    height: 18px;
}
QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                 stop:0 #76ABAE, stop:1 #303841);
    border-radius: 5px;
}

/* ============================================================
   Table + header
   ============================================================ */
QTableWidget {
    background-color: #FFFFFF;
    alternate-background-color: #F0F4F4;
    gridline-color: #DDE4E4;
    border: 1px solid #C7D3D4;
    selection-background-color: #76ABAE;
    selection-color: #F5F5F5;
}
QHeaderView::section {
    background-color: #303841;
    color: #F5F5F5;
    padding: 6px;
    border: none;
    border-right: 1px solid #4A5560;
    font-weight: 600;
}
QTableCornerButton::section {
    background-color: #303841;
    border: none;
}

/* ============================================================
   Inputs
   ============================================================ */
QLineEdit, QComboBox {
    background-color: #FFFFFF;
    border: 1px solid #C7D3D4;
    border-radius: 4px;
    padding: 4px 6px;
    selection-background-color: #76ABAE;
    selection-color: #F5F5F5;
}
QLineEdit:focus, QComboBox:focus {
    border: 1px solid #76ABAE;
}
QLineEdit:disabled, QComboBox:disabled {
    background-color: #ECECEC;
    color: #9AA0A2;
    border: 1px solid #DADADA;
}

/* ============================================================
   Checkboxes are deliberately left fully unstyled -- QSS overrides on
   QCheckBox's indicator were found, through direct testing, to hide
   its native checkmark glyph in this Qt build, not worth fighting.

   QSpinBox only gets a background-color here (to match the white
   QLineEdit/QComboBox fields instead of Fusion's default grayish
   look) and nothing else. Isolated testing confirmed a `border`
   property on QSpinBox is specifically what breaks its native
   up/down arrow glyphs -- background-color alone does not, so it's
   safe on its own.
   ============================================================ */
QSpinBox {
    background-color: #FFFFFF;
}
QSpinBox:disabled {
    background-color: #ECECEC;
    color: #9AA0A2;
}
QComboBox QAbstractItemView {
    background-color: #FFFFFF;
    selection-background-color: #76ABAE;
    selection-color: #F5F5F5;
    outline: none;
}

/* ============================================================
   Log console (QPlainTextEdit, objectName="logConsole")
   ============================================================ */
QPlainTextEdit#logConsole {
    background-color: #303841;
    color: #F5F5F5;
    border: 1px solid #1F252B;
    border-radius: 6px;
    padding: 6px;
    selection-background-color: #76ABAE;
    selection-color: #303841;
}

/* ============================================================
   Status labels -- color stays inline/dynamic for env_status_label
   (see MainWindow._check_environment); objectName only carries
   static styling here.
   ============================================================ */
QLabel#envStatusLabel {
    font-weight: 600;
    padding: 4px 2px;
}
QLabel#summaryLabel {
    font-weight: 600;
    padding-top: 4px;
}

/* ============================================================
   Footer / status bar (author + version)
   ============================================================ */
QStatusBar {
    background-color: #F5F5F5;
    border-top: 1px solid #DADADA;
}
QStatusBar::item {
    border: none;
}
QLabel#footerLabel {
    background: transparent;
    color: #5C6670;
    font-size: 11px;
    padding: 2px 10px;
}
/* Left-hand status bar slot -- which workspace is loaded and when it was
   last written. Teal rather than the footer's grey so "your work is saved"
   reads as live state, not chrome. */
QLabel#workspaceLabel {
    background: transparent;
    color: #185861;
    font-size: 11px;
    font-weight: 600;
    padding: 2px 10px;
}
"""


class QtLogHandler(logging.Handler):
    """logging.Handler that queues formatted records for the GUI thread to drain.

    Deliberately does NOT touch any Qt widgets or emit Qt signals from the
    background thread that calls emit() (prepare_batch/run_mopac_batch/
    export_sdf all use a ThreadPoolExecutor internally, so emit() can be
    called from several worker threads at once). A QTimer on the GUI thread
    (MainWindow._drain_log_queue) polls this plain, thread-safe
    queue.Queue instead. Cross-thread Qt signal emission combined with
    background worker threads was found to intermittently crash or hang
    this app on this PySide6 build, so all direct Qt API calls are kept
    exclusively on the GUI thread.
    """

    def __init__(self, level: int = logging.INFO) -> None:
        super().__init__(level=level)
        self.queue: "queue.Queue[str]" = queue.Queue()
        self.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
        )

    def emit(self, record: logging.LogRecord) -> None:
        self.queue.put(self.format(record))


class PhaseWorker:
    """Runs one long-running mopac_batch_runner phase call on a background thread.

    Reports back purely through plain attributes polled by a QTimer on the
    GUI thread (see MainWindow._poll_worker) rather than Qt signals -- see
    QtLogHandler's docstring for why cross-thread Qt signals are avoided
    here.
    """

    def __init__(self, func, kwargs: dict) -> None:
        self._func = func
        self._kwargs = kwargs
        self.done = False
        self.result: pd.DataFrame | None = None
        self.error: str | None = None

    def run(self) -> None:
        try:
            self.result = self._func(**self._kwargs)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user via QMessageBox
            self.error = f"{type(exc).__name__}: {exc}"
        finally:
            self.done = True


_runner_logger = logging.getLogger("mopac_batch_runner")


def _export_sdf_and_copy(copy_to: Path | None, **export_kwargs) -> pd.DataFrame:
    """Run export_sdf, then optionally copy every produced .sdf file to *copy_to*.

    Executes entirely on the PhaseWorker background thread. Only touches
    plain values (Path/None) and mopac_batch_runner's existing thread-safe
    logger -- no Qt widget or Qt API call happens here, consistent with the
    threading model documented on QtLogHandler/PhaseWorker above.
    """
    result_df = export_sdf(**export_kwargs)

    if copy_to is not None and "sdf_file" in result_df.columns:
        copy_to.mkdir(parents=True, exist_ok=True)
        copied = 0
        for sdf_path in result_df["sdf_file"]:
            if pd.isna(sdf_path):
                continue
            src = Path(sdf_path)
            try:
                shutil.copy2(src, copy_to / src.name)
                copied += 1
            except OSError as exc:
                _runner_logger.warning("Could not copy %s to catalogue folder: %s", src, exc)
        _runner_logger.info("Copied %d SDF file(s) to catalogue folder: %s", copied, copy_to)

    return result_df


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("SMILE and optimize — MOPAC Batch Pipeline")
        self.resize(1500, 850)
        if ICON_PATH.exists():
            self.setWindowIcon(QIcon(str(ICON_PATH)))

        self.df: pd.DataFrame | None = None
        self._worker: PhaseWorker | None = None
        self._current_phase_name: str = ""
        self._prepare_succeeded: bool = False
        self._mopac_succeeded: bool = False
        self._phase_running: bool = False
        self._settings_locked: bool = False
        self._last_prepare_df: pd.DataFrame | None = None
        self._last_batch_df: pd.DataFrame | None = None
        self._cancel_event: threading.Event = threading.Event()

        # One timer for the lifetime of the window, started/stopped per
        # phase. It used to be constructed inside _start_phase, which left a
        # dead QTimer child parented to the window after every run.
        self._poll_timer = QTimer(self)
        self._poll_timer.setInterval(150)
        self._poll_timer.timeout.connect(self._poll_worker)

        # --- Workspace state ------------------------------------------
        # _input_path exists because until now the input file's path lived
        # only as input_path_label's text, which is not something a save
        # routine should be reading back out of a widget.
        self._workspace: dict = ws.new_workspace(PROVENANCE["version"])
        self._workspace_path: Path | None = None
        self._input_path: Path | None = None
        # Set for the duration of a restore so _on_workdir_changed's
        # "invalidate everything" behaviour cannot undo what we just applied.
        self._restoring: bool = False
        self._committing_workdir: bool = False
        # Workspaces the user has already declined to restore this session.
        self._declined_restores: set[Path] = set()
        # Cached results of the cheap */*.mop and */*.arc probes, so
        # _update_button_states (called on every keystroke) never touches
        # the filesystem itself. Refreshed by _refresh_disk_state().
        self._disk_has_mop: bool = False
        self._disk_has_arc: bool = False
        self._qsettings = QSettings(_ORG_NAME, _APP_NAME)

        self._build_ui()
        self._build_menu_bar()
        self._build_footer()
        self._update_button_states()
        self._attach_logging()
        self._check_environment()
        self._refresh_disk_state()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)

        root.addWidget(self._build_header())

        columns = QHBoxLayout()
        columns.addLayout(self._build_left_column(), stretch=2)
        columns.addLayout(self._build_right_column(), stretch=3)
        root.addLayout(columns, stretch=1)

    def _build_left_column(self) -> QVBoxLayout:
        """Configuration side: environment status, inputs, settings, retry."""
        col = QVBoxLayout()
        col.addLayout(self._build_env_row())
        col.addWidget(self._build_input_group())
        col.addWidget(self._build_workdir_group())
        col.addWidget(self._build_settings_group())
        col.addLayout(self._build_save_settings_row())
        col.addLayout(self._build_retry_row())
        col.addStretch(1)
        return col

    def _build_right_column(self) -> QVBoxLayout:
        """Results side: workflow buttons, summary, results table, log console."""
        col = QVBoxLayout()
        col.addLayout(self._build_workflow_row())

        self.summary_label = QLabel("No phase has been run yet.")
        self.summary_label.setObjectName("summaryLabel")
        col.addWidget(self.summary_label)

        self.results_table = QTableWidget()
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.setAlternatingRowColors(True)
        col.addWidget(self.results_table, stretch=3)

        self.log_view = QPlainTextEdit()
        self.log_view.setObjectName("logConsole")
        self.log_view.setReadOnly(True)
        self.log_view.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        col.addWidget(self.log_view, stretch=1)

        return col

    def _build_menu_bar(self) -> None:
        """Application menu bar: File / Run / Tools / Help.

        Attached from __init__ rather than from _build_ui() because
        QMainWindow reserves its own band for it above the central widget --
        same reason _build_footer() is separate.

        The Run menu deliberately reuses the *same slots* as the workflow
        buttons rather than duplicating any logic, and its actions are
        enabled/disabled from _update_button_states(), which stays the sole
        place that decides what is runnable.
        """
        menubar = self.menuBar()

        # -- File ------------------------------------------------------
        file_menu = menubar.addMenu("&File")

        self.act_open_ws = QAction("&Open Workspace…", self)
        self.act_open_ws.setShortcut(QKeySequence.Open)
        self.act_open_ws.triggered.connect(self._on_open_workspace)
        file_menu.addAction(self.act_open_ws)

        self.recent_menu = QMenu("&Recent Workspaces", self)
        # Rebuilt on demand: the list can be changed by another window or by
        # this one between openings, and it is cheap to regenerate.
        self.recent_menu.aboutToShow.connect(self._rebuild_recent_menu)
        file_menu.addMenu(self.recent_menu)

        self.act_save_ws = QAction("&Save Workspace", self)
        self.act_save_ws.setShortcut(QKeySequence.Save)
        self.act_save_ws.triggered.connect(self._on_save_workspace)
        file_menu.addAction(self.act_save_ws)

        self.act_save_ws_as = QAction("Save Workspace &As…", self)
        self.act_save_ws_as.setShortcut(QKeySequence.SaveAs)
        self.act_save_ws_as.triggered.connect(self._on_save_workspace_as)
        file_menu.addAction(self.act_save_ws_as)

        file_menu.addSeparator()

        self.act_load_input = QAction("&Load Input File…", self)
        self.act_load_input.setShortcut(QKeySequence("Ctrl+I"))
        self.act_load_input.triggered.connect(self._on_browse_input)
        file_menu.addAction(self.act_load_input)

        self.act_close_ws = QAction("&Close Workspace", self)
        self.act_close_ws.triggered.connect(self._on_close_workspace)
        file_menu.addAction(self.act_close_ws)

        file_menu.addSeparator()

        self.act_exit = QAction("E&xit", self)
        self.act_exit.setShortcut(QKeySequence("Ctrl+Q"))
        # Routed through close() rather than QApplication.quit() so it goes
        # through closeEvent, which cancels a running phase and arms the
        # forced-exit watchdog.
        self.act_exit.triggered.connect(self.close)
        file_menu.addAction(self.act_exit)

        # -- Run -------------------------------------------------------
        run_menu = menubar.addMenu("&Run")

        self.act_prepare = QAction("&1. Prepare (SMILES → .mop)", self)
        self.act_prepare.setShortcut(QKeySequence("Ctrl+1"))
        self.act_prepare.triggered.connect(self._on_prepare_clicked)
        run_menu.addAction(self.act_prepare)

        self.act_run_mopac = QAction("&2. Run MOPAC", self)
        self.act_run_mopac.setShortcut(QKeySequence("Ctrl+2"))
        self.act_run_mopac.triggered.connect(self._on_run_mopac_clicked)
        run_menu.addAction(self.act_run_mopac)

        self.act_export_sdf = QAction("&3. Export SDF", self)
        self.act_export_sdf.setShortcut(QKeySequence("Ctrl+3"))
        self.act_export_sdf.triggered.connect(self._on_export_sdf_clicked)
        run_menu.addAction(self.act_export_sdf)

        self.act_retry = QAction("Retry &Failed Jobs", self)
        self.act_retry.setShortcut(QKeySequence("Ctrl+4"))
        self.act_retry.triggered.connect(self._on_retry_failed_clicked)
        run_menu.addAction(self.act_retry)

        run_menu.addSeparator()

        self.act_stop = QAction("&Stop", self)
        self.act_stop.setShortcut(QKeySequence("Ctrl+."))
        self.act_stop.triggered.connect(self._on_stop_clicked)
        run_menu.addAction(self.act_stop)

        # -- Tools -----------------------------------------------------
        tools_menu = menubar.addMenu("&Tools")

        self.act_check_env = QAction("Check &Environment", self)
        self.act_check_env.triggered.connect(self._check_environment)
        tools_menu.addAction(self.act_check_env)

        self.act_scan_workdir = QAction("&Scan Working Directory", self)
        self.act_scan_workdir.triggered.connect(self._on_scan_workdir)
        tools_menu.addAction(self.act_scan_workdir)

        tools_menu.addSeparator()

        self.act_open_workdir = QAction("Open &Working Directory", self)
        self.act_open_workdir.triggered.connect(
            lambda: self._open_in_file_manager(self.workdir_edit.text())
        )
        tools_menu.addAction(self.act_open_workdir)

        self.act_open_catalogue = QAction("Open SDF &Catalogue", self)
        self.act_open_catalogue.triggered.connect(
            lambda: self._open_in_file_manager(self.sdf_catalogue_edit.text())
        )
        tools_menu.addAction(self.act_open_catalogue)

        tools_menu.addSeparator()

        self.act_clear_log = QAction("Clear &Log", self)
        self.act_clear_log.triggered.connect(self.log_view.clear)
        tools_menu.addAction(self.act_clear_log)

        # -- Help ------------------------------------------------------
        help_menu = menubar.addMenu("&Help")

        self.act_docs = QAction("&Documentation (README)", self)
        self.act_docs.triggered.connect(
            lambda: self._open_in_file_manager(README_PATH)
        )
        help_menu.addAction(self.act_docs)

        self.act_about = QAction("&About", self)
        self.act_about.triggered.connect(self._on_about)
        help_menu.addAction(self.act_about)

    def _rebuild_recent_menu(self) -> None:
        """Repopulate File ▸ Recent Workspaces from QSettings."""
        self.recent_menu.clear()
        recent = self._recent_workspaces()
        if not recent:
            placeholder = self.recent_menu.addAction("(none yet)")
            placeholder.setEnabled(False)
            return
        for path_text in recent:
            action = self.recent_menu.addAction(path_text)
            # default arg binds the current value, not the loop variable
            action.triggered.connect(
                lambda _checked=False, p=path_text: self._open_workspace(Path(p))
            )

    def _build_footer(self) -> None:
        # Left slot: live workspace state. Right slot (permanent): static
        # provenance. Keeping them on opposite ends stops the workspace text
        # from shifting the credit line around as it updates.
        self.workspace_label = QLabel("Workspace: (none)")
        self.workspace_label.setObjectName("workspaceLabel")
        self.statusBar().addWidget(self.workspace_label)

        footer_label = QLabel(
            f"Created by {PROVENANCE['author']}   |   v{PROVENANCE['version']}   |   "
            f"© {PROVENANCE['year']}   |   {PROVENANCE['contact_email']}"
        )
        footer_label.setObjectName("footerLabel")
        self.statusBar().addPermanentWidget(footer_label)

    def _build_header(self) -> QWidget:
        header = QWidget()
        header.setObjectName("headerBar")
        layout = QHBoxLayout(header)
        layout.setContentsMargins(24, 14, 24, 14)
        layout.setSpacing(16)

        if ICON_PATH.exists():
            icon_label = QLabel()
            icon_label.setObjectName("headerLogo")
            pixmap = QPixmap(str(ICON_PATH))
            icon_label.setPixmap(pixmap.scaledToHeight(68, Qt.SmoothTransformation))
            layout.addWidget(icon_label, alignment=Qt.AlignVCenter)

            divider = QFrame()
            divider.setObjectName("headerDivider")
            divider.setFixedWidth(4)
            divider.setFixedHeight(68)
            layout.addWidget(divider, alignment=Qt.AlignVCenter)

        text_col = QVBoxLayout()
        text_col.setSpacing(2)
        title = QLabel("SMILE and optimize")
        title.setObjectName("headerTitle")
        subtitle = QLabel("Chemical structure optimization via MOPAC")
        subtitle.setObjectName("headerSubtitle")
        text_col.addWidget(title)
        text_col.addWidget(subtitle)
        layout.addLayout(text_col)

        layout.addStretch(1)
        return header

    def _build_env_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        self.env_status_label = QLabel("Checking environment...")
        self.env_status_label.setObjectName("envStatusLabel")
        self.check_env_button = QPushButton("Check environment")
        self.check_env_button.clicked.connect(self._check_environment)
        row.addWidget(self.env_status_label, stretch=1)
        row.addWidget(self.check_env_button)
        return row

    def _build_input_group(self) -> QGroupBox:
        group = QGroupBox("1. Input Data")
        form = QFormLayout(group)

        self.input_path_label = QLabel("(no file loaded)")
        self.browse_input_button = QPushButton("Browse...")
        self.browse_input_button.clicked.connect(self._on_browse_input)
        input_row = QHBoxLayout()
        input_row.addWidget(self.input_path_label, stretch=1)
        input_row.addWidget(self.browse_input_button)
        form.addRow(input_row)

        self.smiles_col_combo = QComboBox()
        self.id_col_combo = QComboBox()
        self.smiles_col_combo.currentIndexChanged.connect(self._update_button_states)
        form.addRow("SMILES column:", self.smiles_col_combo)
        form.addRow("ID column:", self.id_col_combo)

        return group

    def _build_workdir_group(self) -> QGroupBox:
        group = QGroupBox("Working Directory")
        row = QHBoxLayout(group)
        self.workdir_edit = QLineEdit("mopac_jobs")
        # textChanged invalidates session state on every keystroke;
        # editingFinished is where the (comparatively expensive) filesystem
        # probe and the workspace auto-detect happen, so neither runs
        # mid-typing.
        self.workdir_edit.textChanged.connect(self._on_workdir_changed)
        self.workdir_edit.editingFinished.connect(self._on_workdir_committed)
        self.browse_workdir_button = QPushButton("Browse...")
        self.browse_workdir_button.clicked.connect(self._on_browse_workdir)
        row.addWidget(self.workdir_edit, stretch=1)
        row.addWidget(self.browse_workdir_button)
        return group

    def _build_settings_group(self) -> QGroupBox:
        group = QGroupBox("Settings")
        form = QFormLayout(group)

        self.keywords_edit = QLineEdit("PM7 PRECISE GEO-OK")
        form.addRow("MOPAC keywords:", self.keywords_edit)

        self.auto_charge_check = QCheckBox("Auto-detect charge (RDKit)")
        self.auto_charge_check.setChecked(True)
        form.addRow(self.auto_charge_check)

        self.overwrite_check = QCheckBox("Overwrite existing files")
        self.overwrite_check.setChecked(True)
        form.addRow(self.overwrite_check)

        self.prepare_workers_spin = self._make_workers_spin()
        form.addRow("Prepare workers:", self.prepare_workers_spin)

        self.mopac_workers_spin = self._make_workers_spin()
        form.addRow("MOPAC workers:", self.mopac_workers_spin)

        self.export_workers_spin = self._make_workers_spin()
        form.addRow("Export workers:", self.export_workers_spin)

        self.obabel_timeout_spin = QSpinBox()
        self.obabel_timeout_spin.setRange(1, 100000)
        self.obabel_timeout_spin.setValue(120)
        self.obabel_timeout_spin.setSuffix(" s")
        form.addRow("Prepare (obabel) timeout:", self.obabel_timeout_spin)

        self.mopac_timeout_spin = QSpinBox()
        self.mopac_timeout_spin.setRange(1, 100000)
        self.mopac_timeout_spin.setValue(900)
        self.mopac_timeout_spin.setSuffix(" s")
        form.addRow("MOPAC timeout:", self.mopac_timeout_spin)

        self.export_timeout_spin = QSpinBox()
        self.export_timeout_spin.setRange(1, 100000)
        self.export_timeout_spin.setValue(60)
        self.export_timeout_spin.setSuffix(" s")
        form.addRow("Export timeout:", self.export_timeout_spin)

        retry_separator = QLabel("Retry settings (used only by “Retry Failed Jobs”)")
        retry_separator.setStyleSheet("font-weight: 600; padding-top: 6px;")
        form.addRow(retry_separator)

        self.retry_obabel_timeout_spin = QSpinBox()
        self.retry_obabel_timeout_spin.setRange(1, 100000)
        self.retry_obabel_timeout_spin.setValue(300)
        self.retry_obabel_timeout_spin.setSuffix(" s")
        form.addRow("Retry obabel timeout:", self.retry_obabel_timeout_spin)

        self.retry_mopac_timeout_spin = QSpinBox()
        self.retry_mopac_timeout_spin.setRange(1, 100000)
        self.retry_mopac_timeout_spin.setValue(1800)
        self.retry_mopac_timeout_spin.setSuffix(" s")
        form.addRow("Retry MOPAC timeout:", self.retry_mopac_timeout_spin)

        self.retry_keywords_edit = QLineEdit()
        self.retry_keywords_edit.setPlaceholderText("(leave blank to reuse the MOPAC keywords above)")
        form.addRow("Retry MOPAC keywords:", self.retry_keywords_edit)

        self.copy_sdf_check = QCheckBox("Also copy all SDFs to a separate results folder")
        self.copy_sdf_check.setChecked(True)
        form.addRow(self.copy_sdf_check)

        self.sdf_catalogue_edit = QLineEdit("sdf_catalogue")
        self.browse_catalogue_button = QPushButton("Browse...")
        self.browse_catalogue_button.clicked.connect(self._on_browse_catalogue)
        catalogue_row = QHBoxLayout()
        catalogue_row.addWidget(self.sdf_catalogue_edit, stretch=1)
        catalogue_row.addWidget(self.browse_catalogue_button)
        form.addRow("SDF catalogue folder:", catalogue_row)

        self.copy_sdf_check.toggled.connect(self.sdf_catalogue_edit.setEnabled)
        self.copy_sdf_check.toggled.connect(self.browse_catalogue_button.setEnabled)

        return group

    def _build_save_settings_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        row.addStretch(1)
        self.save_settings_button = QPushButton("SAVE SETTINGS")
        self.save_settings_button.setObjectName("lockToggleButton")
        self.save_settings_button.clicked.connect(self._on_toggle_settings_lock_clicked)
        row.addWidget(self.save_settings_button)
        row.addStretch(1)
        return row

    @staticmethod
    def _make_workers_spin() -> QSpinBox:
        spin = QSpinBox()
        spin.setRange(1, 64)
        spin.setValue(DEFAULT_WORKERS)
        return spin

    def _build_workflow_row(self) -> QHBoxLayout:
        """The three sequential phase buttons + STOP + progress bar.

        Lives at the top of the right column, above the results table/log,
        since that's what the user is watching while a phase runs. Retry
        Failed Jobs is deliberately NOT here -- see _build_retry_row(), kept
        on the left near the retry settings it reads from.
        """
        row = QHBoxLayout()

        self.prepare_button = QPushButton("1. Prepare (SMILES -> .mop)")
        self.prepare_button.setObjectName("primaryButton")
        self.prepare_button.clicked.connect(self._on_prepare_clicked)

        self.run_mopac_button = QPushButton("2. Run MOPAC")
        self.run_mopac_button.setObjectName("primaryButton")
        self.run_mopac_button.clicked.connect(self._on_run_mopac_clicked)

        self.export_sdf_button = QPushButton("3. Export SDF")
        self.export_sdf_button.setObjectName("primaryButton")
        self.export_sdf_button.clicked.connect(self._on_export_sdf_clicked)

        self.stop_button = QPushButton("STOP")
        self.stop_button.setObjectName("stopButton")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._on_stop_clicked)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(False)

        row.addWidget(self.prepare_button)
        row.addWidget(self.run_mopac_button)
        row.addWidget(self.export_sdf_button)
        row.addWidget(self.stop_button)
        row.addWidget(self.progress_bar, stretch=1)
        return row

    def _build_retry_row(self) -> QHBoxLayout:
        """Retry Failed Jobs -- kept away from the main workflow row on
        purpose, and placed right under Settings since it reads the retry
        timeout/keyword fields defined there.
        """
        row = QHBoxLayout()
        row.addStretch(1)
        self.retry_failed_button = QPushButton("Retry Failed Jobs")
        self.retry_failed_button.setObjectName("primaryButton")
        self.retry_failed_button.clicked.connect(self._on_retry_failed_clicked)
        row.addWidget(self.retry_failed_button)
        row.addStretch(1)
        return row

    def _attach_logging(self) -> None:
        self._log_handler = QtLogHandler(level=logging.INFO)

        runner_logger = logging.getLogger("mopac_batch_runner")
        runner_logger.addHandler(self._log_handler)
        runner_logger.setLevel(logging.INFO)
        runner_logger.propagate = False

        self._log_timer = QTimer(self)
        self._log_timer.setInterval(150)
        self._log_timer.timeout.connect(self._drain_log_queue)
        self._log_timer.start()

    def _drain_log_queue(self) -> None:
        while True:
            try:
                line = self._log_handler.queue.get_nowait()
            except queue.Empty:
                break
            self.log_view.appendPlainText(line)

    # ------------------------------------------------------------------
    # Environment check
    # ------------------------------------------------------------------
    def _check_environment(self) -> None:
        try:
            validate_environment()
        except EnvironmentError as exc:
            self.env_status_label.setText(f"Environment problem: {exc}")
            self.env_status_label.setStyleSheet("color: #b00020;")
        else:
            self.env_status_label.setText("Environment OK - obabel and mopac found in PATH.")
            self.env_status_label.setStyleSheet("color: #1b7a1b;")

    # ------------------------------------------------------------------
    # Input file / workdir selection
    # ------------------------------------------------------------------
    def _on_browse_input(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select input file",
            str(self._input_path.parent) if self._input_path else "",
            "Excel/CSV Files (*.xlsx *.csv);;Excel Files (*.xlsx);;"
            "CSV Files (*.csv);;All Files (*)",
        )
        if path:
            self._load_input_file(Path(path))

    def _load_input_file(self, path: Path, reset_phases: bool = True) -> bool:
        """Read an .xlsx/.csv into self.df and repopulate the column combos.

        Split out of _on_browse_input so the browse button, File ▸ Load Input
        File… and workspace restore all go through exactly one implementation
        of "load this spreadsheet".

        reset_phases=False is used by restore, which sets the phase flags
        itself from the workspace and the results CSVs -- clearing them here
        would immediately undo that.

        Returns True when the file loaded.
        """
        try:
            if path.suffix.lower() == ".csv":
                df = pd.read_csv(path)
            else:
                df = pd.read_excel(path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Failed to load file", str(exc))
            return False

        self.df = df
        self._input_path = path
        self.input_path_label.setText(str(path))
        if reset_phases:
            self._prepare_succeeded = False
            self._mopac_succeeded = False
            self._last_prepare_df = None
            self._last_batch_df = None

        self.smiles_col_combo.clear()
        self.smiles_col_combo.addItems([str(c) for c in df.columns])

        self.id_col_combo.clear()
        self.id_col_combo.addItem("(none - auto-generate sequential IDs)", None)
        for col in df.columns:
            self.id_col_combo.addItem(str(col), str(col))

        self._update_button_states()
        return True

    def _on_browse_workdir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select working directory")
        if path:
            self.workdir_edit.setText(path)
            self._on_workdir_committed()

    def _on_workdir_changed(self, _text: str) -> None:
        """Invalidate in-session phase progress when the workdir is retyped.

        Fires on every keystroke, so it stays deliberately cheap and does not
        touch the filesystem -- the disk probe lives in _on_workdir_committed,
        wired to editingFinished instead.
        """
        if self._restoring:
            return
        self._prepare_succeeded = False
        self._mopac_succeeded = False
        self._last_prepare_df = None
        self._last_batch_df = None
        self._disk_has_mop = False
        self._disk_has_arc = False
        self._update_button_states()

    def _on_workdir_committed(self) -> None:
        """The user finished typing/picking a working directory.

        Re-probes what is on disk (which is what re-enables the later phase
        buttons for a folder prepared in an earlier session) and offers to
        restore a workspace file if the folder has one.

        Guarded against re-entry: editingFinished fires on focus-out, and
        opening a modal from inside it can move focus enough to emit the
        signal a second time -- which would stack a duplicate dialog.
        """
        if self._restoring or self._committing_workdir:
            return
        self._committing_workdir = True
        try:
            self._refresh_disk_state()
            self._update_button_states()
            self._maybe_offer_restore()
        finally:
            self._committing_workdir = False

    def _on_browse_catalogue(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select SDF catalogue folder")
        if path:
            self.sdf_catalogue_edit.setText(path)

    # ------------------------------------------------------------------
    # Settings lock
    # ------------------------------------------------------------------
    def _on_toggle_settings_lock_clicked(self) -> None:
        """Lock the Settings section so its values can't be changed by accident.

        Toggles between "SAVE SETTINGS" (editable) and "EDIT SETTINGS"
        (locked). sdf_catalogue_edit/browse_catalogue_button are also
        conditionally enabled by copy_sdf_check itself, so on unlock they
        only re-enable if that checkbox is actually checked, instead of
        blindly overriding that dependency.
        """
        lock = not self._settings_locked
        self._settings_locked = lock

        for attr in _LOCKABLE_SETTINGS:
            getattr(self, attr).setEnabled(not lock)

        if lock:
            self.sdf_catalogue_edit.setEnabled(False)
            self.browse_catalogue_button.setEnabled(False)
        else:
            self.sdf_catalogue_edit.setEnabled(self.copy_sdf_check.isChecked())
            self.browse_catalogue_button.setEnabled(self.copy_sdf_check.isChecked())

        if lock:
            self.save_settings_button.setText("EDIT SETTINGS")
            self.save_settings_button.setStyleSheet(
                "QPushButton#lockToggleButton {"
                "   background-color: #FF5722;"
                "   color: #F5F5F5;"
                "   border: 2px solid #FF5722;"
                "}"
                "QPushButton#lockToggleButton:hover {"
                "   background-color: #E64A19;"
                "   border: 2px solid #E64A19;"
                "}"
                "QPushButton#lockToggleButton:pressed {"
                "   background-color: #C43D14;"
                "}"
            )
        else:
            self.save_settings_button.setText("SAVE SETTINGS")
            self.save_settings_button.setStyleSheet("")

        self._update_button_states()
        if lock:
            # Confirming the settings is the natural moment to persist them:
            # from here on they are what the phases will actually run with.
            self._save_workspace("settings confirmed")

    # ------------------------------------------------------------------
    # Workspace: collecting and applying state
    # ------------------------------------------------------------------
    def _collect_settings(self) -> dict:
        """Read the Settings widgets into plain JSON-able values."""
        values: dict = {}
        for key, attr, kind in _SETTINGS_FIELDS:
            widget = getattr(self, attr)
            if kind == "text":
                values[key] = widget.text()
            elif kind == "bool":
                values[key] = widget.isChecked()
            else:
                values[key] = widget.value()
        return values

    def _apply_settings(self, values: dict) -> None:
        """Push saved values back into the Settings widgets.

        Signals are blocked per widget so that restoring copy_sdf triggers
        neither of the two setEnabled connections wired up in
        _build_settings_group -- the lock state is applied separately, and
        afterwards, by _apply_workspace.
        """
        for key, attr, kind in _SETTINGS_FIELDS:
            if key not in values:
                continue
            widget = getattr(self, attr)
            blocked = widget.blockSignals(True)
            try:
                if kind == "text":
                    widget.setText(str(values[key]))
                elif kind == "bool":
                    widget.setChecked(bool(values[key]))
                else:
                    widget.setValue(int(values[key]))
            except (TypeError, ValueError):
                _runner_logger.warning(
                    "Ignoring unusable saved value for %r: %r", key, values[key]
                )
            finally:
                widget.blockSignals(blocked)

    def _set_settings_locked(self, locked: bool) -> None:
        """Force the lock into a given state, reusing the toggle's own logic."""
        if self._settings_locked != locked:
            self._on_toggle_settings_lock_clicked()

    def _collect_workspace(self) -> dict:
        """Snapshot everything worth persisting, on top of the current one.

        Built by updating self._workspace rather than from scratch so the
        per-phase history (states, timestamps, counts) accumulated across the
        session survives each save.
        """
        data = self._workspace
        data["app_version"] = PROVENANCE["version"]
        data["workdir"] = self.workdir_edit.text()
        data["settings"] = self._collect_settings()
        data["settings_locked"] = self._settings_locked

        source = dict(data.get("input") or {})
        if self._input_path is not None:
            source.update(ws.fingerprint_input(self._input_path))
            source["path"] = str(self._input_path)
            source["n_rows"] = int(len(self.df)) if self.df is not None else None
            source["columns"] = (
                [str(c) for c in self.df.columns] if self.df is not None else []
            )
        source["smiles_col"] = self.smiles_col_combo.currentText()
        source["id_col"] = self.id_col_combo.currentData()
        data["input"] = source
        return data

    def _record_phase(self, phase_name: str, state: str, df: pd.DataFrame | None = None) -> None:
        """Update one phase's entry in the in-memory workspace.

        Called at phase start with state="running" and at phase end with the
        outcome. The "running" marker is the whole point of saving at the
        start: a process killed 90 minutes into a two-hour Prepare reads back
        as interrupted rather than as never having run.
        """
        key = _PHASE_TO_WORKSPACE_KEY.get(phase_name)
        if key is None:
            return
        phase = self._workspace["phases"].setdefault(key, ws.new_phase(key))
        phase["state"] = state
        if state == ws.STATE_RUNNING:
            phase["started_at"] = ws.now_iso()
            phase["finished_at"] = None
        else:
            phase["finished_at"] = ws.now_iso()
        if df is not None and "status" in df.columns:
            phase["counts"] = {
                str(status): int(count) for status, count in df["status"].value_counts().items()
            }

    # ------------------------------------------------------------------
    # Workspace: saving
    # ------------------------------------------------------------------
    def _default_workspace_path(self) -> Path | None:
        text = self.workdir_edit.text().strip()
        return ws.workspace_path_for(text) if text else None

    def _save_workspace(self, reason: str, path: Path | None = None) -> bool:
        """Write the workspace, silently by default.

        Autosave must never be able to interrupt the user, so a failure here
        is a warning in the log console (and a note in the status bar) rather
        than a modal -- losing the sidecar is recoverable, since the working
        directory still holds every .mop/.out/.arc and the results CSVs.

        Restoring is never a save: applying a workspace re-locks the Settings
        section, and the lock autosaves, which would otherwise write the
        half-applied state straight back over the file being read.
        """
        if self._restoring:
            return False
        target = path or self._workspace_path or self._default_workspace_path()
        if target is None:
            return False
        try:
            ws.save_workspace(target, self._collect_workspace())
        except OSError as exc:
            _runner_logger.warning("Could not save workspace to %s: %s", target, exc)
            self.workspace_label.setText(f"Workspace: NOT saved — {exc.strerror or exc}")
            return False

        self._workspace_path = target
        self._remember_recent(target)
        self.workspace_label.setText(
            f"Workspace: {target.parent.name or target.parent} — "
            f"saved {ws.now_iso()[11:16]}"
        )
        _runner_logger.info("Workspace saved (%s): %s", reason, target)
        return True

    def _on_save_workspace(self) -> None:
        if not self.workdir_edit.text().strip():
            QMessageBox.warning(
                self, "No working directory",
                "Set a working directory first — the workspace file is saved inside it.",
            )
            return
        self._save_workspace("manual save")

    def _on_save_workspace_as(self) -> None:
        """Export a copy of the workspace somewhere else.

        Deliberately does NOT change where autosave writes: the canonical
        copy stays at <workdir>/workspace.json so the working directory
        remains self-describing. This is for sharing or archiving a setup.
        """
        start = self._default_workspace_path() or Path(ws.WORKSPACE_FILENAME)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save workspace copy as", str(start), "Workspace JSON (*.json);;All Files (*)"
        )
        if not path:
            return
        try:
            ws.save_workspace(Path(path), self._collect_workspace())
        except OSError as exc:
            QMessageBox.critical(self, "Could not save workspace", str(exc))
            return
        _runner_logger.info("Workspace copy written to %s", path)

    # ------------------------------------------------------------------
    # Workspace: opening and restoring
    # ------------------------------------------------------------------
    def _on_open_workspace(self) -> None:
        start = self.workdir_edit.text().strip() or ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open workspace", start, "Workspace JSON (*.json);;All Files (*)"
        )
        if path:
            self._open_workspace(Path(path))

    def _on_close_workspace(self) -> None:
        """Detach from the current workspace without touching the disk."""
        self._workspace = ws.new_workspace(PROVENANCE["version"])
        self._workspace_path = None
        self.workspace_label.setText("Workspace: (none)")
        _runner_logger.info("Workspace closed (nothing was deleted).")

    def _maybe_offer_restore(self) -> None:
        """Offer to restore if the chosen working directory already has one."""
        candidate = self._default_workspace_path()
        if candidate is None or not candidate.is_file():
            return
        if self._workspace_path is not None and candidate == self._workspace_path:
            return
        # Saying no once means no: without this, simply clicking back into
        # the working-directory field would re-ask every time.
        if candidate in self._declined_restores:
            return
        try:
            data = ws.load_workspace(candidate)
        except ws.WorkspaceError as exc:
            _runner_logger.warning("Ignoring unreadable workspace at %s: %s", candidate, exc)
            return

        answer = QMessageBox.question(
            self,
            "Saved workspace found",
            f"{candidate} describes an earlier session.\n\n"
            f"{ws.describe_workspace(data)}\n\nRestore it?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        if answer == QMessageBox.Yes:
            self._open_workspace(candidate, preloaded=data)
        else:
            self._declined_restores.add(candidate)

    def _open_workspace(self, path: Path, preloaded: dict | None = None) -> None:
        """Load a workspace file and put the window back into that state."""
        if self._phase_running:
            QMessageBox.warning(
                self, "Phase running",
                "Stop the running phase before opening another workspace.",
            )
            return
        try:
            data = preloaded if preloaded is not None else ws.load_workspace(path)
        except ws.WorkspaceError as exc:
            QMessageBox.critical(self, "Could not open workspace", str(exc))
            return

        if ws.is_future_version(data):
            QMessageBox.warning(
                self, "Newer workspace",
                f"{path.name} was written by a newer version of this app "
                f"(schema {data.get('schema_version')} vs {ws.SCHEMA_VERSION}). "
                "Restoring it anyway; some settings may be ignored.",
            )

        resume = self._ask_resume_mode(data)
        if resume is None:
            return

        self._restoring = True
        try:
            self._apply_workspace(data, resume=resume)
        finally:
            self._restoring = False

        self._workspace = data
        self._workspace_path = path
        self._remember_recent(path)
        self._refresh_disk_state()
        self._update_button_states()
        self.workspace_label.setText(
            f"Workspace: {path.parent.name or path.parent} — restored"
        )
        _runner_logger.info("Workspace restored from %s", path)
        self._start_scan(silent=True)

    def _ask_resume_mode(self, data: dict) -> bool | None:
        """Ask whether to resume (skip finished work) or restore verbatim.

        This exists because of a genuine footgun: overwrite defaults to True
        in this GUI while the library default is False. Restoring an
        interrupted 20,000-compound Prepare and pressing Prepare again with
        overwrite on would redo every one of the two hours already spent.
        Resuming turns it off so finished compounds are skipped instead.

        Returns True to resume, False to restore settings as saved, or None
        if the user cancelled.
        """
        if not ws.has_unfinished_work(data):
            return False

        box = QMessageBox(self)
        box.setIcon(QMessageBox.Question)
        box.setWindowTitle("Resume this workspace?")
        box.setText("This workspace has unfinished or failed work.")
        box.setInformativeText(
            f"{ws.describe_workspace(data)}\n\n"
            "“Resume” turns Overwrite off, so compounds that already have "
            "their .mop/.out/.sdf files are skipped instead of being redone."
        )
        resume_button = box.addButton("Resume — skip completed", QMessageBox.AcceptRole)
        as_saved_button = box.addButton("Restore settings as saved", QMessageBox.AcceptRole)
        box.addButton(QMessageBox.Cancel)
        box.setDefaultButton(resume_button)
        box.exec()

        clicked = box.clickedButton()
        if clicked is resume_button:
            return True
        if clicked is as_saved_button:
            return False
        return None

    def _apply_workspace(self, data: dict, resume: bool) -> None:
        """Put every restorable piece of *data* back into the window.

        Runs with self._restoring set, so the workdir field's per-keystroke
        invalidation cannot wipe the phase flags this sets.
        """
        self.workdir_edit.setText(str(data.get("workdir") or ""))
        self._apply_settings(data.get("settings") or {})
        if resume:
            self.overwrite_check.setChecked(False)
            _runner_logger.info(
                "Resuming: Overwrite turned off so finished compounds are skipped."
            )
        # Re-apply copy_sdf's dependency by hand -- _apply_settings blocked
        # the signal that normally does this.
        catalogue_enabled = self.copy_sdf_check.isChecked()
        self.sdf_catalogue_edit.setEnabled(catalogue_enabled)
        self.browse_catalogue_button.setEnabled(catalogue_enabled)

        self._restore_input(data)
        self._set_settings_locked(bool(data.get("settings_locked")))
        self._restore_phase_results(data)

    def _restore_input(self, data: dict) -> None:
        """Reload the recorded input file and re-select its columns.

        Three outcomes worth calling out:
          * file present and unchanged  -> silent, everything works;
          * file moved/renamed          -> offer to locate it;
          * file edited since the save  -> warn, because compound IDs are
            derived from the whole DataFrame in order (_build_job_list in
            mopac_batch_runner.py) and retry_failed_batch re-derives them, so
            a reordered sheet maps compounds onto the wrong directories.

        Declining to supply the file is a supported outcome, not a failure:
        Run MOPAC and Export SDF read the working directory alone, so they
        stay usable. Only Prepare and Retry need the spreadsheet.
        """
        source = data.get("input") or {}
        path_text = source.get("path")
        if not path_text:
            return

        path = Path(path_text)
        if not path.is_file():
            answer = QMessageBox.question(
                self,
                "Input file not found",
                f"The workspace refers to:\n{path}\n\nIt is no longer there. "
                "Locate it now?\n\n(Without it you can still run MOPAC and "
                "export SDFs — only Prepare and Retry need the spreadsheet.)",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes,
            )
            if answer != QMessageBox.Yes:
                return
            picked, _ = QFileDialog.getOpenFileName(
                self, "Locate input file", str(path.parent),
                "Excel/CSV Files (*.xlsx *.csv);;All Files (*)",
            )
            if not picked:
                return
            path = Path(picked)

        if not self._load_input_file(path, reset_phases=False):
            return

        if not ws.input_matches(source, ws.fingerprint_input(path)):
            QMessageBox.warning(
                self,
                "Input file has changed",
                f"{path.name} is not the file this workspace was saved with.\n\n"
                "Compound IDs are rebuilt from the whole spreadsheet in row "
                "order, so if rows were added, removed or reordered, Retry "
                "Failed Jobs can match compounds to the wrong job folders.\n\n"
                "Preparing into a fresh working directory is the safe option.",
            )

        smiles_col = source.get("smiles_col")
        if smiles_col:
            index = self.smiles_col_combo.findText(str(smiles_col))
            if index >= 0:
                self.smiles_col_combo.setCurrentIndex(index)
        id_index = self.id_col_combo.findData(source.get("id_col"))
        if id_index >= 0:
            self.id_col_combo.setCurrentIndex(id_index)

    def _restore_phase_results(self, data: dict) -> None:
        """Reload the per-phase results CSVs and re-enable the phases they prove.

        The CSVs are the authority on *what happened*; the workspace is the
        authority on *what was configured*. A phase counts as having run if
        its CSV is on disk, whatever the workspace claims -- that way a
        directory processed through the Python API still restores usefully.
        """
        workdir = Path(self.workdir_edit.text().strip() or ".")
        loaded: list[tuple[str, pd.DataFrame]] = []

        for key, label in (("prepare", "Prepare"), ("mopac", "Run MOPAC"), ("export", "Export SDF")):
            csv_path = workdir / ws.PHASE_RESULT_CSV[key]
            if not csv_path.is_file():
                continue
            try:
                df = pd.read_csv(csv_path)
            except (OSError, ValueError) as exc:
                _runner_logger.warning("Could not read %s: %s", csv_path, exc)
                continue
            loaded.append((label, df))
            if key == "prepare":
                self._prepare_succeeded = True
                self._last_prepare_df = df
            elif key == "mopac":
                self._mopac_succeeded = True
                self._last_batch_df = df

        if loaded:
            label, df = loaded[-1]
            self._populate_table(df)
            self._update_summary_label(f"{label} (restored)", df)
        else:
            self.summary_label.setText(
                "Workspace restored — no results CSV found in the working directory yet."
            )

    # ------------------------------------------------------------------
    # Workspace: recent list
    # ------------------------------------------------------------------
    def _recent_workspaces(self) -> list[str]:
        stored = self._qsettings.value(_RECENT_KEY, [])
        if isinstance(stored, str):  # QSettings collapses one-element lists
            stored = [stored]
        return [p for p in (stored or []) if p]

    def _remember_recent(self, path: Path) -> None:
        text = str(path)
        recent = [p for p in self._recent_workspaces() if p != text]
        recent.insert(0, text)
        self._qsettings.setValue(_RECENT_KEY, recent[:_MAX_RECENT])

    # ------------------------------------------------------------------
    # Working-directory inspection
    # ------------------------------------------------------------------
    def _refresh_disk_state(self) -> None:
        """Probe whether the working directory already holds .mop / .arc files.

        This is what lets a folder prepared in an earlier session (or by the
        notebook, or by the Python API) unlock Run MOPAC and Export SDF
        without re-running Prepare -- the library has always supported it,
        since run_mopac_batch globs */*.mop and export_sdf globs */*.arc.

        next(glob) short-circuits on the first match, so this stays cheap
        even against a 20,000-folder working directory.
        """
        text = self.workdir_edit.text().strip()
        workdir = Path(text) if text else None
        if workdir is None or not workdir.is_dir():
            self._disk_has_mop = False
            self._disk_has_arc = False
            return
        try:
            self._disk_has_mop = next(workdir.glob("*/*.mop"), None) is not None
            self._disk_has_arc = next(workdir.glob("*/*.arc"), None) is not None
        except OSError:
            self._disk_has_mop = False
            self._disk_has_arc = False

    def _on_scan_workdir(self) -> None:
        self._start_scan(silent=False)

    def _start_scan(self, silent: bool) -> None:
        """Count the artefacts actually on disk, on a background thread.

        Worth doing on restore rather than trusting the results CSVs: those
        are only rewritten every 200 completed jobs (_CHECKPOINT_EVERY in
        mopac_batch_runner.py), so after a hard exit they can be missing up
        to 199 compounds whose files are sitting in the working directory.

        Runs as a pseudo-phase through the existing _start_phase machinery so
        it inherits the progress bar, the STOP button and the button gating,
        instead of introducing a second background-work pattern.
        """
        if self._phase_running:
            if not silent:
                QMessageBox.information(
                    self, "Busy", "Wait for the current phase to finish, then scan."
                )
            return
        text = self.workdir_edit.text().strip()
        if not text or not Path(text).is_dir():
            if not silent:
                QMessageBox.warning(
                    self, "No working directory", f"{text or '(empty)'} is not a directory."
                )
            return
        self._start_phase(ws.scan_workdir, {"workdir": text}, _SCAN_PHASE_NAME)

    def _apply_scan_result(self, scan: dict) -> None:
        self._workspace["workdir_scan"] = scan
        self.summary_label.setText(ws.describe_scan(scan))
        _runner_logger.info(ws.describe_scan(scan))

    # ------------------------------------------------------------------
    # Misc menu actions
    # ------------------------------------------------------------------
    @staticmethod
    def _open_in_file_manager(target: "str | Path") -> None:
        path = Path(str(target).strip())
        if not str(path) or not path.exists():
            _runner_logger.warning("Cannot open %s — it does not exist yet.", path)
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path.resolve())))

    def _on_about(self) -> None:
        QMessageBox.about(
            self,
            "About SMILE and optimize",
            f"<b>SMILE and optimize</b> v{PROVENANCE['version']}<br><br>"
            "Chemical structure optimization via MOPAC.<br>"
            "SMILES → 3D (Open Babel) → MOPAC → SDF.<br><br>"
            f"Created by {PROVENANCE['author']}<br>"
            f"© {PROVENANCE['year']} — {PROVENANCE['contact_email']}",
        )

    # ------------------------------------------------------------------
    # Phase actions
    # ------------------------------------------------------------------
    def _on_prepare_clicked(self) -> None:
        if self.df is None:
            QMessageBox.warning(self, "No input file", "Load an input file first.")
            return
        smiles_col = self.smiles_col_combo.currentText()
        if not smiles_col:
            QMessageBox.warning(self, "No SMILES column", "Choose a SMILES column first.")
            return
        if not self._settings_locked:
            QMessageBox.warning(
                self, "Settings not saved", "Click SAVE SETTINGS to confirm settings first."
            )
            return
        id_col = self.id_col_combo.currentData()

        kwargs = dict(
            df=self.df,
            smiles_col=smiles_col,
            id_col=id_col,
            keywords=self.keywords_edit.text(),
            workdir=self.workdir_edit.text(),
            auto_charge=self.auto_charge_check.isChecked(),
            overwrite=self.overwrite_check.isChecked(),
            max_workers=self.prepare_workers_spin.value(),
            obabel_timeout=self.obabel_timeout_spin.value(),
        )
        self._start_phase(prepare_batch, kwargs, "Prepare")

    def _on_run_mopac_clicked(self) -> None:
        kwargs = dict(
            workdir=self.workdir_edit.text(),
            overwrite=self.overwrite_check.isChecked(),
            timeout=self.mopac_timeout_spin.value(),
            max_workers=self.mopac_workers_spin.value(),
        )
        self._start_phase(run_mopac_batch, kwargs, "Run MOPAC")

    def _on_export_sdf_clicked(self) -> None:
        copy_to: Path | None = None
        if self.copy_sdf_check.isChecked():
            catalogue_text = self.sdf_catalogue_edit.text().strip()
            if catalogue_text:
                copy_to = Path(catalogue_text)

        kwargs = dict(
            workdir=self.workdir_edit.text(),
            overwrite=self.overwrite_check.isChecked(),
            timeout=self.export_timeout_spin.value(),
            max_workers=self.export_workers_spin.value(),
        )
        export_func = functools.partial(_export_sdf_and_copy, copy_to)
        self._start_phase(export_func, kwargs, "Export SDF")

    def _on_retry_failed_clicked(self) -> None:
        if self.df is None:
            QMessageBox.warning(self, "No input file", "Load an input file first.")
            return
        smiles_col = self.smiles_col_combo.currentText()
        if not smiles_col:
            QMessageBox.warning(self, "No SMILES column", "Choose a SMILES column first.")
            return
        id_col = self.id_col_combo.currentData()
        retry_keywords = self.retry_keywords_edit.text().strip() or None

        kwargs = dict(
            df=self.df,
            smiles_col=smiles_col,
            id_col=id_col,
            keywords=self.keywords_edit.text(),
            workdir=self.workdir_edit.text(),
            auto_charge=self.auto_charge_check.isChecked(),
            obabel_timeout=self.retry_obabel_timeout_spin.value(),
            mopac_timeout=self.retry_mopac_timeout_spin.value(),
            retry_keywords=retry_keywords,
            max_workers=self.mopac_workers_spin.value(),
        )
        self._start_phase(retry_failed_batch, kwargs, "Retry Failed")

    def _on_stop_clicked(self) -> None:
        """Ask the running phase to stop immediately.

        Sets the shared cancel_event, which every mopac_batch_runner batch
        function polls -- any obabel/MOPAC process currently running is
        killed within about one poll interval (see run_command's docstring
        in mopac_batch_runner.py), and no further compounds are started.
        The phase's own worker thread notices it finished (whether by
        completion or cancellation) via the normal _poll_worker path.
        """
        if not self._phase_running:
            return
        self._cancel_event.set()
        self.stop_button.setEnabled(False)
        self.stop_button.setText("Stopping...")
        _runner_logger.info("Stop requested by user -- killing the current job(s).")

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt override
        """Make closing the window actually end the process.

        The PhaseWorker thread is a daemon, but the ThreadPoolExecutor that
        the batch functions run inside is not: concurrent.futures registers an
        atexit hook that joins its worker threads, so a single worker stuck on
        an unresponsive external process kept the whole app alive after the
        window was gone -- it had to be killed from the terminal.

        So: ask nicely first, then leave via os._exit, which skips the atexit
        machinery. Results are safe to abandon this way because every batch
        function checkpoints its CSV as jobs complete rather than only at the
        end (see _write_partial_csv in mopac_batch_runner.py).
        """
        if self._phase_running:
            self._cancel_event.set()
            _runner_logger.info("Window closed -- stopping the running phase.")
            # Whatever was running did not finish; record that before the
            # snapshot, so a restore reports it as interrupted.
            self._record_phase(self._current_phase_name, ws.STATE_STOPPED)

        # A few kilobytes of JSON, written atomically -- cheap enough to do
        # on the way out even when a phase is being force-killed behind us.
        if self.workdir_edit.text().strip():
            self._save_workspace("window closed")

        event.accept()
        super().closeEvent(event)

        if self._phase_running:
            # A plain daemon timer rather than QTimer.singleShot: closing the
            # last window stops the Qt event loop, so a queued Qt timer would
            # never fire -- and the interpreter would then block in the very
            # atexit join this is here to escape.
            watchdog = threading.Timer(_FORCED_EXIT_GRACE_S, self._force_exit)
            watchdog.daemon = True
            watchdog.start()

    @staticmethod
    def _force_exit() -> None:
        _runner_logger.warning(
            "Worker threads did not stop within %.0fs of closing the window; "
            "exiting anyway.", _FORCED_EXIT_GRACE_S,
        )
        logging.shutdown()
        os._exit(0)

    # ------------------------------------------------------------------
    # Button gating
    # ------------------------------------------------------------------
    def _set_phase_buttons_enabled(self, running: bool) -> None:
        """running=True while a phase is active; False once it finishes."""
        self._phase_running = running
        if not running:
            self.stop_button.setText("STOP")
        self._update_button_states()

    @staticmethod
    def _has_failures(df: pd.DataFrame | None) -> bool:
        if df is None or df.empty or "status" not in df.columns:
            return False
        return bool(df["status"].isin(["failed", "timeout", "cancelled"]).any())

    def _update_button_states(self) -> None:
        """Sole place that computes final enabled/disabled for the phase buttons.

        Phases 2 and 3 unlock either from having run this session or from
        what is already in the working directory (see _refresh_disk_state).
        The library never needed the earlier phase to have run in-process --
        run_mopac_batch discovers its work by globbing */*.mop and export_sdf
        by globbing */*.arc -- so a directory left over from a previous
        session, the notebook, or the Python API is picked up as-is.
        """
        can_prepare = (
            self.df is not None
            and bool(self.smiles_col_combo.currentText())
            and self._settings_locked
        )
        self.prepare_button.setEnabled(can_prepare and not self._phase_running)
        self.run_mopac_button.setEnabled(
            (self._prepare_succeeded or self._disk_has_mop) and not self._phase_running
        )
        self.export_sdf_button.setEnabled(
            (self._mopac_succeeded or self._disk_has_arc) and not self._phase_running
        )
        self.stop_button.setEnabled(self._phase_running)

        # Retry rebuilds compound IDs from the input DataFrame, so unlike the
        # two phases above it genuinely cannot run without the spreadsheet.
        can_retry = (
            self.df is not None
            and (self._mopac_succeeded or self._disk_has_mop)
            and (
                self._has_failures(self._last_prepare_df)
                or self._has_failures(self._last_batch_df)
            )
        )
        self.retry_failed_button.setEnabled(can_retry and not self._phase_running)

        # Mirror onto the Run menu so the menu items and the buttons can
        # never disagree about what is runnable.
        for action, button in (
            (self.act_prepare, self.prepare_button),
            (self.act_run_mopac, self.run_mopac_button),
            (self.act_export_sdf, self.export_sdf_button),
            (self.act_retry, self.retry_failed_button),
            (self.act_stop, self.stop_button),
        ):
            action.setEnabled(button.isEnabled())
        self.act_scan_workdir.setEnabled(not self._phase_running)

    # ------------------------------------------------------------------
    # Threaded phase execution
    # ------------------------------------------------------------------
    def _start_phase(self, func, kwargs: dict, phase_name: str) -> None:
        self._cancel_event.clear()
        kwargs = {**kwargs, "cancel_event": self._cancel_event}

        self._set_phase_buttons_enabled(True)
        self.progress_bar.setVisible(True)
        self._current_phase_name = phase_name

        # Mark the phase as running *before* it starts, and persist that.
        # This is what makes a 2-hour Prepare recoverable after a crash or a
        # forced exit: the workspace on disk says "running", so the next
        # session restores it as interrupted rather than as never-run.
        if phase_name != _SCAN_PHASE_NAME:
            self._record_phase(phase_name, ws.STATE_RUNNING)
            self._save_workspace(f"{phase_name} started")

        # PhaseWorker.run() executes on a plain background thread and never
        # touches Qt directly; a QTimer on the GUI thread (_poll_worker)
        # checks worker.done instead of the worker emitting a signal back.
        # An earlier QThread-based version, and later a raw-thread-plus-
        # Qt-signal version, both intermittently crashed or hung the whole
        # process under concurrent ThreadPoolExecutor activity inside
        # run_mopac_batch -- polling sidesteps cross-thread Qt calls
        # entirely and was the version that held up under repeated testing.
        self._worker = PhaseWorker(func, kwargs)
        threading.Thread(target=self._worker.run, daemon=True).start()
        self._poll_timer.start()

    def _poll_worker(self) -> None:
        worker = self._worker
        if worker is None or not worker.done:
            return

        self._poll_timer.stop()
        self._worker = None
        self.progress_bar.setVisible(False)
        phase_name = self._current_phase_name
        stopped = self._cancel_event.is_set()

        if worker.error is not None:
            QMessageBox.critical(self, f"{phase_name} failed", worker.error)
        elif phase_name == _SCAN_PHASE_NAME:
            # The scan is not a pipeline phase: it returns a counts dict, not
            # a results DataFrame, so it never touches the results table.
            self._apply_scan_result(worker.result or {})
        else:
            # "Succeeded" means the phase call itself didn't raise -- not
            # that every row succeeded. Partial per-row failure is normal
            # pipeline behavior and must not block progressing to the next
            # phase for the rows that did succeed.
            if phase_name == "Prepare":
                self._prepare_succeeded = True
                self._last_prepare_df = worker.result
            elif phase_name == "Run MOPAC":
                self._mopac_succeeded = True
                self._last_batch_df = worker.result
            elif phase_name == "Retry Failed":
                self._last_batch_df = worker.result
            self._populate_table(worker.result)
            self._update_summary_label(phase_name, worker.result)
            if stopped:
                self.summary_label.setText(self.summary_label.text() + "  (stopped by user)")

        if phase_name != _SCAN_PHASE_NAME:
            # Anything short of "ran to the end without being cancelled"
            # counts as unfinished, so a restore offers to resume it.
            finished_cleanly = worker.error is None and not stopped
            self._record_phase(
                phase_name,
                ws.STATE_COMPLETED if finished_cleanly else ws.STATE_STOPPED,
                worker.result,
            )

        # A finished phase has usually changed what is on disk -- Prepare
        # creates the .mop files that unlock Run MOPAC, MOPAC creates the
        # .arc files that unlock Export SDF.
        self._refresh_disk_state()
        self._save_workspace(f"{phase_name} finished")
        self._set_phase_buttons_enabled(False)

    # ------------------------------------------------------------------
    # Results display
    # ------------------------------------------------------------------
    def _populate_table(self, df: pd.DataFrame) -> None:
        self.results_table.clear()
        self.results_table.setColumnCount(len(df.columns))
        self.results_table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        self.results_table.setRowCount(len(df))

        status_idx = list(df.columns).index("status") if "status" in df.columns else None
        fail_statuses = {"failed", "timeout"}
        cancelled_color = QColor("#9AA0A2")  # matches disabled-widget grey elsewhere
        failed_color = QColor("#FF5722")

        for row_idx, (_, row) in enumerate(df.iterrows()):
            row_status = str(row.iloc[status_idx]).lower() if status_idx is not None else ""
            row_color = None
            if row_status == "cancelled":
                row_color = cancelled_color
            elif row_status in fail_statuses:
                row_color = failed_color
            for col_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                if row_color is not None:
                    item.setForeground(row_color)
                self.results_table.setItem(row_idx, col_idx, item)

        self.results_table.resizeColumnsToContents()
        self.results_table.horizontalHeader().setStretchLastSection(True)

    def _update_summary_label(self, phase_name: str, df: pd.DataFrame) -> None:
        if "status" in df.columns:
            counts = df["status"].value_counts()
            counts_text = ", ".join(f"{status}: {count}" for status, count in counts.items())
        else:
            counts_text = "n/a"
        self.summary_label.setText(f"{phase_name} - {len(df)} row(s) - {counts_text}")


def main() -> None:
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(APP_STYLESHEET)
    if ICON_PATH.exists():
        app.setWindowIcon(QIcon(str(ICON_PATH)))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
