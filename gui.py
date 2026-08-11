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
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QFontDatabase, QIcon, QPixmap
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

# Seconds a running phase gets to wind down after the window is closed before
# the process exits regardless. See MainWindow.closeEvent.
_FORCED_EXIT_GRACE_S: float = 10.0

PROVENANCE_PATH = Path(__file__).resolve().parent / "provenance.json"
# Kept in sync with provenance.json; used only if that file is missing or
# malformed so a broken/deleted provenance file can't stop the app launching.
_DEFAULT_PROVENANCE = {
    "author": "Maciej Gromelski",
    "version": "1.2.0",
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

        self._build_ui()
        self._build_footer()
        self._update_button_states()
        self._attach_logging()
        self._check_environment()

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

    def _build_footer(self) -> None:
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
        self.workdir_edit.textChanged.connect(self._on_workdir_changed)
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
            "",
            "Excel/CSV Files (*.xlsx *.csv);;Excel Files (*.xlsx);;"
            "CSV Files (*.csv);;All Files (*)",
        )
        if not path:
            return
        try:
            if path.lower().endswith(".csv"):
                df = pd.read_csv(path)
            else:
                df = pd.read_excel(path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Failed to load file", str(exc))
            return

        self.df = df
        self.input_path_label.setText(path)
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

    def _on_browse_workdir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select working directory")
        if path:
            self.workdir_edit.setText(path)

    def _on_workdir_changed(self, _text: str) -> None:
        self._prepare_succeeded = False
        self._mopac_succeeded = False
        self._last_prepare_df = None
        self._last_batch_df = None
        self._update_button_states()

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

        widgets = [
            self.keywords_edit,
            self.auto_charge_check,
            self.overwrite_check,
            self.prepare_workers_spin,
            self.mopac_workers_spin,
            self.export_workers_spin,
            self.obabel_timeout_spin,
            self.mopac_timeout_spin,
            self.export_timeout_spin,
            self.retry_obabel_timeout_spin,
            self.retry_mopac_timeout_spin,
            self.retry_keywords_edit,
            self.copy_sdf_check,
        ]
        for widget in widgets:
            widget.setEnabled(not lock)

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
        """Sole place that computes final enabled/disabled for the phase buttons."""
        can_prepare = (
            self.df is not None
            and bool(self.smiles_col_combo.currentText())
            and self._settings_locked
        )
        self.prepare_button.setEnabled(can_prepare and not self._phase_running)
        self.run_mopac_button.setEnabled(self._prepare_succeeded and not self._phase_running)
        self.export_sdf_button.setEnabled(self._mopac_succeeded and not self._phase_running)
        self.stop_button.setEnabled(self._phase_running)

        can_retry = self._mopac_succeeded and (
            self._has_failures(self._last_prepare_df) or self._has_failures(self._last_batch_df)
        )
        self.retry_failed_button.setEnabled(can_retry and not self._phase_running)

    # ------------------------------------------------------------------
    # Threaded phase execution
    # ------------------------------------------------------------------
    def _start_phase(self, func, kwargs: dict, phase_name: str) -> None:
        self._cancel_event.clear()
        kwargs = {**kwargs, "cancel_event": self._cancel_event}

        self._set_phase_buttons_enabled(True)
        self.progress_bar.setVisible(True)
        self._current_phase_name = phase_name

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

        self._poll_timer = QTimer(self)
        self._poll_timer.setInterval(150)
        self._poll_timer.timeout.connect(self._poll_worker)
        self._poll_timer.start()

    def _poll_worker(self) -> None:
        worker = self._worker
        if worker is None or not worker.done:
            return

        self._poll_timer.stop()
        self._worker = None
        self.progress_bar.setVisible(False)

        if worker.error is not None:
            QMessageBox.critical(self, f"{self._current_phase_name} failed", worker.error)
        else:
            # "Succeeded" means the phase call itself didn't raise -- not
            # that every row succeeded. Partial per-row failure is normal
            # pipeline behavior and must not block progressing to the next
            # phase for the rows that did succeed.
            if self._current_phase_name == "Prepare":
                self._prepare_succeeded = True
                self._last_prepare_df = worker.result
            elif self._current_phase_name == "Run MOPAC":
                self._mopac_succeeded = True
                self._last_batch_df = worker.result
            elif self._current_phase_name == "Retry Failed":
                self._last_batch_df = worker.result
            self._populate_table(worker.result)
            self._update_summary_label(self._current_phase_name, worker.result)
            if self._cancel_event.is_set():
                self.summary_label.setText(self.summary_label.text() + "  (stopped by user)")

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
