"""
Import dialog for CSV/Excel files.
"""

from pathlib import Path
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QLineEdit, QPushButton, QComboBox, QSpinBox, QFileDialog,
    QTableWidget, QTableWidgetItem, QDialogButtonBox, QMessageBox
)
from PySide6.QtCore import Signal

from core.data_io import DataImporter


class ImportDialog(QDialog):
    """Dialog for importing data from CSV/Excel files."""
    
    dataImported = Signal(object)  # Emits DataFrame
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.importer = DataImporter()
        self.current_data = None
        self._setup_ui()
    
    def _setup_ui(self):
        self.setWindowTitle("Import Data")
        self.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(self)
        
        # File selection
        file_layout = QHBoxLayout()
        self.file_path = QLineEdit()
        self.file_path.setReadOnly(True)
        file_layout.addWidget(self.file_path)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_file)
        file_layout.addWidget(browse_btn)
        
        layout.addLayout(file_layout)
        
        # Options
        options_layout = QFormLayout()
        
        self.sheet_combo = QComboBox()
        self.sheet_combo.setEnabled(False)
        options_layout.addRow("Sheet:", self.sheet_combo)
        
        self.header_spin = QSpinBox()
        self.header_spin.setRange(0, 100)
        self.header_spin.setValue(0)
        options_layout.addRow("Header Row:", self.header_spin)
        
        self.skip_spin = QSpinBox()
        self.skip_spin.setRange(0, 100)
        self.skip_spin.setValue(0)
        options_layout.addRow("Skip Rows:", self.skip_spin)
        
        layout.addLayout(options_layout)
        
        # Preview button
        preview_btn = QPushButton("Preview")
        preview_btn.clicked.connect(self._preview)
        layout.addWidget(preview_btn)
        
        # Preview table
        self.preview_table = QTableWidget()
        self.preview_table.setMaximumHeight(200)
        layout.addWidget(self.preview_table)
        
        # Status
        self.status_label = QLabel()
        layout.addWidget(self.status_label)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept_import)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
    
    def _browse_file(self):
        """Open file browser."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Select Data File",
            "",
            "Data Files (*.csv *.xlsx *.xls);;All Files (*)"
        )
        
        if filepath:
            self.file_path.setText(filepath)
            
            # Update sheet selection for Excel
            if filepath.endswith(('.xlsx', '.xls')):
                sheets = self.importer.get_excel_sheets(filepath)
                self.sheet_combo.clear()
                self.sheet_combo.addItems(sheets)
                self.sheet_combo.setEnabled(True)
            else:
                self.sheet_combo.setEnabled(False)
    
    def _preview(self):
        """Preview the file contents."""
        filepath = self.file_path.text()
        if not filepath:
            return
        
        sheet = self.sheet_combo.currentText() if self.sheet_combo.isEnabled() else None
        
        result = self.importer.import_file(
            filepath,
            sheet_name=sheet,
            header_row=self.header_spin.value(),
            skip_rows=self.skip_spin.value()
        )
        
        if result.success:
            self.current_data = result.data
            self._show_preview(result.data)
            
            status = f"Loaded {result.row_count} rows, {result.column_count} columns"
            if result.warnings:
                status += f"\nWarnings: {'; '.join(result.warnings)}"
            self.status_label.setText(status)
        else:
            self.status_label.setText(f"Error: {result.error_message}")
            QMessageBox.warning(self, "Import Error", result.error_message)
    
    def _show_preview(self, df):
        """Show DataFrame preview in table."""
        self.preview_table.clear()
        
        # Show first 10 rows
        preview = df.head(10)
        
        self.preview_table.setColumnCount(len(preview.columns))
        self.preview_table.setRowCount(len(preview))
        self.preview_table.setHorizontalHeaderLabels([str(c) for c in preview.columns])
        
        for i, row in enumerate(preview.itertuples(index=False)):
            for j, value in enumerate(row):
                self.preview_table.setItem(i, j, QTableWidgetItem(str(value)))
    
    def _accept_import(self):
        """Accept and emit the imported data."""
        if self.current_data is not None:
            self.dataImported.emit(self.current_data)
            self.accept()
        else:
            QMessageBox.warning(self, "No Data", "Please preview data before importing.")
