"""
Export dialog for PDF reports.
"""

from pathlib import Path
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QLabel, QLineEdit,
    QTextEdit, QPushButton, QFileDialog, QDialogButtonBox,
    QCheckBox, QMessageBox
)
from PySide6.QtCore import Signal

from core.data_io import PDFExporter


class ExportDialog(QDialog):
    """Dialog for exporting results to PDF."""
    
    def __init__(self, equation_name: str, inputs: dict, outputs: dict,
                 validation: dict = None, uncertainty: dict = None, parent=None):
        super().__init__(parent)
        self.equation_name = equation_name
        self.inputs = inputs
        self.outputs = outputs
        self.validation = validation
        self.uncertainty = uncertainty
        self.exporter = PDFExporter()
        self._setup_ui()
    
    def _setup_ui(self):
        self.setWindowTitle("Export to PDF")
        self.setMinimumWidth(400)
        
        layout = QVBoxLayout(self)
        
        form = QFormLayout()
        
        self.title_input = QLineEdit()
        self.title_input.setText(f"{self.equation_name} Calculation Report")
        form.addRow("Report Title:", self.title_input)
        
        self.notes_input = QTextEdit()
        self.notes_input.setMaximumHeight(100)
        self.notes_input.setPlaceholderText("Optional notes to include in report...")
        form.addRow("Notes:", self.notes_input)
        
        layout.addLayout(form)
        
        # Options
        self.include_validation = QCheckBox("Include Validation Results")
        self.include_validation.setChecked(True)
        self.include_validation.setEnabled(self.validation is not None)
        layout.addWidget(self.include_validation)
        
        self.include_uncertainty = QCheckBox("Include Uncertainty Analysis")
        self.include_uncertainty.setChecked(True)
        self.include_uncertainty.setEnabled(self.uncertainty is not None)
        layout.addWidget(self.include_uncertainty)
        
        # File selection
        file_layout = QFormLayout()
        self.file_path = QLineEdit()
        self.file_path.setPlaceholderText("Select output file...")
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_file)
        
        file_h = QFormLayout()
        file_h.addRow(self.file_path, browse_btn)
        file_layout.addRow("Save As:", self.file_path)
        
        layout.addLayout(file_layout)
        layout.addWidget(browse_btn)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._export)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
    
    def _browse_file(self):
        """Select output file."""
        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save PDF Report",
            f"{self.equation_name}_report.pdf",
            "PDF Files (*.pdf)"
        )
        if filepath:
            self.file_path.setText(filepath)
    
    def _export(self):
        """Export the report."""
        filepath = self.file_path.text()
        if not filepath:
            QMessageBox.warning(self, "No File", "Please select output file.")
            return
        
        validation = self.validation if self.include_validation.isChecked() else None
        uncertainty = self.uncertainty if self.include_uncertainty.isChecked() else None
        
        success = self.exporter.export(
            filepath=filepath,
            title=self.title_input.text(),
            equation_name=self.equation_name,
            inputs=self.inputs,
            outputs=self.outputs,
            validation_results=validation,
            uncertainty_results=uncertainty,
            notes=self.notes_input.toPlainText()
        )
        
        if success:
            QMessageBox.information(self, "Export Complete", f"Report saved to:\n{filepath}")
            self.accept()
        else:
            QMessageBox.warning(self, "Export Failed", "Failed to create PDF report.")
