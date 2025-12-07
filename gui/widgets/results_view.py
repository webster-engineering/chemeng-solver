"""
Results display widget with validation status.
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QGroupBox, QScrollArea, QFrame, QTextEdit
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from equations.base import EquationResult
from core.validation import ValidationStatus
from core.units import get_converter


class ValidationBadge(QLabel):
    """Status badge for validation results."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumWidth(80)
        self.set_status(ValidationStatus.SKIPPED)
    
    def set_status(self, status: ValidationStatus):
        self.status = status
        
        if status == ValidationStatus.PASS:
            self.setText("✓ PASS")
            self.setObjectName("validation_pass")
        elif status == ValidationStatus.WARNING:
            self.setText("⚠ WARNING")
            self.setObjectName("validation_warning")
        elif status == ValidationStatus.FAIL:
            self.setText("✗ FAIL")
            self.setObjectName("validation_fail")
        else:
            self.setText("— SKIPPED")
            self.setObjectName("")
        
        # Force style refresh
        self.style().unpolish(self)
        self.style().polish(self)


class OutputValueWidget(QFrame):
    """Widget to display a single output value."""
    
    def __init__(self, name: str, value: str, unit: str, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 8)
        
        # Parameter name
        name_label = QLabel(name)
        name_label.setStyleSheet("color: #89b4fa; font-weight: bold;")
        layout.addWidget(name_label)
        
        # Value with unit
        value_layout = QHBoxLayout()
        
        value_label = QLabel(value)
        font = value_label.font()
        font.setPointSize(14)
        font.setBold(True)
        value_label.setFont(font)
        value_layout.addWidget(value_label)
        
        if unit:
            unit_label = QLabel(unit)
            unit_label.setStyleSheet("color: #a6adc8;")
            value_layout.addWidget(unit_label)
        
        value_layout.addStretch()
        layout.addLayout(value_layout)


class ResultsView(QWidget):
    """Widget to display calculation results with validation."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Validation status section
        validation_group = QGroupBox("Validation Status")
        validation_layout = QHBoxLayout(validation_group)
        
        self.validation_badge = ValidationBadge()
        validation_layout.addWidget(self.validation_badge)
        
        self.validation_message = QLabel()
        self.validation_message.setWordWrap(True)
        validation_layout.addWidget(self.validation_message, stretch=1)
        
        layout.addWidget(validation_group)
        
        # Results section
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)
        
        self.results_scroll = QScrollArea()
        self.results_scroll.setWidgetResizable(True)
        self.results_scroll.setFrameShape(QScrollArea.NoFrame)
        
        self.results_container = QWidget()
        self.results_flow = QVBoxLayout(self.results_container)
        self.results_flow.setSpacing(8)
        
        self.results_scroll.setWidget(self.results_container)
        results_layout.addWidget(self.results_scroll)
        
        layout.addWidget(results_group, stretch=1)
        
        # Uncertainty section (collapsible)
        self.uncertainty_group = QGroupBox("Uncertainty Analysis")
        self.uncertainty_group.setVisible(False)
        uncertainty_layout = QVBoxLayout(self.uncertainty_group)
        
        self.uncertainty_text = QTextEdit()
        self.uncertainty_text.setReadOnly(True)
        self.uncertainty_text.setMaximumHeight(150)
        uncertainty_layout.addWidget(self.uncertainty_text)
        
        layout.addWidget(self.uncertainty_group)
        
        # Details section
        self.details_group = QGroupBox("Validation Details")
        self.details_group.setVisible(False)
        details_layout = QVBoxLayout(self.details_group)
        
        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        self.details_text.setMaximumHeight(100)
        details_layout.addWidget(self.details_text)
        
        layout.addWidget(self.details_group)
    
    def display_result(self, result: EquationResult):
        """Display calculation results."""
        self._clear_results()
        
        if not result.success:
            self.validation_badge.set_status(ValidationStatus.FAIL)
            self.validation_message.setText(f"Calculation failed: {result.error_message}")
            return
        
        converter = get_converter()
        
        # Display outputs
        for name, quantity in result.outputs.items():
            value_str = f"{quantity.magnitude:.6g}"
            unit_str = str(quantity.units) if not quantity.dimensionless else ""
            
            output_widget = OutputValueWidget(name, value_str, unit_str)
            self.results_flow.addWidget(output_widget)
        
        self.results_flow.addStretch()
        
        # Display validation
        if result.validation:
            self.validation_badge.set_status(result.validation.overall_status)
            
            messages = [msg.message for msg in result.validation.messages]
            self.validation_message.setText("; ".join(messages) if messages else "All checks passed")
            
            # Show details if any
            details = []
            for msg in result.validation.messages:
                if msg.details:
                    details.append(f"{msg.message}: {msg.details}")
            
            if details:
                self.details_text.setText("\n".join(details))
                self.details_group.setVisible(True)
        else:
            self.validation_badge.set_status(ValidationStatus.SKIPPED)
            self.validation_message.setText("Validation not performed")
        
        # Display uncertainty
        if result.uncertainty:
            uncertainty_lines = []
            for name, unc_result in result.uncertainty.items():
                uncertainty_lines.append(
                    f"{name}: {unc_result.mean:.4g} ± {unc_result.std_dev:.4g} "
                    f"(95% CI: [{unc_result.ci_lower:.4g}, {unc_result.ci_upper:.4g}])"
                )
            self.uncertainty_text.setText("\n".join(uncertainty_lines))
            self.uncertainty_group.setVisible(True)
    
    def _clear_results(self):
        """Clear all displayed results."""
        # Clear output widgets
        while self.results_flow.count():
            item = self.results_flow.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.validation_badge.set_status(ValidationStatus.SKIPPED)
        self.validation_message.clear()
        self.uncertainty_text.clear()
        self.uncertainty_group.setVisible(False)
        self.details_text.clear()
        self.details_group.setVisible(False)
    
    def clear(self):
        """Clear the view."""
        self._clear_results()
