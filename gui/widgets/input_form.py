"""
Dynamic input form widget for equations.
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QLineEdit, QComboBox, QCheckBox, QGroupBox, QScrollArea,
    QDoubleSpinBox
)
from PySide6.QtCore import Signal, Qt

from equations.base import EquationParameter, ParameterType
from core.units import UNIT_CATEGORIES, get_converter


class ParameterInput(QWidget):
    """Widget for a single parameter input."""
    
    valueChanged = Signal()
    
    def __init__(self, param: EquationParameter, parent=None):
        super().__init__(parent)
        self.param = param
        self.converter = get_converter()
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Value input
        self.value_input = QLineEdit()
        self.value_input.setPlaceholderText(f"Enter {self.param.symbol or self.param.name}")
        self.value_input.textChanged.connect(self.valueChanged)
        layout.addWidget(self.value_input, stretch=2)
        
        # Unit selector
        self.unit_combo = QComboBox()
        self.unit_combo.setMinimumWidth(80)
        
        if self.param.unit_category and self.param.unit_category in UNIT_CATEGORIES:
            units = UNIT_CATEGORIES[self.param.unit_category].common_units
            self.unit_combo.addItems([u for u in units if u])
            
            # Set default
            default = self.param.default_unit
            if default in units:
                self.unit_combo.setCurrentText(default)
        elif self.param.default_unit:
            self.unit_combo.addItem(self.param.default_unit)
        else:
            self.unit_combo.addItem("—")
            self.unit_combo.setEnabled(False)
        
        self.unit_combo.currentTextChanged.connect(self.valueChanged)
        layout.addWidget(self.unit_combo)
        
        # Uncertainty checkbox and input
        self.uncertainty_check = QCheckBox("±")
        self.uncertainty_check.setToolTip("Include uncertainty")
        self.uncertainty_check.toggled.connect(self._toggle_uncertainty)
        layout.addWidget(self.uncertainty_check)
        
        self.uncertainty_input = QLineEdit()
        self.uncertainty_input.setMaximumWidth(60)
        self.uncertainty_input.setPlaceholderText("%")
        self.uncertainty_input.setVisible(False)
        layout.addWidget(self.uncertainty_input)
        
        self.setToolTip(self.param.tooltip or self.param.description)
    
    def _toggle_uncertainty(self, checked):
        self.uncertainty_input.setVisible(checked)
    
    def get_value(self):
        """Get current value as (value, unit) tuple."""
        try:
            value = float(self.value_input.text())
            unit = self.unit_combo.currentText()
            if unit == "—":
                unit = ""
            return value, unit
        except ValueError:
            return None, None
    
    def get_uncertainty(self):
        """Get uncertainty as relative fraction."""
        if self.uncertainty_check.isChecked():
            try:
                return float(self.uncertainty_input.text()) / 100
            except ValueError:
                return 0.05  # Default 5%
        return None
    
    def set_value(self, value, unit=None):
        """Set the input value."""
        self.value_input.setText(str(value))
        if unit and self.unit_combo.findText(unit) >= 0:
            self.unit_combo.setCurrentText(unit)
    
    def clear(self):
        """Clear the input."""
        self.value_input.clear()
        self.uncertainty_check.setChecked(False)
        self.uncertainty_input.clear()


class InputForm(QWidget):
    """Dynamic input form generated from equation parameters."""
    
    valuesChanged = Signal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.param_inputs = {}
        self._setup_ui()
    
    def _setup_ui(self):
        self.main_layout = QVBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        
        # Scroll area for many parameters
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        
        self.form_container = QWidget()
        self.form_layout = QFormLayout(self.form_container)
        self.form_layout.setLabelAlignment(Qt.AlignRight)
        self.form_layout.setSpacing(10)
        
        scroll.setWidget(self.form_container)
        self.main_layout.addWidget(scroll)
    
    def load_parameters(self, parameters: list):
        """Load parameters into the form."""
        # Clear existing
        self.clear()
        
        for param in parameters:
            if param.param_type != ParameterType.INPUT:
                continue
            
            # Create label
            label_text = f"{param.symbol or param.name}"
            if param.required:
                label_text += " *"
            label = QLabel(label_text)
            label.setToolTip(param.description)
            
            # Create input widget
            input_widget = ParameterInput(param)
            input_widget.valueChanged.connect(self.valuesChanged)
            
            self.form_layout.addRow(label, input_widget)
            self.param_inputs[param.name] = input_widget
    
    def get_values(self):
        """Get all input values as dict."""
        values = {}
        for name, widget in self.param_inputs.items():
            val, unit = widget.get_value()
            if val is not None:
                values[name] = (val, unit) if unit else val
        return values
    
    def get_uncertainties(self):
        """Get uncertainty values for parameters that have them."""
        uncertainties = {}
        for name, widget in self.param_inputs.items():
            unc = widget.get_uncertainty()
            if unc is not None:
                uncertainties[name] = unc
        return uncertainties
    
    def set_values(self, values: dict):
        """Set input values from dict."""
        for name, value in values.items():
            if name in self.param_inputs:
                if isinstance(value, dict):
                    self.param_inputs[name].set_value(value.get('value'), value.get('unit'))
                elif isinstance(value, tuple):
                    self.param_inputs[name].set_value(value[0], value[1] if len(value) > 1 else None)
                else:
                    self.param_inputs[name].set_value(value)
    
    def clear(self):
        """Clear all parameters."""
        for widget in self.param_inputs.values():
            widget.deleteLater()
        self.param_inputs.clear()
        
        # Clear form layout
        while self.form_layout.count():
            item = self.form_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
    
    def validate(self):
        """Check if all required inputs are filled."""
        for name, widget in self.param_inputs.items():
            if widget.param.required:
                val, _ = widget.get_value()
                if val is None:
                    return False, f"Missing required input: {name}"
        return True, ""
