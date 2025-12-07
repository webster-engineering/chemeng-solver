"""
Settings dialog for application configuration.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QComboBox, QSpinBox, QDoubleSpinBox, QPushButton, QGroupBox,
    QDialogButtonBox
)
from PySide6.QtCore import Signal

from config import get_settings, save_settings


class SettingsDialog(QDialog):
    """Dialog for editing application settings."""
    
    settingsChanged = Signal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.settings = get_settings()
        self._setup_ui()
        self._load_settings()
    
    def _setup_ui(self):
        self.setWindowTitle("Settings")
        self.setMinimumWidth(400)
        
        layout = QVBoxLayout(self)
        
        # Display settings
        display_group = QGroupBox("Display")
        display_layout = QFormLayout(display_group)
        
        self.theme_combo = QComboBox()
        self.theme_combo.addItems(["dark", "light"])
        display_layout.addRow("Theme:", self.theme_combo)
        
        self.font_size_spin = QSpinBox()
        self.font_size_spin.setRange(8, 16)
        display_layout.addRow("Font Size:", self.font_size_spin)
        
        layout.addWidget(display_group)
        
        # Unit settings
        unit_group = QGroupBox("Units")
        unit_layout = QFormLayout(unit_group)
        
        self.unit_system_combo = QComboBox()
        self.unit_system_combo.addItems(["imperial", "si"])
        unit_layout.addRow("Default Unit System:", self.unit_system_combo)
        
        layout.addWidget(unit_group)
        
        # Uncertainty settings
        unc_group = QGroupBox("Uncertainty Analysis")
        unc_layout = QFormLayout(unc_group)
        
        self.mc_samples_spin = QSpinBox()
        self.mc_samples_spin.setRange(1000, 100000)
        self.mc_samples_spin.setSingleStep(1000)
        unc_layout.addRow("Monte Carlo Samples:", self.mc_samples_spin)
        
        self.confidence_spin = QDoubleSpinBox()
        self.confidence_spin.setRange(0.8, 0.99)
        self.confidence_spin.setSingleStep(0.01)
        self.confidence_spin.setDecimals(2)
        unc_layout.addRow("Confidence Level:", self.confidence_spin)
        
        layout.addWidget(unc_group)
        
        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.RestoreDefaults
        )
        buttons.accepted.connect(self._save_and_accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self._restore_defaults)
        
        layout.addWidget(buttons)
    
    def _load_settings(self):
        """Load current settings into widgets."""
        self.theme_combo.setCurrentText(self.settings.theme)
        self.font_size_spin.setValue(self.settings.font_size)
        self.unit_system_combo.setCurrentText(self.settings.default_unit_system)
        self.mc_samples_spin.setValue(self.settings.monte_carlo_samples)
        self.confidence_spin.setValue(self.settings.confidence_level)
    
    def _save_and_accept(self):
        """Save settings and close dialog."""
        self.settings.theme = self.theme_combo.currentText()
        self.settings.font_size = self.font_size_spin.value()
        self.settings.default_unit_system = self.unit_system_combo.currentText()
        self.settings.monte_carlo_samples = self.mc_samples_spin.value()
        self.settings.confidence_level = self.confidence_spin.value()
        
        save_settings()
        self.settingsChanged.emit()
        self.accept()
    
    def _restore_defaults(self):
        """Restore default settings."""
        from config.settings import AppSettings
        defaults = AppSettings()
        
        self.theme_combo.setCurrentText(defaults.theme)
        self.font_size_spin.setValue(defaults.font_size)
        self.unit_system_combo.setCurrentText(defaults.default_unit_system)
        self.mc_samples_spin.setValue(defaults.monte_carlo_samples)
        self.confidence_spin.setValue(defaults.confidence_level)
