"""
Main application window.
"""

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QTreeWidget, QTreeWidgetItem, QLabel, QPushButton, QGroupBox,
    QToolBar, QStatusBar, QMessageBox, QCheckBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QIcon

from .themes import get_theme_stylesheet
from .widgets.input_form import InputForm
from .widgets.results_view import ResultsView
from .widgets.chart_widget import ChartWidget
from .widgets.profile_manager import ProfileManagerWidget
from .dialogs.settings_dialog import SettingsDialog
from .dialogs.import_dialog import ImportDialog
from .dialogs.export_dialog import ExportDialog

from config import get_settings
from equations.process_control import PROCESS_CONTROL_EQUATIONS
from equations.fluid_dynamics import FLUID_DYNAMICS_EQUATIONS
from equations.heat_transfer import HEAT_TRANSFER_EQUATIONS
from equations.mass_transfer import MASS_TRANSFER_EQUATIONS
from equations.reaction_kinetics import REACTION_KINETICS_EQUATIONS
from equations.thermodynamics import THERMODYNAMICS_EQUATIONS
from equations.basic_math import BASIC_MATH_EQUATIONS
from equations.piping import PIPING_EQUATIONS
from equations.instrumentation import INSTRUMENTATION_EQUATIONS
from equations.distillation import DISTILLATION_EQUATIONS
from equations.pumps import PUMPS_EQUATIONS
from equations.vessels import VESSELS_EQUATIONS
from equations.safety import SAFETY_EQUATIONS
from equations.economics import ECONOMICS_EQUATIONS
from core.uncertainty import UncertainParameter, DistributionType


# All equation categories
EQUATION_REGISTRY = {
    "Process Control": PROCESS_CONTROL_EQUATIONS,
    "Fluid Dynamics": FLUID_DYNAMICS_EQUATIONS,
    "Piping": PIPING_EQUATIONS,
    "Heat Transfer": HEAT_TRANSFER_EQUATIONS,
    "Thermodynamics": THERMODYNAMICS_EQUATIONS,
    "Mass Transfer": MASS_TRANSFER_EQUATIONS,
    "Distillation": DISTILLATION_EQUATIONS,
    "Reaction Kinetics": REACTION_KINETICS_EQUATIONS,
    "Pumps & Compressors": PUMPS_EQUATIONS,
    "Vessels & Tanks": VESSELS_EQUATIONS,
    "Instrumentation": INSTRUMENTATION_EQUATIONS,
    "Safety & Relief": SAFETY_EQUATIONS,
    "Economics": ECONOMICS_EQUATIONS,
    "Basic Math": BASIC_MATH_EQUATIONS,
}


class MainWindow(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        super().__init__()
        self.settings = get_settings()
        self.current_equation = None
        self.current_result = None
        self._setup_ui()
        self._apply_theme()
    
    def _setup_ui(self):
        self.setWindowTitle("Chemical Engineering Equation Solver")
        self.setMinimumSize(1200, 800)
        
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)
        
        # Main splitter
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - Equation navigation
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        
        nav_label = QLabel("Equations")
        nav_label.setObjectName("header")
        left_layout.addWidget(nav_label)
        
        self.equation_tree = QTreeWidget()
        self.equation_tree.setHeaderHidden(True)
        self.equation_tree.itemClicked.connect(self._on_equation_selected)
        self._populate_equation_tree()
        left_layout.addWidget(self.equation_tree)
        
        # Profile section
        profile_label = QLabel("Saved Profiles")
        left_layout.addWidget(profile_label)
        
        self.profile_manager = ProfileManagerWidget()
        self.profile_manager.profileLoaded.connect(self._on_profile_loaded)
        self.profile_manager.setMaximumHeight(200)
        left_layout.addWidget(self.profile_manager)
        
        splitter.addWidget(left_panel)
        
        # Center panel - Input form
        center_panel = QWidget()
        center_layout = QVBoxLayout(center_panel)
        
        self.equation_header = QLabel("Select an equation")
        self.equation_header.setObjectName("header")
        center_layout.addWidget(self.equation_header)
        
        self.equation_desc = QLabel()
        self.equation_desc.setWordWrap(True)
        center_layout.addWidget(self.equation_desc)
        
        input_group = QGroupBox("Inputs")
        input_layout = QVBoxLayout(input_group)
        self.input_form = InputForm()
        input_layout.addWidget(self.input_form)
        center_layout.addWidget(input_group)
        
        # Options
        options_layout = QHBoxLayout()
        self.uncertainty_check = QCheckBox("Run Uncertainty Analysis")
        options_layout.addWidget(self.uncertainty_check)
        options_layout.addStretch()
        center_layout.addLayout(options_layout)
        
        # Calculate button
        self.calculate_btn = QPushButton("Calculate")
        self.calculate_btn.setObjectName("calculate_btn")
        self.calculate_btn.clicked.connect(self._on_calculate)
        self.calculate_btn.setEnabled(False)
        center_layout.addWidget(self.calculate_btn)
        
        center_layout.addStretch()
        splitter.addWidget(center_panel)
        
        # Right panel - Results
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        self.results_view = ResultsView()
        right_layout.addWidget(self.results_view)
        
        self.chart_widget = ChartWidget()
        self.chart_widget.setMinimumHeight(300)
        right_layout.addWidget(self.chart_widget)
        
        splitter.addWidget(right_panel)
        
        # Set splitter proportions
        splitter.setSizes([250, 400, 450])
        layout.addWidget(splitter)
        
        # Setup toolbar and menu
        self._setup_toolbar()
        self._setup_menu()
        self._setup_statusbar()
    
    def _populate_equation_tree(self):
        """Populate the equation navigation tree."""
        for category, equations in EQUATION_REGISTRY.items():
            category_item = QTreeWidgetItem([category])
            category_item.setFlags(category_item.flags() & ~Qt.ItemIsSelectable)
            
            for eq_id, eq_class in equations.items():
                eq = eq_class()
                eq_item = QTreeWidgetItem([eq.name])
                eq_item.setData(0, Qt.UserRole, (category, eq_id))
                eq_item.setToolTip(0, eq.description)
                category_item.addChild(eq_item)
            
            self.equation_tree.addTopLevelItem(category_item)
        
        self.equation_tree.expandAll()
    
    def _setup_toolbar(self):
        toolbar = QToolBar("Main Toolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        
        # Save profile
        save_action = QAction("Save Profile", self)
        save_action.triggered.connect(self._on_save_profile)
        toolbar.addAction(save_action)
        
        toolbar.addSeparator()
        
        # Import
        import_action = QAction("Import Data", self)
        import_action.triggered.connect(self._on_import)
        toolbar.addAction(import_action)
        
        # Export
        export_action = QAction("Export PDF", self)
        export_action.triggered.connect(self._on_export)
        toolbar.addAction(export_action)
        
        toolbar.addSeparator()
        
        # Settings
        settings_action = QAction("Settings", self)
        settings_action.triggered.connect(self._on_settings)
        toolbar.addAction(settings_action)
    
    def _setup_menu(self):
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        file_menu.addAction("Import Data...", self._on_import)
        file_menu.addAction("Export PDF...", self._on_export)
        file_menu.addSeparator()
        file_menu.addAction("Exit", self.close)
        
        # Profile menu
        profile_menu = menubar.addMenu("Profile")
        profile_menu.addAction("Save Profile...", self._on_save_profile)
        
        # Settings menu
        settings_menu = menubar.addMenu("Settings")
        settings_menu.addAction("Preferences...", self._on_settings)
    
    def _setup_statusbar(self):
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        
        self.unit_label = QLabel(f"Units: {self.settings.default_unit_system.upper()}")
        self.statusbar.addPermanentWidget(self.unit_label)
    
    def _apply_theme(self):
        """Apply the current theme."""
        stylesheet = get_theme_stylesheet(self.settings.theme)
        self.setStyleSheet(stylesheet)
    
    def _on_equation_selected(self, item: QTreeWidgetItem, column: int):
        """Handle equation selection."""
        data = item.data(0, Qt.UserRole)
        if data is None:
            return
        
        category, eq_id = data
        eq_class = EQUATION_REGISTRY[category][eq_id]
        self.current_equation = eq_class()
        
        # Update header
        self.equation_header.setText(self.current_equation.name)
        self.equation_desc.setText(self.current_equation.description)
        
        # Load parameters into form
        self.input_form.load_parameters(self.current_equation.parameters)
        
        # Clear results
        self.results_view.clear()
        self.chart_widget.clear()
        
        self.calculate_btn.setEnabled(True)
        self.statusbar.showMessage(f"Selected: {self.current_equation.name}")
    
    def _on_calculate(self):
        """Run the calculation."""
        if self.current_equation is None:
            return
        
        # Validate inputs
        valid, error = self.input_form.validate()
        if not valid:
            QMessageBox.warning(self, "Invalid Input", error)
            return
        
        # Get inputs
        inputs = self.input_form.get_values()
        run_uncertainty = self.uncertainty_check.isChecked()
        
        # Build uncertainty parameters if requested
        uncertain_params = None
        if run_uncertainty:
            uncertainties = self.input_form.get_uncertainties()
            uncertain_params = []
            for name, (val, unit) in inputs.items():
                rel_std = uncertainties.get(name, 0.05)
                uncertain_params.append(UncertainParameter(
                    name=name, nominal_value=val, unit=unit or "",
                    distribution=DistributionType.NORMAL, relative_std=rel_std
                ))
        
        # Run calculation
        self.statusbar.showMessage("Calculating...")
        result = self.current_equation.calculate(
            inputs, validate=True, with_uncertainty=run_uncertainty,
            uncertain_params=uncertain_params
        )
        
        self.current_result = result
        
        # Display results
        self.results_view.display_result(result)
        
        # Plot chart for process control equations
        if result.success and "response" in self.current_equation.equation_id:
            self._plot_response(inputs)
        elif result.success and result.uncertainty:
            # Plot first output histogram
            first_output = list(result.uncertainty.keys())[0]
            unc = result.uncertainty[first_output]
            if unc.samples is not None:
                self.chart_widget.plot_histogram(unc.samples, first_output, "Uncertainty Distribution")
        
        status = "Calculation complete" if result.success else f"Error: {result.error_message}"
        self.statusbar.showMessage(status)
    
    def _plot_response(self, inputs):
        """Plot step response for control equations."""
        try:
            K = inputs.get('K', (1, ''))[0] if isinstance(inputs.get('K'), tuple) else inputs.get('K', 1)
            tau = inputs.get('tau', (1, 'min'))[0] if isinstance(inputs.get('tau'), tuple) else inputs.get('tau', 1)
            zeta = inputs.get('zeta', (1, ''))[0] if isinstance(inputs.get('zeta'), tuple) else 1.0
            theta = inputs.get('theta', (0, 'min'))[0] if isinstance(inputs.get('theta'), tuple) else 0
            
            self.chart_widget.plot_step_response(K, tau, zeta, theta)
        except Exception:
            pass
    
    def _on_profile_loaded(self, profile):
        """Handle profile load."""
        # Find and select the equation
        category = profile.category
        eq_id = profile.equation_id
        
        if category in EQUATION_REGISTRY and eq_id in EQUATION_REGISTRY[category]:
            eq_class = EQUATION_REGISTRY[category][eq_id]
            self.current_equation = eq_class()
            
            self.equation_header.setText(self.current_equation.name)
            self.equation_desc.setText(self.current_equation.description)
            self.input_form.load_parameters(self.current_equation.parameters)
            self.input_form.set_values(profile.inputs)
            
            self.calculate_btn.setEnabled(True)
            self.statusbar.showMessage(f"Loaded profile: {profile.name}")
    
    def _on_save_profile(self):
        """Save current calculation as profile."""
        if self.current_equation is None:
            QMessageBox.warning(self, "No Equation", "Please select an equation first.")
            return
        
        inputs = {}
        for name, (val, unit) in self.input_form.get_values().items():
            inputs[name] = {'value': val, 'unit': unit}
        
        category = self.current_equation.category
        self.profile_manager.save_profile(
            self.current_equation.equation_id,
            self.current_equation.name,
            category, inputs
        )
    
    def _on_import(self):
        """Show import dialog."""
        dialog = ImportDialog(self)
        dialog.dataImported.connect(self._handle_imported_data)
        dialog.exec()
    
    def _handle_imported_data(self, df):
        """Handle imported data."""
        self.statusbar.showMessage(f"Imported {len(df)} rows")
    
    def _on_export(self):
        """Show export dialog."""
        if self.current_equation is None or self.current_result is None:
            QMessageBox.warning(self, "No Results", "Please run a calculation first.")
            return
        
        inputs = {}
        for name, (val, unit) in self.input_form.get_values().items():
            inputs[name] = {'value': val, 'unit': unit}
        
        outputs = {}
        for name, qty in self.current_result.outputs.items():
            outputs[name] = {'value': qty.magnitude, 'unit': str(qty.units)}
        
        validation = None
        if self.current_result.validation:
            validation = {
                'status': self.current_result.validation.overall_status.value,
                'messages': [msg.message for msg in self.current_result.validation.messages]
            }
        
        uncertainty = None
        if self.current_result.uncertainty:
            uncertainty = {name: unc.to_dict() for name, unc in self.current_result.uncertainty.items()}
        
        dialog = ExportDialog(
            self.current_equation.name, inputs, outputs,
            validation, uncertainty, self
        )
        dialog.exec()
    
    def _on_settings(self):
        """Show settings dialog."""
        dialog = SettingsDialog(self)
        dialog.settingsChanged.connect(self._apply_theme)
        dialog.settingsChanged.connect(self._update_status)
        dialog.exec()
    
    def _update_status(self):
        """Update status bar after settings change."""
        self.settings = get_settings()
        self.unit_label.setText(f"Units: {self.settings.default_unit_system.upper()}")
