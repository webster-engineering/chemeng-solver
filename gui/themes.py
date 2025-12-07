"""
Theme definitions for the application.
"""

DARK_THEME = """
QMainWindow, QWidget {
    background-color: #1e1e2e;
    color: #cdd6f4;
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 10pt;
}

QMenuBar {
    background-color: #181825;
    color: #cdd6f4;
    border-bottom: 1px solid #313244;
}

QMenuBar::item:selected {
    background-color: #45475a;
}

QMenu {
    background-color: #1e1e2e;
    color: #cdd6f4;
    border: 1px solid #313244;
}

QMenu::item:selected {
    background-color: #45475a;
}

QToolBar {
    background-color: #181825;
    border: none;
    spacing: 5px;
    padding: 5px;
}

QPushButton {
    background-color: #45475a;
    color: #cdd6f4;
    border: 1px solid #585b70;
    border-radius: 6px;
    padding: 8px 16px;
    min-width: 80px;
}

QPushButton:hover {
    background-color: #585b70;
    border-color: #89b4fa;
}

QPushButton:pressed {
    background-color: #313244;
}

QPushButton:disabled {
    background-color: #313244;
    color: #6c7086;
}

QPushButton#calculate_btn {
    background-color: #89b4fa;
    color: #1e1e2e;
    font-weight: bold;
}

QPushButton#calculate_btn:hover {
    background-color: #b4befe;
}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background-color: #313244;
    color: #cdd6f4;
    border: 1px solid #45475a;
    border-radius: 4px;
    padding: 6px 10px;
}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border-color: #89b4fa;
}

QComboBox::drop-down {
    border: none;
    width: 30px;
}

QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid #cdd6f4;
    margin-right: 10px;
}

QComboBox QAbstractItemView {
    background-color: #313244;
    color: #cdd6f4;
    selection-background-color: #45475a;
}

QListWidget, QTreeWidget, QTableWidget {
    background-color: #1e1e2e;
    color: #cdd6f4;
    border: 1px solid #313244;
    alternate-background-color: #181825;
}

QListWidget::item:selected, QTreeWidget::item:selected {
    background-color: #45475a;
}

QTreeWidget::item:hover {
    background-color: #313244;
}

QHeaderView::section {
    background-color: #181825;
    color: #cdd6f4;
    padding: 8px;
    border: none;
    border-right: 1px solid #313244;
    border-bottom: 1px solid #313244;
}

QTabWidget::pane {
    border: 1px solid #313244;
    background-color: #1e1e2e;
}

QTabBar::tab {
    background-color: #181825;
    color: #a6adc8;
    padding: 10px 20px;
    border: 1px solid #313244;
    border-bottom: none;
    margin-right: 2px;
}

QTabBar::tab:selected {
    background-color: #1e1e2e;
    color: #89b4fa;
    border-bottom: 2px solid #89b4fa;
}

QTabBar::tab:hover:!selected {
    background-color: #313244;
}

QGroupBox {
    font-weight: bold;
    border: 1px solid #313244;
    border-radius: 6px;
    margin-top: 12px;
    padding-top: 10px;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 5px;
    color: #89b4fa;
}

QLabel {
    color: #cdd6f4;
}

QLabel#header {
    font-size: 14pt;
    font-weight: bold;
    color: #89b4fa;
}

QLabel#validation_pass {
    color: #a6e3a1;
    font-weight: bold;
}

QLabel#validation_warning {
    color: #f9e2af;
    font-weight: bold;
}

QLabel#validation_fail {
    color: #f38ba8;
    font-weight: bold;
}

QScrollBar:vertical {
    background-color: #1e1e2e;
    width: 12px;
    border-radius: 6px;
}

QScrollBar::handle:vertical {
    background-color: #45475a;
    border-radius: 6px;
    min-height: 30px;
}

QScrollBar::handle:vertical:hover {
    background-color: #585b70;
}

QScrollBar:horizontal {
    background-color: #1e1e2e;
    height: 12px;
    border-radius: 6px;
}

QScrollBar::handle:horizontal {
    background-color: #45475a;
    border-radius: 6px;
    min-width: 30px;
}

QScrollArea {
    border: none;
}

QStatusBar {
    background-color: #181825;
    color: #a6adc8;
    border-top: 1px solid #313244;
}

QSplitter::handle {
    background-color: #313244;
}

QSplitter::handle:horizontal {
    width: 2px;
}

QSplitter::handle:vertical {
    height: 2px;
}

QCheckBox {
    spacing: 8px;
}

QCheckBox::indicator {
    width: 18px;
    height: 18px;
    border: 2px solid #45475a;
    border-radius: 4px;
    background-color: #313244;
}

QCheckBox::indicator:checked {
    background-color: #89b4fa;
    border-color: #89b4fa;
}

QToolTip {
    background-color: #313244;
    color: #cdd6f4;
    border: 1px solid #45475a;
    padding: 5px;
}

QProgressBar {
    background-color: #313244;
    border: none;
    border-radius: 4px;
    text-align: center;
    color: #cdd6f4;
}

QProgressBar::chunk {
    background-color: #89b4fa;
    border-radius: 4px;
}
"""

LIGHT_THEME = """
QMainWindow, QWidget {
    background-color: #eff1f5;
    color: #4c4f69;
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 10pt;
}

QPushButton {
    background-color: #dce0e8;
    color: #4c4f69;
    border: 1px solid #bcc0cc;
    border-radius: 6px;
    padding: 8px 16px;
}

QPushButton:hover {
    background-color: #ccd0da;
    border-color: #1e66f5;
}

QPushButton#calculate_btn {
    background-color: #1e66f5;
    color: white;
    font-weight: bold;
}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background-color: #ffffff;
    color: #4c4f69;
    border: 1px solid #bcc0cc;
    border-radius: 4px;
    padding: 6px 10px;
}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border-color: #1e66f5;
}

QListWidget, QTreeWidget, QTableWidget {
    background-color: #ffffff;
    color: #4c4f69;
    border: 1px solid #ccd0da;
    alternate-background-color: #e6e9ef;
}

QGroupBox {
    border: 1px solid #ccd0da;
    border-radius: 6px;
}

QGroupBox::title {
    color: #1e66f5;
}

QLabel#validation_pass {
    color: #40a02b;
    font-weight: bold;
}

QLabel#validation_warning {
    color: #df8e1d;
    font-weight: bold;
}

QLabel#validation_fail {
    color: #d20f39;
    font-weight: bold;
}
"""

THEMES = {
    'dark': DARK_THEME,
    'light': LIGHT_THEME
}


def get_theme_stylesheet(theme_name: str = 'dark') -> str:
    """Get the stylesheet for a theme."""
    return THEMES.get(theme_name, DARK_THEME)
