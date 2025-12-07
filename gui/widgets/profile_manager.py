"""
Profile manager widget for saving and loading calculation profiles.
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QPushButton,
    QLineEdit, QLabel, QMessageBox, QInputDialog, QListWidgetItem
)
from PySide6.QtCore import Signal

from core.data_io import ProfileManager, CalculationProfile


class ProfileManagerWidget(QWidget):
    """Widget for managing saved profiles."""
    
    profileLoaded = Signal(CalculationProfile)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.manager = ProfileManager()
        self._setup_ui()
        self._refresh_list()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Profile list
        self.profile_list = QListWidget()
        self.profile_list.itemDoubleClicked.connect(self._on_load)
        layout.addWidget(self.profile_list)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.load_btn = QPushButton("Load")
        self.load_btn.clicked.connect(self._on_load)
        btn_layout.addWidget(self.load_btn)
        
        self.delete_btn = QPushButton("Delete")
        self.delete_btn.clicked.connect(self._on_delete)
        btn_layout.addWidget(self.delete_btn)
        
        layout.addLayout(btn_layout)
    
    def _refresh_list(self):
        """Refresh the list of profiles."""
        self.profile_list.clear()
        for name in self.manager.list_profiles():
            self.profile_list.addItem(name)
    
    def _on_load(self, item=None):
        """Load selected profile."""
        if item is None:
            item = self.profile_list.currentItem()
        
        if item:
            profile = self.manager.load_profile(item.text())
            if profile:
                self.profileLoaded.emit(profile)
    
    def _on_delete(self):
        """Delete selected profile."""
        item = self.profile_list.currentItem()
        if item:
            reply = QMessageBox.question(
                self, "Delete Profile",
                f"Delete profile '{item.text()}'?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.manager.delete_profile(item.text())
                self._refresh_list()
    
    def save_profile(self, equation_id: str, equation_name: str, 
                     category: str, inputs: dict) -> bool:
        """Save current calculation as a profile."""
        name, ok = QInputDialog.getText(self, "Save Profile", "Profile name:")
        
        if ok and name:
            profile = CalculationProfile(
                name=name,
                equation_id=equation_id,
                equation_name=equation_name,
                category=category,
                inputs=inputs
            )
            success = self.manager.save_profile(profile)
            if success:
                self._refresh_list()
            return success
        return False
    
    def refresh(self):
        """Refresh the profile list."""
        self._refresh_list()
