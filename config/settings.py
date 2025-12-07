"""
Application configuration and settings management.
"""

import json
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional


@dataclass
class AppSettings:
    """Application-wide settings."""
    
    # Display
    theme: str = "dark"
    accent_color: str = "#00bcd4"  # Cyan accent
    font_size: int = 10
    
    # Units
    default_unit_system: str = "imperial"  # "imperial" or "si"
    
    # Uncertainty
    monte_carlo_samples: int = 10000
    confidence_level: float = 0.95
    
    # Paths
    profiles_dir: str = ""
    last_import_dir: str = ""
    last_export_dir: str = ""
    
    def __post_init__(self):
        if not self.profiles_dir:
            self.profiles_dir = str(Path(__file__).parent.parent / "profiles")


class SettingsManager:
    """Manages loading and saving application settings."""
    
    def __init__(self, settings_file: Optional[Path] = None):
        if settings_file is None:
            settings_file = Path(__file__).parent / "user_settings.json"
        self.settings_file = settings_file
        self.settings = self.load()
    
    def load(self) -> AppSettings:
        """Load settings from file, or return defaults."""
        if self.settings_file.exists():
            try:
                with open(self.settings_file, 'r') as f:
                    data = json.load(f)
                return AppSettings(**data)
            except (json.JSONDecodeError, TypeError):
                pass
        return AppSettings()
    
    def save(self) -> None:
        """Save current settings to file."""
        self.settings_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.settings_file, 'w') as f:
            json.dump(asdict(self.settings), f, indent=2)
    
    def reset_to_defaults(self) -> None:
        """Reset all settings to defaults."""
        self.settings = AppSettings()
        self.save()


# Global settings instance
_settings_manager: Optional[SettingsManager] = None


def get_settings() -> AppSettings:
    """Get the global settings instance."""
    global _settings_manager
    if _settings_manager is None:
        _settings_manager = SettingsManager()
    return _settings_manager.settings


def save_settings() -> None:
    """Save the global settings."""
    global _settings_manager
    if _settings_manager is not None:
        _settings_manager.save()
