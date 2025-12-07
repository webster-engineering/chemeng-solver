"""
Data import/export and profile management functionality.
"""

import json
import csv
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Any, Optional, Union
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.platypus import PageBreak, KeepTogether

from .units import ureg
from config import get_settings


@dataclass
class DataImportResult:
    """Result of a data import operation."""
    success: bool
    data: Optional[pd.DataFrame] = None
    error_message: str = ""
    warnings: List[str] = field(default_factory=list)
    row_count: int = 0
    column_count: int = 0


class DataImporter:
    """
    Handles importing data from CSV and Excel files.
    """
    
    def __init__(self):
        self.supported_extensions = {'.csv', '.xlsx', '.xls'}
    
    def import_file(
        self,
        filepath: Union[str, Path],
        sheet_name: Optional[str] = None,
        header_row: int = 0,
        skip_rows: int = 0
    ) -> DataImportResult:
        """
        Import data from a file.
        
        Args:
            filepath: Path to CSV or Excel file
            sheet_name: Sheet name for Excel files (optional)
            header_row: Row index containing headers (0-based)
            skip_rows: Number of rows to skip at beginning
        
        Returns:
            DataImportResult with DataFrame or error
        """
        filepath = Path(filepath)
        
        if not filepath.exists():
            return DataImportResult(
                success=False,
                error_message=f"File not found: {filepath}"
            )
        
        ext = filepath.suffix.lower()
        
        if ext not in self.supported_extensions:
            return DataImportResult(
                success=False,
                error_message=f"Unsupported file type: {ext}. Supported: {self.supported_extensions}"
            )
        
        try:
            if ext == '.csv':
                df = self._import_csv(filepath, header_row, skip_rows)
            else:
                df = self._import_excel(filepath, sheet_name, header_row, skip_rows)
            
            warnings = self._validate_data(df)
            
            return DataImportResult(
                success=True,
                data=df,
                warnings=warnings,
                row_count=len(df),
                column_count=len(df.columns)
            )
        
        except Exception as e:
            return DataImportResult(
                success=False,
                error_message=str(e)
            )
    
    def _import_csv(self, filepath: Path, header_row: int, skip_rows: int) -> pd.DataFrame:
        """Import CSV file with auto-detection of delimiter."""
        # Try to detect delimiter
        with open(filepath, 'r') as f:
            sample = f.read(4096)
            sniffer = csv.Sniffer()
            try:
                dialect = sniffer.sniff(sample)
                delimiter = dialect.delimiter
            except csv.Error:
                delimiter = ','
        
        return pd.read_csv(
            filepath,
            delimiter=delimiter,
            header=header_row,
            skiprows=list(range(skip_rows)) if skip_rows > 0 else None
        )
    
    def _import_excel(
        self,
        filepath: Path,
        sheet_name: Optional[str],
        header_row: int,
        skip_rows: int
    ) -> pd.DataFrame:
        """Import Excel file."""
        return pd.read_excel(
            filepath,
            sheet_name=sheet_name or 0,
            header=header_row,
            skiprows=list(range(skip_rows)) if skip_rows > 0 else None,
            engine='openpyxl'
        )
    
    def _validate_data(self, df: pd.DataFrame) -> List[str]:
        """Validate imported data and return warnings."""
        warnings = []
        
        # Check for empty columns
        empty_cols = df.columns[df.isna().all()].tolist()
        if empty_cols:
            warnings.append(f"Empty columns detected: {empty_cols}")
        
        # Check for missing values
        missing_count = df.isna().sum().sum()
        if missing_count > 0:
            warnings.append(f"Data contains {missing_count} missing values")
        
        # Check for numeric columns
        numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
        if not numeric_cols:
            warnings.append("No numeric columns detected")
        
        return warnings
    
    def get_excel_sheets(self, filepath: Union[str, Path]) -> List[str]:
        """Get list of sheet names from an Excel file."""
        filepath = Path(filepath)
        if filepath.suffix.lower() not in {'.xlsx', '.xls'}:
            return []
        
        try:
            xl = pd.ExcelFile(filepath, engine='openpyxl')
            return xl.sheet_names
        except Exception:
            return []


@dataclass
class CalculationProfile:
    """A saved calculation profile."""
    name: str
    equation_id: str
    equation_name: str
    category: str
    inputs: Dict[str, Dict[str, Any]]  # {param_name: {value, unit, uncertainty}}
    settings: Dict[str, Any] = field(default_factory=dict)
    created_at: str = ""
    modified_at: str = ""
    notes: str = ""
    
    def __post_init__(self):
        if not self.created_at:
            self.created_at = datetime.now().isoformat()
        self.modified_at = datetime.now().isoformat()


class ProfileManager:
    """
    Manages saving and loading calculation profiles.
    """
    
    def __init__(self, profiles_dir: Optional[Path] = None):
        """
        Initialize profile manager.
        
        Args:
            profiles_dir: Directory for storing profiles
        """
        if profiles_dir is None:
            settings = get_settings()
            profiles_dir = Path(settings.profiles_dir)
        
        self.profiles_dir = Path(profiles_dir)
        self.profiles_dir.mkdir(parents=True, exist_ok=True)
    
    def save_profile(self, profile: CalculationProfile) -> bool:
        """
        Save a profile to disk.
        
        Args:
            profile: Profile to save
        
        Returns:
            True if successful
        """
        profile.modified_at = datetime.now().isoformat()
        filename = self._sanitize_filename(profile.name) + '.json'
        filepath = self.profiles_dir / filename
        
        try:
            with open(filepath, 'w') as f:
                json.dump(asdict(profile), f, indent=2)
            return True
        except Exception as e:
            print(f"Error saving profile: {e}")
            return False
    
    def load_profile(self, name: str) -> Optional[CalculationProfile]:
        """
        Load a profile from disk.
        
        Args:
            name: Profile name
        
        Returns:
            CalculationProfile or None if not found
        """
        filename = self._sanitize_filename(name) + '.json'
        filepath = self.profiles_dir / filename
        
        if not filepath.exists():
            return None
        
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            return CalculationProfile(**data)
        except Exception as e:
            print(f"Error loading profile: {e}")
            return None
    
    def list_profiles(self) -> List[str]:
        """Get list of available profile names."""
        profiles = []
        for filepath in self.profiles_dir.glob('*.json'):
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                profiles.append(data.get('name', filepath.stem))
            except Exception:
                continue
        return sorted(profiles)
    
    def delete_profile(self, name: str) -> bool:
        """Delete a profile."""
        filename = self._sanitize_filename(name) + '.json'
        filepath = self.profiles_dir / filename
        
        if filepath.exists():
            filepath.unlink()
            return True
        return False
    
    def _sanitize_filename(self, name: str) -> str:
        """Sanitize a name for use as filename."""
        # Remove/replace invalid characters
        invalid_chars = '<>:"/\\|?*'
        for char in invalid_chars:
            name = name.replace(char, '_')
        return name.strip()


class PDFExporter:
    """
    Exports calculation results to PDF format.
    """
    
    def __init__(self):
        self.styles = getSampleStyleSheet()
        self._setup_custom_styles()
    
    def _setup_custom_styles(self):
        """Setup custom paragraph styles."""
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=18,
            spaceAfter=20,
            textColor=colors.HexColor('#1a237e')
        ))
        
        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            parent=self.styles['Heading2'],
            fontSize=14,
            spaceBefore=15,
            spaceAfter=10,
            textColor=colors.HexColor('#303f9f')
        ))
        
        self.styles.add(ParagraphStyle(
            name='SubSection',
            parent=self.styles['Heading3'],
            fontSize=12,
            spaceBefore=10,
            spaceAfter=5
        ))
    
    def export(
        self,
        filepath: Union[str, Path],
        title: str,
        equation_name: str,
        inputs: Dict[str, Dict[str, Any]],
        outputs: Dict[str, Dict[str, Any]],
        validation_results: Optional[Dict[str, Any]] = None,
        uncertainty_results: Optional[Dict[str, Any]] = None,
        notes: str = "",
        chart_images: Optional[List[Path]] = None
    ) -> bool:
        """
        Export calculation results to PDF.
        
        Args:
            filepath: Output PDF path
            title: Report title
            equation_name: Name of the equation used
            inputs: Input parameters {name: {value, unit, description}}
            outputs: Output results {name: {value, unit, description}}
            validation_results: Optional validation information
            uncertainty_results: Optional uncertainty analysis results
            notes: Optional notes
            chart_images: Optional list of chart image paths to include
        
        Returns:
            True if successful
        """
        filepath = Path(filepath)
        
        try:
            doc = SimpleDocTemplate(
                str(filepath),
                pagesize=letter,
                rightMargin=0.75*inch,
                leftMargin=0.75*inch,
                topMargin=0.75*inch,
                bottomMargin=0.75*inch
            )
            
            story = []
            
            # Title
            story.append(Paragraph(title, self.styles['CustomTitle']))
            story.append(Paragraph(f"Equation: {equation_name}", self.styles['Normal']))
            story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", self.styles['Normal']))
            story.append(Spacer(1, 20))
            
            # Inputs section
            story.append(Paragraph("Input Parameters", self.styles['SectionHeader']))
            input_data = [['Parameter', 'Value', 'Unit']]
            for name, info in inputs.items():
                input_data.append([
                    name,
                    f"{info.get('value', 'N/A'):.6g}" if isinstance(info.get('value'), (int, float)) else str(info.get('value', 'N/A')),
                    info.get('unit', '')
                ])
            
            input_table = Table(input_data, colWidths=[2.5*inch, 2*inch, 1.5*inch])
            input_table.setStyle(self._get_table_style())
            story.append(input_table)
            story.append(Spacer(1, 20))
            
            # Outputs section
            story.append(Paragraph("Results", self.styles['SectionHeader']))
            output_data = [['Parameter', 'Value', 'Unit']]
            for name, info in outputs.items():
                output_data.append([
                    name,
                    f"{info.get('value', 'N/A'):.6g}" if isinstance(info.get('value'), (int, float)) else str(info.get('value', 'N/A')),
                    info.get('unit', '')
                ])
            
            output_table = Table(output_data, colWidths=[2.5*inch, 2*inch, 1.5*inch])
            output_table.setStyle(self._get_table_style(header_color=colors.HexColor('#1b5e20')))
            story.append(output_table)
            story.append(Spacer(1, 20))
            
            # Validation section
            if validation_results:
                story.append(Paragraph("Validation", self.styles['SectionHeader']))
                
                status = validation_results.get('status', 'Unknown')
                status_color = {
                    'pass': colors.green,
                    'warning': colors.orange,
                    'fail': colors.red
                }.get(status.lower(), colors.gray)
                
                story.append(Paragraph(
                    f"<b>Status:</b> <font color='{status_color}'>{status.upper()}</font>",
                    self.styles['Normal']
                ))
                
                messages = validation_results.get('messages', [])
                for msg in messages:
                    story.append(Paragraph(f"• {msg}", self.styles['Normal']))
                
                story.append(Spacer(1, 20))
            
            # Uncertainty section
            if uncertainty_results:
                story.append(Paragraph("Uncertainty Analysis", self.styles['SectionHeader']))
                
                uncertainty_data = [['Parameter', 'Mean ± Std Dev', '95% CI', 'Unit']]
                for name, info in uncertainty_results.items():
                    uncertainty_data.append([
                        name,
                        f"{info.get('mean', 0):.4g} ± {info.get('std_dev', 0):.4g}",
                        f"[{info.get('ci_lower', 0):.4g}, {info.get('ci_upper', 0):.4g}]",
                        info.get('unit', '')
                    ])
                
                unc_table = Table(uncertainty_data, colWidths=[1.5*inch, 1.75*inch, 1.75*inch, 1*inch])
                unc_table.setStyle(self._get_table_style(header_color=colors.HexColor('#e65100')))
                story.append(unc_table)
                story.append(Spacer(1, 20))
            
            # Charts
            if chart_images:
                story.append(Paragraph("Charts", self.styles['SectionHeader']))
                for img_path in chart_images:
                    if Path(img_path).exists():
                        img = Image(str(img_path), width=5*inch, height=3.5*inch)
                        story.append(img)
                        story.append(Spacer(1, 10))
            
            # Notes
            if notes:
                story.append(Paragraph("Notes", self.styles['SectionHeader']))
                story.append(Paragraph(notes, self.styles['Normal']))
            
            # Build PDF
            doc.build(story)
            return True
        
        except Exception as e:
            print(f"Error exporting PDF: {e}")
            return False
    
    def _get_table_style(self, header_color=None) -> TableStyle:
        """Get standard table style."""
        if header_color is None:
            header_color = colors.HexColor('#1a237e')
        
        return TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), header_color),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 11),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 10),
            ('TOPPADDING', (0, 0), (-1, 0), 10),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#f5f5f5')),
            ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 1), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
        ])
