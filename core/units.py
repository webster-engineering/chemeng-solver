"""
Comprehensive unit conversion system using pint library.
Supports seamless SI/Imperial conversion with default to Imperial.
"""

import pint
from typing import Union, Tuple, Optional, List
from dataclasses import dataclass


# Create custom unit registry
ureg = pint.UnitRegistry()

# Define custom units commonly used in chemical engineering
ureg.define('gpm = gallon / minute')
ureg.define('scfm = ft**3 / minute')  # Standard cubic feet per minute
ureg.define('acfm = ft**3 / minute')  # Actual cubic feet per minute
ureg.define('barg = bar')  # Gauge pressure (just alias, user handles offset)
ureg.define('psig = psi')  # Gauge pressure (just alias, user handles offset)
ureg.define('BTU_per_hr = BTU / hour')
ureg.define('BTU_per_hr_ft2_degF = BTU / hour / ft**2 / degF')


@dataclass
class UnitCategory:
    """Defines a category of units with common conversions."""
    name: str
    dimension: str
    imperial_default: str
    si_default: str
    common_units: List[str]


# Comprehensive unit categories for chemical engineering
UNIT_CATEGORIES = {
    'length': UnitCategory(
        name='Length',
        dimension='[length]',
        imperial_default='ft',
        si_default='m',
        common_units=['ft', 'm', 'in', 'cm', 'mm', 'mile', 'km']
    ),
    'mass': UnitCategory(
        name='Mass',
        dimension='[mass]',
        imperial_default='lb',
        si_default='kg',
        common_units=['lb', 'kg', 'g', 'ton', 'tonne', 'oz']
    ),
    'time': UnitCategory(
        name='Time',
        dimension='[time]',
        imperial_default='s',
        si_default='s',
        common_units=['s', 'min', 'hr', 'day']
    ),
    'temperature': UnitCategory(
        name='Temperature',
        dimension='[temperature]',
        imperial_default='degF',
        si_default='degC',
        common_units=['degF', 'degC', 'K', 'degR']
    ),
    'temperature_difference': UnitCategory(
        name='Temperature Difference',
        dimension='[temperature]',
        imperial_default='delta_degF',
        si_default='delta_degC',
        common_units=['delta_degF', 'delta_degC', 'K', 'delta_degR']
    ),
    'pressure': UnitCategory(
        name='Pressure',
        dimension='[pressure]',
        imperial_default='psi',
        si_default='bar',
        common_units=['psi', 'psig', 'bar', 'barg', 'kPa', 'MPa', 'atm', 'Pa', 'inH2O', 'mmHg']
    ),
    'volumetric_flow': UnitCategory(
        name='Volumetric Flow Rate',
        dimension='[length] ** 3 / [time]',
        imperial_default='gpm',
        si_default='m**3/hr',
        common_units=['gpm', 'm**3/hr', 'L/min', 'L/s', 'scfm', 'acfm', 'ft**3/s']
    ),
    'mass_flow': UnitCategory(
        name='Mass Flow Rate',
        dimension='[mass] / [time]',
        imperial_default='lb/hr',
        si_default='kg/hr',
        common_units=['lb/hr', 'lb/s', 'kg/hr', 'kg/s', 'ton/hr', 'tonne/hr']
    ),
    'velocity': UnitCategory(
        name='Velocity',
        dimension='[length] / [time]',
        imperial_default='ft/s',
        si_default='m/s',
        common_units=['ft/s', 'm/s', 'ft/min', 'm/min', 'mph', 'km/hr']
    ),
    'density': UnitCategory(
        name='Density',
        dimension='[mass] / [length] ** 3',
        imperial_default='lb/ft**3',
        si_default='kg/m**3',
        common_units=['lb/ft**3', 'kg/m**3', 'lb/gal', 'g/cm**3', 'g/L']
    ),
    'viscosity_dynamic': UnitCategory(
        name='Dynamic Viscosity',
        dimension='[mass] / [length] / [time]',
        imperial_default='cP',
        si_default='Pa*s',
        common_units=['cP', 'Pa*s', 'mPa*s', 'lb/(ft*s)', 'poise']
    ),
    'viscosity_kinematic': UnitCategory(
        name='Kinematic Viscosity',
        dimension='[length] ** 2 / [time]',
        imperial_default='cSt',
        si_default='m**2/s',
        common_units=['cSt', 'm**2/s', 'ft**2/s', 'stokes']
    ),
    'energy': UnitCategory(
        name='Energy',
        dimension='[energy]',
        imperial_default='BTU',
        si_default='kJ',
        common_units=['BTU', 'kJ', 'J', 'cal', 'kcal', 'kWh']
    ),
    'power': UnitCategory(
        name='Power',
        dimension='[power]',
        imperial_default='hp',
        si_default='kW',
        common_units=['hp', 'kW', 'W', 'BTU/hr', 'MW']
    ),
    'heat_transfer_coefficient': UnitCategory(
        name='Heat Transfer Coefficient',
        dimension='[mass] / [time] ** 3 / [temperature]',
        imperial_default='BTU/(hr*ft**2*delta_degF)',
        si_default='W/(m**2*K)',
        common_units=['BTU/(hr*ft**2*delta_degF)', 'W/(m**2*K)', 'kW/(m**2*K)']
    ),
    'thermal_conductivity': UnitCategory(
        name='Thermal Conductivity',
        dimension='[power] / [length] / [temperature]',
        imperial_default='BTU/(hr*ft*delta_degF)',
        si_default='W/(m*K)',
        common_units=['BTU/(hr*ft*delta_degF)', 'W/(m*K)']
    ),
    'specific_heat': UnitCategory(
        name='Specific Heat',
        dimension='[energy] / [mass] / [temperature]',
        imperial_default='BTU/(lb*delta_degF)',
        si_default='kJ/(kg*K)',
        common_units=['BTU/(lb*delta_degF)', 'kJ/(kg*K)', 'J/(g*K)', 'cal/(g*degC)']
    ),
    'area': UnitCategory(
        name='Area',
        dimension='[length] ** 2',
        imperial_default='ft**2',
        si_default='m**2',
        common_units=['ft**2', 'm**2', 'in**2', 'cm**2', 'mm**2']
    ),
    'volume': UnitCategory(
        name='Volume',
        dimension='[length] ** 3',
        imperial_default='gal',
        si_default='L',
        common_units=['gal', 'L', 'ft**3', 'm**3', 'bbl', 'mL']
    ),
    'concentration': UnitCategory(
        name='Concentration',
        dimension='[substance] / [length] ** 3',
        imperial_default='mol/ft**3',
        si_default='mol/L',
        common_units=['mol/L', 'mol/m**3', 'mol/ft**3', 'kmol/m**3']
    ),
    'dimensionless': UnitCategory(
        name='Dimensionless',
        dimension='',
        imperial_default='',
        si_default='',
        common_units=['']  # dimensionless
    ),
}


class UnitConverter:
    """
    Main unit conversion class for chemical engineering calculations.
    """
    
    def __init__(self, default_system: str = 'imperial'):
        """
        Initialize converter with default unit system.
        
        Args:
            default_system: 'imperial' or 'si'
        """
        self.default_system = default_system
        self._ureg = ureg
    
    @property
    def Q_(self):
        """Shortcut to create quantities."""
        return self._ureg.Quantity
    
    def parse(self, value: Union[str, float], unit: str = '') -> pint.Quantity:
        """
        Parse a value with unit into a pint Quantity.
        
        Args:
            value: Numeric value or string like "100 gpm"
            unit: Unit string (if value is numeric)
        
        Returns:
            pint Quantity
        """
        if isinstance(value, str):
            return self._ureg.parse_expression(value)
        elif unit:
            return self.Q_(value, unit)
        else:
            return self.Q_(value, '')  # dimensionless
    
    def convert(self, quantity: pint.Quantity, target_unit: str) -> pint.Quantity:
        """
        Convert a quantity to a target unit.
        
        Args:
            quantity: pint Quantity to convert
            target_unit: Target unit string
        
        Returns:
            Converted pint Quantity
        """
        return quantity.to(target_unit)
    
    def convert_value(self, value: float, from_unit: str, to_unit: str) -> float:
        """
        Convert a numeric value between units.
        
        Args:
            value: Numeric value
            from_unit: Source unit
            to_unit: Target unit
        
        Returns:
            Converted numeric value
        """
        q = self.Q_(value, from_unit)
        return q.to(to_unit).magnitude
    
    def get_default_unit(self, category: str) -> str:
        """
        Get the default unit for a category based on current system.
        
        Args:
            category: Unit category name
        
        Returns:
            Default unit string
        """
        if category not in UNIT_CATEGORIES:
            return ''
        cat = UNIT_CATEGORIES[category]
        if self.default_system == 'imperial':
            return cat.imperial_default
        return cat.si_default
    
    def get_common_units(self, category: str) -> List[str]:
        """
        Get list of common units for a category.
        
        Args:
            category: Unit category name
        
        Returns:
            List of unit strings
        """
        if category not in UNIT_CATEGORIES:
            return []
        return UNIT_CATEGORIES[category].common_units
    
    def check_dimensionality(self, quantity: pint.Quantity, expected_dimension: str) -> bool:
        """
        Check if a quantity has the expected dimensionality.
        
        Args:
            quantity: pint Quantity to check
            expected_dimension: Expected dimension string like '[length] / [time]'
        
        Returns:
            True if dimensions match
        """
        if not expected_dimension:  # dimensionless
            return quantity.dimensionless
        
        expected = self._ureg.parse_expression(expected_dimension).dimensionality
        return quantity.dimensionality == expected
    
    def get_dimension_string(self, quantity: pint.Quantity) -> str:
        """
        Get a human-readable dimension string for a quantity.
        
        Args:
            quantity: pint Quantity
        
        Returns:
            Dimension string like "[length] / [time]"
        """
        return str(quantity.dimensionality)
    
    def to_system(self, quantity: pint.Quantity, system: str = None) -> pint.Quantity:
        """
        Convert a quantity to the default unit of its dimension in the specified system.
        
        Args:
            quantity: pint Quantity
            system: 'imperial' or 'si' (defaults to converter's default)
        
        Returns:
            Converted quantity
        """
        if system is None:
            system = self.default_system
        
        # Find matching category
        for cat_name, cat in UNIT_CATEGORIES.items():
            if cat.dimension and self.check_dimensionality(quantity, cat.dimension):
                target = cat.imperial_default if system == 'imperial' else cat.si_default
                try:
                    return quantity.to(target)
                except pint.DimensionalityError:
                    pass
        
        return quantity  # Return unchanged if no matching category
    
    def format_quantity(self, quantity: pint.Quantity, precision: int = 4) -> str:
        """
        Format a quantity as a string with appropriate precision.
        
        Args:
            quantity: pint Quantity
            precision: Number of significant figures
        
        Returns:
            Formatted string
        """
        return f"{quantity.magnitude:.{precision}g} {quantity.units:~P}"


# Global converter instance
_converter: Optional[UnitConverter] = None


def get_converter() -> UnitConverter:
    """Get the global unit converter instance."""
    global _converter
    if _converter is None:
        from config import get_settings
        settings = get_settings()
        _converter = UnitConverter(default_system=settings.default_unit_system)
    return _converter


def set_default_system(system: str) -> None:
    """Set the default unit system for the global converter."""
    global _converter
    if _converter is not None:
        _converter.default_system = system
