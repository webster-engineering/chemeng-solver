"""
Tests for unit conversion system.
"""

import pytest
from core.units import UnitConverter, ureg, UNIT_CATEGORIES


class TestUnitConverter:
    """Test the UnitConverter class."""
    
    def setup_method(self):
        self.converter = UnitConverter(default_system='imperial')
    
    def test_parse_value_with_unit(self):
        """Test parsing a value with unit string."""
        q = self.converter.parse("100 gpm")
        assert q.magnitude == 100
    
    def test_parse_numeric_with_unit(self):
        """Test parsing numeric value with separate unit."""
        q = self.converter.parse(50, "psi")
        assert q.magnitude == 50
    
    def test_convert_temperature(self):
        """Test temperature conversion."""
        q = ureg.Quantity(100, "degF")
        result = self.converter.convert(q, "degC")
        assert abs(result.magnitude - 37.78) < 0.1
    
    def test_convert_pressure(self):
        """Test pressure conversion."""
        q = ureg.Quantity(14.7, "psi")
        result = self.converter.convert(q, "bar")
        assert abs(result.magnitude - 1.013) < 0.01
    
    def test_convert_flow(self):
        """Test volumetric flow conversion."""
        q = ureg.Quantity(100, "gpm")
        result = self.converter.convert(q, "L/min")
        assert abs(result.magnitude - 378.5) < 1.0
    
    def test_convert_value(self):
        """Test simple value conversion."""
        result = self.converter.convert_value(1.0, "ft", "m")
        assert abs(result - 0.3048) < 0.001
    
    def test_get_default_unit_imperial(self):
        """Test getting default unit for imperial system."""
        self.converter.default_system = 'imperial'
        assert self.converter.get_default_unit('pressure') == 'psi'
        assert self.converter.get_default_unit('temperature') == 'degF'
    
    def test_get_default_unit_si(self):
        """Test getting default unit for SI system."""
        self.converter.default_system = 'si'
        assert self.converter.get_default_unit('pressure') == 'bar'
        assert self.converter.get_default_unit('temperature') == 'degC'
    
    def test_common_units(self):
        """Test getting common units for a category."""
        units = self.converter.get_common_units('pressure')
        assert 'psi' in units
        assert 'bar' in units
        assert 'kPa' in units


class TestUnitCategories:
    """Test unit category definitions."""
    
    def test_all_categories_defined(self):
        """Ensure all expected categories exist."""
        expected = ['length', 'mass', 'time', 'temperature', 'pressure',
                    'volumetric_flow', 'mass_flow', 'density', 'viscosity_dynamic']
        for cat in expected:
            assert cat in UNIT_CATEGORIES
    
    def test_categories_have_defaults(self):
        """Check all categories have imperial and SI defaults."""
        for name, cat in UNIT_CATEGORIES.items():
            if cat.dimension:  # Skip dimensionless
                assert cat.imperial_default, f"{name} missing imperial default"
                assert cat.si_default, f"{name} missing SI default"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
