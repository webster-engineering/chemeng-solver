"""
Tests for equation calculations.
"""

import pytest
from core.units import ureg
from equations.process_control import (
    ZieglerNicholsPID, CohenCoonPID, ControlValveCvLiquid,
    FirstOrderResponse, SecondOrderResponse
)
from equations.fluid_dynamics import DarcyWeisbach, HazenWilliams
from equations.heat_transfer import LMTD


class TestZieglerNichols:
    """Test Ziegler-Nichols PID tuning."""
    
    def test_pid_tuning(self):
        """Test standard PID tuning calculation."""
        eq = ZieglerNicholsPID()
        result = eq.calculate({'Ku': 2.0, 'Pu': (5.0, 'min')})
        
        assert result.success
        assert abs(result.get_output_value('Kp') - 1.2) < 0.01
        assert abs(result.get_output_value('Ti', 'min') - 2.5) < 0.01
        assert abs(result.get_output_value('Td', 'min') - 0.625) < 0.01


class TestControlValve:
    """Test control valve Cv calculations."""
    
    def test_liquid_cv(self):
        """Test liquid Cv calculation."""
        eq = ControlValveCvLiquid()
        result = eq.calculate({
            'Q': (100, 'gpm'),
            'dP': (25, 'psi'),
            'SG': 1.0
        })
        
        assert result.success
        assert abs(result.get_output_value('Cv') - 20.0) < 0.1


class TestFirstOrderResponse:
    """Test first-order response analysis."""
    
    def test_response_times(self):
        """Test time constant calculations."""
        eq = FirstOrderResponse()
        result = eq.calculate({
            'K': 2.0,
            'tau': (10, 'min'),
            'step_size': 1.0
        })
        
        assert result.success
        assert abs(result.get_output_value('t63', 'min') - 10) < 0.1
        assert abs(result.get_output_value('t95', 'min') - 30) < 0.1
        assert abs(result.get_output_value('final_value') - 2.0) < 0.01


class TestSecondOrderResponse:
    """Test second-order response analysis."""
    
    def test_underdamped(self):
        """Test underdamped response."""
        eq = SecondOrderResponse()
        result = eq.calculate({
            'K': 1.0, 'tau': (1, 'min'), 'zeta': 0.5
        })
        
        assert result.success
        assert result.get_output_value('overshoot') > 0


class TestLMTD:
    """Test LMTD calculation."""
    
    def test_counterflow(self):
        """Test LMTD for counterflow heat exchanger."""
        eq = LMTD()
        result = eq.calculate({
            'T_hot_in': (300, 'degF'),
            'T_hot_out': (200, 'degF'),
            'T_cold_in': (100, 'degF'),
            'T_cold_out': (150, 'degF')
        })
        
        assert result.success
        lmtd = result.get_output_value('LMTD')
        # dt1 = 300-150 = 150, dt2 = 200-100 = 100
        # LMTD = (150-100)/ln(150/100) = 50/0.405 = 123.3
        assert abs(lmtd - 123.3) < 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
