"""
Reaction Kinetics equations - rate equations, reactor sizing.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class Arrhenius(BaseEquation):
    """Arrhenius equation for temperature dependence of rate constant."""
    
    equation_id = "arrhenius"
    name = "Arrhenius Equation"
    category = "Reaction Kinetics"
    description = "Calculate rate constant at given temperature"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("A", "Pre-exponential factor", "dimensionless", "1/s",
                              ParameterType.INPUT, symbol="A"),
            EquationParameter("Ea", "Activation energy", "energy", "kJ/mol",
                              ParameterType.INPUT, symbol="Eₐ"),
            EquationParameter("T", "Temperature", "temperature", "degF",
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("k", "Rate constant", "dimensionless", "1/s",
                              ParameterType.OUTPUT, symbol="k"),
        ]
    
    def _calculate(self, A: pint.Quantity, Ea: pint.Quantity,
                   T: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        a_val = A.magnitude if hasattr(A, 'magnitude') else float(A)
        ea = Ea.to('J/mol').magnitude
        t = T.to('K').magnitude
        
        R = 8.314  # J/(mol*K)
        k = a_val * np.exp(-ea / (R * t))
        
        return {"k": ureg.Quantity(k, "1/s")}


class CSTRDesign(BaseEquation):
    """CSTR volume calculation for first-order reaction."""
    
    equation_id = "cstr_design"
    name = "CSTR Design (First Order)"
    category = "Reaction Kinetics"
    description = "Calculate CSTR volume for desired conversion"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("F_A0", "Inlet molar flow rate", "mass_flow", "mol/s",
                              ParameterType.INPUT, symbol="FA₀"),
            EquationParameter("C_A0", "Inlet concentration", "concentration", "mol/L",
                              ParameterType.INPUT, symbol="CA₀"),
            EquationParameter("X", "Desired conversion", "dimensionless", "",
                              ParameterType.INPUT, symbol="X", typical_range=(0, 0.99)),
            EquationParameter("k", "Rate constant", "dimensionless", "1/s",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("V", "Reactor volume", "volume", "L",
                              ParameterType.OUTPUT, symbol="V"),
            EquationParameter("tau", "Space time", "time", "s",
                              ParameterType.OUTPUT, symbol="τ"),
        ]
    
    def _calculate(self, F_A0: pint.Quantity, C_A0: pint.Quantity,
                   X: pint.Quantity, k: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        fa0 = F_A0.to('mol/s').magnitude
        ca0 = C_A0.to('mol/L').magnitude
        x = X.magnitude if hasattr(X, 'magnitude') else float(X)
        k_val = k.to('1/s').magnitude
        
        v0 = fa0 / ca0  # volumetric flow rate L/s
        tau = x / (k_val * (1 - x))
        v = v0 * tau
        
        return {
            "V": ureg.Quantity(v, "L"),
            "tau": ureg.Quantity(tau, "s")
        }


class PFRDesign(BaseEquation):
    """PFR volume calculation for first-order reaction."""
    
    equation_id = "pfr_design"
    name = "PFR Design (First Order)"
    category = "Reaction Kinetics"
    description = "Calculate PFR volume for desired conversion"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("F_A0", "Inlet molar flow rate", "mass_flow", "mol/s",
                              ParameterType.INPUT, symbol="FA₀"),
            EquationParameter("C_A0", "Inlet concentration", "concentration", "mol/L",
                              ParameterType.INPUT, symbol="CA₀"),
            EquationParameter("X", "Desired conversion", "dimensionless", "",
                              ParameterType.INPUT, symbol="X", typical_range=(0, 0.99)),
            EquationParameter("k", "Rate constant", "dimensionless", "1/s",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("V", "Reactor volume", "volume", "L",
                              ParameterType.OUTPUT, symbol="V"),
            EquationParameter("tau", "Space time", "time", "s",
                              ParameterType.OUTPUT, symbol="τ"),
        ]
    
    def _calculate(self, F_A0: pint.Quantity, C_A0: pint.Quantity,
                   X: pint.Quantity, k: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        fa0 = F_A0.to('mol/s').magnitude
        ca0 = C_A0.to('mol/L').magnitude
        x = X.magnitude if hasattr(X, 'magnitude') else float(X)
        k_val = k.to('1/s').magnitude
        
        v0 = fa0 / ca0
        tau = -np.log(1 - x) / k_val
        v = v0 * tau
        
        return {
            "V": ureg.Quantity(v, "L"),
            "tau": ureg.Quantity(tau, "s")
        }


REACTION_KINETICS_EQUATIONS = {
    'arrhenius': Arrhenius,
    'cstr_design': CSTRDesign,
    'pfr_design': PFRDesign,
}


class HalfLife(BaseEquation):
    """Half-life calculation for first-order reactions."""
    
    equation_id = "half_life"
    name = "Reaction Half-Life"
    category = "Reaction Kinetics"
    description = "Calculate half-life for first or second order reactions"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("k", "Rate constant", "1/s", "1/min",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("order", "Reaction order (1 or 2)", "dimensionless", "",
                              ParameterType.INPUT, symbol="n", typical_range=(1, 2)),
            EquationParameter("C_A0", "Initial concentration (for 2nd order)", "mol/L", "",
                              ParameterType.INPUT, symbol="CA₀", required=False),
            EquationParameter("t_half", "Half-life", "s", "min",
                              ParameterType.OUTPUT, symbol="t½"),
        ]
    
    def _calculate(self, k: pint.Quantity, order: pint.Quantity, 
                   C_A0: pint.Quantity = None, **kwargs) -> Dict[str, pint.Quantity]:
        k_val = k.to('1/s').magnitude
        n = int(order.magnitude if hasattr(order, 'magnitude') else float(order))
        
        if n == 1:
            t_half = np.log(2) / k_val
        elif n == 2:
            ca0 = C_A0.to('mol/L').magnitude if C_A0 else 1.0
            t_half = 1 / (k_val * ca0)
        else:
            t_half = np.log(2) / k_val  # default to first order
        
        return {"t_half": ureg.Quantity(t_half, "s")}


class BatchReactor(BaseEquation):
    """Batch reactor time calculation."""
    
    equation_id = "batch_reactor"
    name = "Batch Reactor Time"
    category = "Reaction Kinetics"
    description = "Calculate batch time for desired conversion"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("C_A0", "Initial concentration", "mol/L", "",
                              ParameterType.INPUT, symbol="CA₀"),
            EquationParameter("X", "Desired conversion", "dimensionless", "",
                              ParameterType.INPUT, symbol="X", typical_range=(0, 0.99)),
            EquationParameter("k", "Rate constant", "1/s", "1/min",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("t_batch", "Batch time", "s", "min",
                              ParameterType.OUTPUT, symbol="t"),
            EquationParameter("C_A", "Final concentration", "mol/L", "",
                              ParameterType.OUTPUT, symbol="CA"),
        ]
    
    def _calculate(self, C_A0: pint.Quantity, X: pint.Quantity, 
                   k: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        ca0 = C_A0.to('mol/L').magnitude
        x = X.magnitude if hasattr(X, 'magnitude') else float(X)
        k_val = k.to('1/s').magnitude
        
        # First order: t = -ln(1-X)/k
        t = -np.log(1 - x) / k_val
        ca = ca0 * (1 - x)
        
        return {
            "t_batch": ureg.Quantity(t, "s"),
            "C_A": ureg.Quantity(ca, "mol/L")
        }


class SecondOrderReaction(BaseEquation):
    """Second order reaction conversion."""
    
    equation_id = "second_order"
    name = "Second Order Reaction"
    category = "Reaction Kinetics"
    description = "Calculate conversion for 2nd order irreversible reaction (2A → P)"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("C_A0", "Initial concentration", "mol/L", "",
                              ParameterType.INPUT, symbol="CA₀"),
            EquationParameter("k", "Rate constant", "L/(mol*s)", "",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("t", "Reaction time", "s", "min",
                              ParameterType.INPUT, symbol="t"),
            EquationParameter("X", "Conversion", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="X"),
            EquationParameter("C_A", "Final concentration", "mol/L", "",
                              ParameterType.OUTPUT, symbol="CA"),
        ]
    
    def _calculate(self, C_A0: pint.Quantity, k: pint.Quantity, 
                   t: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        ca0 = C_A0.to('mol/L').magnitude
        k_val = k.magnitude if hasattr(k, 'magnitude') else float(k)
        time = t.to('s').magnitude
        
        # For 2nd order: 1/CA - 1/CA0 = kt
        ca = 1 / (1/ca0 + k_val * time)
        x = (ca0 - ca) / ca0
        
        return {
            "X": ureg.Quantity(x, ""),
            "C_A": ureg.Quantity(ca, "mol/L")
        }


class CatalystDeactivation(BaseEquation):
    """Catalyst activity decay with time-on-stream."""
    
    equation_id = "catalyst_deactivation"
    name = "Catalyst Deactivation"
    category = "Reaction Kinetics"
    description = "First-order catalyst decay: a = exp(-kd*t)"
    reference = "Fogler, Elements of Chemical Reaction Engineering"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("kd", "Deactivation rate constant", "1/hr", "1/day",
                              ParameterType.INPUT, symbol="kd"),
            EquationParameter("t", "Time on stream", "hr", "day",
                              ParameterType.INPUT, symbol="t"),
            EquationParameter("a0", "Initial activity (usually 1)", "dimensionless", "",
                              ParameterType.INPUT, symbol="a₀", typical_range=(0, 1)),
            EquationParameter("a", "Activity at time t", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="a"),
            EquationParameter("t_half", "Catalyst half-life", "hr", "day",
                              ParameterType.OUTPUT, symbol="t½"),
        ]
    
    def _calculate(self, kd: pint.Quantity, t: pint.Quantity, 
                   a0: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        kd_val = kd.to('1/hr').magnitude
        time = t.to('hr').magnitude
        a0_val = a0.magnitude if hasattr(a0, 'magnitude') else float(a0)
        
        a = a0_val * np.exp(-kd_val * time)
        t_half = np.log(2) / kd_val
        
        return {
            "a": ureg.Quantity(a, ""),
            "t_half": ureg.Quantity(t_half, "hr")
        }


# Update registry
REACTION_KINETICS_EQUATIONS = {
    'arrhenius': Arrhenius,
    'cstr_design': CSTRDesign,
    'pfr_design': PFRDesign,
    'half_life': HalfLife,
    'batch_reactor': BatchReactor,
    'second_order': SecondOrderReaction,
    'catalyst_deactivation': CatalystDeactivation,
}
