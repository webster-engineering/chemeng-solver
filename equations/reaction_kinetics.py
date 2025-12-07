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
