"""
Mass Transfer equations - absorption, distillation, extraction.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class McCabeThiele(BaseEquation):
    """McCabe-Thiele method for distillation stages."""
    
    equation_id = "mccabe_thiele"
    name = "McCabe-Thiele Stages"
    category = "Mass Transfer"
    description = "Estimate theoretical stages for binary distillation"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("xD", "Distillate composition (mole fraction)", "dimensionless", "",
                              ParameterType.INPUT, symbol="xD", typical_range=(0, 1)),
            EquationParameter("xB", "Bottoms composition (mole fraction)", "dimensionless", "",
                              ParameterType.INPUT, symbol="xB", typical_range=(0, 1)),
            EquationParameter("xF", "Feed composition (mole fraction)", "dimensionless", "",
                              ParameterType.INPUT, symbol="xF", typical_range=(0, 1)),
            EquationParameter("R", "Reflux ratio (L/D)", "dimensionless", "",
                              ParameterType.INPUT, symbol="R", typical_range=(0.5, 20)),
            EquationParameter("alpha", "Relative volatility", "dimensionless", "",
                              ParameterType.INPUT, symbol="α", typical_range=(1.1, 10)),
            EquationParameter("N_theoretical", "Theoretical stages", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="N"),
        ]
    
    def _calculate(self, xD: pint.Quantity, xB: pint.Quantity, xF: pint.Quantity,
                   R: pint.Quantity, alpha: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        xd = xD.magnitude if hasattr(xD, 'magnitude') else float(xD)
        xb = xB.magnitude if hasattr(xB, 'magnitude') else float(xB)
        xf = xF.magnitude if hasattr(xF, 'magnitude') else float(xF)
        r = R.magnitude if hasattr(R, 'magnitude') else float(R)
        a = alpha.magnitude if hasattr(alpha, 'magnitude') else float(alpha)
        
        # Fenske equation for minimum stages
        n_min = np.log((xd/(1-xd)) * ((1-xb)/xb)) / np.log(a)
        
        # Underwood for minimum reflux (simplified)
        r_min = (1/(a-1)) * (xd/xf - a*(1-xd)/(1-xf))
        
        # Gilliland correlation
        x = (r - r_min) / (r + 1)
        y = 1 - np.exp((1 + 54.4*x) / (11 + 117.2*x) * (x - 1) / np.sqrt(x))
        
        n_theo = (n_min + y) / (1 - y)
        
        return {"N_theoretical": ureg.Quantity(max(n_theo, n_min), "")}


class Kremser(BaseEquation):
    """Kremser equation for absorption/stripping."""
    
    equation_id = "kremser"
    name = "Kremser Equation"
    category = "Mass Transfer"
    description = "Calculate stages for absorption or stripping"
    reference = "Seader & Henley, Separation Process Principles"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("A", "Absorption factor (L/mG)", "dimensionless", "",
                              ParameterType.INPUT, symbol="A", typical_range=(0.5, 5)),
            EquationParameter("y_in", "Inlet gas composition", "dimensionless", "",
                              ParameterType.INPUT, symbol="yᵢₙ"),
            EquationParameter("y_out", "Outlet gas composition", "dimensionless", "",
                              ParameterType.INPUT, symbol="yₒᵤₜ"),
            EquationParameter("y_star", "Equilibrium composition", "dimensionless", "",
                              ParameterType.INPUT, symbol="y*"),
            EquationParameter("N", "Number of stages", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="N"),
        ]
    
    def _calculate(self, A: pint.Quantity, y_in: pint.Quantity, y_out: pint.Quantity,
                   y_star: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        a = A.magnitude if hasattr(A, 'magnitude') else float(A)
        yin = y_in.magnitude if hasattr(y_in, 'magnitude') else float(y_in)
        yout = y_out.magnitude if hasattr(y_out, 'magnitude') else float(y_out)
        ys = y_star.magnitude if hasattr(y_star, 'magnitude') else float(y_star)
        
        if abs(a - 1) < 0.001:
            n = (yin - yout) / (yin - ys)
        else:
            phi = (yin - ys) / (yout - ys)
            n = np.log((phi - 1/a) / (1 - 1/a)) / np.log(a)
        
        return {"N": ureg.Quantity(max(n, 1), "")}


MASS_TRANSFER_EQUATIONS = {
    'mccabe_thiele': McCabeThiele,
    'kremser': Kremser,
}
