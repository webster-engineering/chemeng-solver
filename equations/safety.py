"""
Safety and Relief equations - relief valve, flare, overpressure.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class ReliefValveLiquid(BaseEquation):
    """Size relief valve for liquid service per API 520."""
    
    equation_id = "relief_liquid"
    name = "Relief Valve Sizing (Liquid)"
    category = "Safety"
    description = "Calculate required orifice area for liquid relief per API 520"
    reference = "API 520 Part I"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Required relief flow", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("P1", "Set pressure + overpressure", "pressure", "psig",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Back pressure", "pressure", "psig",
                              ParameterType.INPUT, symbol="P₂"),
            EquationParameter("SG", "Specific gravity", "dimensionless", "",
                              ParameterType.INPUT, symbol="SG"),
            EquationParameter("Kd", "Discharge coefficient", "dimensionless", "",
                              ParameterType.INPUT, symbol="Kd",
                              typical_range=(0.62, 0.65), tooltip="0.65 typical for preliminary"),
            EquationParameter("Kw", "Back pressure correction", "dimensionless", "",
                              ParameterType.INPUT, symbol="Kw",
                              typical_range=(0.9, 1.0)),
            EquationParameter("A", "Required orifice area", "area", "in**2",
                              ParameterType.OUTPUT, symbol="A"),
            EquationParameter("orifice", "Standard orifice letter", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Orifice"),
        ]
    
    # Standard API orifice sizes (in²)
    ORIFICES = [
        ('D', 0.110), ('E', 0.196), ('F', 0.307), ('G', 0.503),
        ('H', 0.785), ('J', 1.287), ('K', 1.838), ('L', 2.853),
        ('M', 3.60), ('N', 4.34), ('P', 6.38), ('Q', 11.05),
        ('R', 16.0), ('T', 26.0)
    ]
    
    def _calculate(self, Q: pint.Quantity, P1: pint.Quantity, P2: pint.Quantity,
                   SG: pint.Quantity, Kd: pint.Quantity, Kw: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gpm').magnitude
        p1 = P1.to('psig').magnitude + 14.7  # Convert to psia
        p2 = P2.to('psig').magnitude + 14.7
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        kd = Kd.magnitude if hasattr(Kd, 'magnitude') else float(Kd)
        kw = Kw.magnitude if hasattr(Kw, 'magnitude') else float(Kw)
        
        dp = p1 - p2
        
        # API 520 liquid formula: A = Q / (38 * Kd * Kw * sqrt(dP/SG))
        a = q / (38 * kd * kw * np.sqrt(dp / sg))
        
        # Find next standard orifice
        orifice = 'T'  # Default to largest
        for letter, size in self.ORIFICES:
            if size >= a:
                orifice = letter
                break
        
        return {
            "A": ureg.Quantity(a, "in**2"),
            "orifice": ureg.Quantity(ord(orifice) - ord('A'), "")  # Numeric code
        }


class ReliefValveGas(BaseEquation):
    """Size relief valve for gas/vapor service per API 520."""
    
    equation_id = "relief_gas"
    name = "Relief Valve Sizing (Gas)"
    category = "Safety"
    description = "Calculate required orifice area for gas relief per API 520"
    reference = "API 520 Part I"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("W", "Required relief flow", "mass_flow", "lb/hr",
                              ParameterType.INPUT, symbol="W"),
            EquationParameter("P1", "Set pressure + overpressure (abs)", "pressure", "psia",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("T", "Relieving temperature", "temperature", "degF",
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("MW", "Molecular weight", "dimensionless", "",
                              ParameterType.INPUT, symbol="MW"),
            EquationParameter("k", "Specific heat ratio Cp/Cv", "dimensionless", "",
                              ParameterType.INPUT, symbol="k"),
            EquationParameter("Z", "Compressibility factor", "dimensionless", "",
                              ParameterType.INPUT, symbol="Z", typical_range=(0.9, 1.0)),
            EquationParameter("Kd", "Discharge coefficient", "dimensionless", "",
                              ParameterType.INPUT, symbol="Kd", typical_range=(0.95, 0.975)),
            EquationParameter("A", "Required orifice area", "area", "in**2",
                              ParameterType.OUTPUT, symbol="A"),
        ]
    
    def _calculate(self, W: pint.Quantity, P1: pint.Quantity, T: pint.Quantity,
                   MW: pint.Quantity, k: pint.Quantity, Z: pint.Quantity,
                   Kd: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        w = W.to('lb/hr').magnitude
        p1 = P1.to('psia').magnitude
        t = T.to('degR').magnitude
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        k_val = k.magnitude if hasattr(k, 'magnitude') else float(k)
        z = Z.magnitude if hasattr(Z, 'magnitude') else float(Z)
        kd = Kd.magnitude if hasattr(Kd, 'magnitude') else float(Kd)
        
        # C coefficient
        c = 520 * np.sqrt(k_val * (2/(k_val+1))**((k_val+1)/(k_val-1)))
        
        # API 520: A = W * sqrt(T*Z) / (C * Kd * P1 * Kb * sqrt(MW))
        kb = 1.0  # Back pressure factor for conventional valve
        a = w * np.sqrt(t * z) / (c * kd * p1 * kb * np.sqrt(mw))
        
        return {"A": ureg.Quantity(a, "in**2")}


class FireCaseHeatInput(BaseEquation):
    """Calculate heat input for fire case relief sizing."""
    
    equation_id = "fire_case"
    name = "Fire Case Heat Input"
    category = "Safety"
    description = "Calculate heat absorption from pool fire per API 521"
    reference = "API 521"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("A_wet", "Wetted surface area exposed to fire", "area", "ft**2",
                              ParameterType.INPUT, symbol="Awet",
                              tooltip="Area below liquid level, up to 25 ft height"),
            EquationParameter("F", "Environmental factor", "dimensionless", "",
                              ParameterType.INPUT, symbol="F",
                              tooltip="1.0=bare, 0.3=insulated, 0.15=water spray"),
            EquationParameter("Q", "Heat absorbed", "power", "BTU/hr",
                              ParameterType.OUTPUT, symbol="Q"),
            EquationParameter("W_relief", "Approximate relief rate (for water)", "mass_flow", "lb/hr",
                              ParameterType.OUTPUT, symbol="W"),
        ]
    
    def _calculate(self, A_wet: pint.Quantity, F: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        a = A_wet.to('ft**2').magnitude
        f = F.magnitude if hasattr(F, 'magnitude') else float(F)
        
        # API 521 formula: Q = 21000 * F * A^0.82 (for A < 200 ft²)
        # Q = 34500 * F * A^0.82 for A >= 200 ft² (with drainage)
        if a < 200:
            q = 21000 * f * a**0.82
        else:
            q = 34500 * f * a**0.82
        
        # Approximate relief rate using water latent heat
        latent_heat = 970  # BTU/lb for water
        w_relief = q / latent_heat
        
        return {
            "Q": ureg.Quantity(q, "BTU/hr"),
            "W_relief": ureg.Quantity(w_relief, "lb/hr")
        }


class FlareHeaderSizing(BaseEquation):
    """Size flare header for gas flow."""
    
    equation_id = "flare_header"
    name = "Flare Header Sizing"
    category = "Safety"
    description = "Calculate flare header diameter for given flow and pressure drop"
    reference = "API 521"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("W", "Mass flow rate", "mass_flow", "lb/hr",
                              ParameterType.INPUT, symbol="W"),
            EquationParameter("P1", "Inlet pressure", "pressure", "psig",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("T", "Temperature", "temperature", "degF",
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("MW", "Average molecular weight", "dimensionless", "",
                              ParameterType.INPUT, symbol="MW"),
            EquationParameter("L", "Header length", "length", "ft",
                              ParameterType.INPUT, symbol="L"),
            EquationParameter("dP_allow", "Allowable pressure drop", "pressure", "psi",
                              ParameterType.INPUT, symbol="ΔP"),
            EquationParameter("D", "Required diameter", "length", "in",
                              ParameterType.OUTPUT, symbol="D"),
            EquationParameter("velocity", "Gas velocity", "velocity", "ft/s",
                              ParameterType.OUTPUT, symbol="V"),
        ]
    
    def _calculate(self, W: pint.Quantity, P1: pint.Quantity, T: pint.Quantity,
                   MW: pint.Quantity, L: pint.Quantity, dP_allow: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        w = W.to('lb/hr').magnitude
        p1 = P1.to('psia').magnitude + 14.7
        t = T.to('degR').magnitude
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        length = L.to('ft').magnitude
        dp = dP_allow.to('psi').magnitude
        
        # Gas density
        rho = p1 * mw / (10.73 * t)  # lb/ft³
        
        # Simplified flow equation solving for D
        # Using Darcy with f=0.02 approximation
        f = 0.02
        w_lbs = w / 3600  # lb/s
        
        # Iterative solution for D
        d = 6  # Initial guess (inches)
        for _ in range(10):
            d_ft = d / 12
            a = np.pi * (d_ft/2)**2
            v = w_lbs / (rho * a)
            dp_calc = f * length * rho * v**2 / (2 * 32.174 * d_ft) / 144
            d = d * (dp_calc / dp)**0.2  # Adjust
        
        velocity = v
        
        return {
            "D": ureg.Quantity(d, "in"),
            "velocity": ureg.Quantity(velocity, "ft/s")
        }


SAFETY_EQUATIONS = {
    'relief_liquid': ReliefValveLiquid,
    'relief_gas': ReliefValveGas,
    'fire_case': FireCaseHeatInput,
    'flare_header': FlareHeaderSizing,
}
