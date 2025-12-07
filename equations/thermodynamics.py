"""
Thermodynamics equations - ideal gas, flash, VLE calculations.
Each equation includes detailed variable descriptions and typical units.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class IdealGas(BaseEquation):
    """Ideal gas law calculations."""
    
    equation_id = "ideal_gas"
    name = "Ideal Gas Law (PV=nRT)"
    category = "Thermodynamics"
    description = "Calculate any one unknown from P, V, n, T using PV = nRT"
    reference = "Standard thermodynamics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("n", "Moles of gas - amount of substance", 
                              "dimensionless", "mol", ParameterType.INPUT, symbol="n",
                              tooltip="n = mass / molecular weight"),
            EquationParameter("T", "Absolute temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="T",
                              tooltip="Must use absolute scale internally (R or K)"),
            EquationParameter("P", "Absolute pressure", "pressure", "psia", 
                              ParameterType.INPUT, symbol="P",
                              tooltip="Absolute, not gauge pressure"),
            EquationParameter("V", "Volume occupied by gas", "volume", "ft**3", 
                              ParameterType.OUTPUT, symbol="V"),
        ]
    
    def _calculate(self, n: pint.Quantity, T: pint.Quantity,
                   P: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        n_val = n.to('mol').magnitude
        t = T.to('K').magnitude
        p = P.to('Pa').magnitude
        
        R = 8.314  # J/(mol·K)
        v = n_val * R * t / p  # m³
        v_ft3 = v * 35.3147
        
        return {"V": ureg.Quantity(v_ft3, "ft**3")}


class AntoineVaporPressure(BaseEquation):
    """Antoine equation for vapor pressure."""
    
    equation_id = "antoine"
    name = "Antoine Vapor Pressure"
    category = "Thermodynamics"
    description = "Calculate saturation pressure using Antoine equation: log₁₀(P) = A - B/(C+T)"
    reference = "NIST Chemistry WebBook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("A", "Antoine A coefficient - compound-specific constant", 
                              "dimensionless", "", ParameterType.INPUT, symbol="A",
                              tooltip="Look up in Antoine tables (NIST, Perry's)"),
            EquationParameter("B", "Antoine B coefficient", "dimensionless", "", 
                              ParameterType.INPUT, symbol="B"),
            EquationParameter("C", "Antoine C coefficient", "dimensionless", "", 
                              ParameterType.INPUT, symbol="C"),
            EquationParameter("T", "Temperature - must match units used for coefficients", 
                              "temperature", "degC", ParameterType.INPUT, symbol="T",
                              tooltip="Most Antoine tables use °C"),
            EquationParameter("P_vap", "Vapor (saturation) pressure", "pressure", "mmHg", 
                              ParameterType.OUTPUT, symbol="Pᵛᵃᵖ",
                              tooltip="Pressure at which liquid and vapor are in equilibrium"),
        ]
    
    def _calculate(self, A: pint.Quantity, B: pint.Quantity, C: pint.Quantity,
                   T: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        a = A.magnitude if hasattr(A, 'magnitude') else float(A)
        b = B.magnitude if hasattr(B, 'magnitude') else float(B)
        c = C.magnitude if hasattr(C, 'magnitude') else float(C)
        t = T.to('degC').magnitude
        
        log_p = a - b / (c + t)
        p = 10 ** log_p
        
        return {"P_vap": ureg.Quantity(p, "mmHg")}


class RaoultsLaw(BaseEquation):
    """Raoult's law for ideal VLE."""
    
    equation_id = "raoults_law"
    name = "Raoult's Law (Ideal VLE)"
    category = "Thermodynamics"
    description = "Calculate vapor composition for ideal solution: yᵢ·P = xᵢ·Pᵢˢᵃᵗ"
    reference = "Smith, Van Ness & Abbott"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x", "Liquid mole fraction of component", 
                              "dimensionless", "", ParameterType.INPUT, symbol="x", 
                              typical_range=(0, 1),
                              tooltip="Moles of component / total moles in liquid"),
            EquationParameter("P_sat", "Saturation pressure of pure component at T", 
                              "pressure", "psia", ParameterType.INPUT, symbol="Pˢᵃᵗ",
                              tooltip="From Antoine equation or tables"),
            EquationParameter("P_total", "Total system pressure", "pressure", "psia", 
                              ParameterType.INPUT, symbol="P"),
            EquationParameter("y", "Vapor mole fraction of component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y", 
                              typical_range=(0, 1)),
            EquationParameter("K", "K-value (equilibrium ratio) - y/x", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="K",
                              tooltip="K > 1: component concentrates in vapor"),
        ]
    
    def _calculate(self, x: pint.Quantity, P_sat: pint.Quantity,
                   P_total: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        x_val = x.magnitude if hasattr(x, 'magnitude') else float(x)
        psat = P_sat.to('psia').magnitude
        ptot = P_total.to('psia').magnitude
        
        k = psat / ptot
        y = x_val * k
        
        return {
            "y": ureg.Quantity(min(y, 1.0), ""),
            "K": ureg.Quantity(k, "")
        }


class ClausiusClapeyron(BaseEquation):
    """Clausius-Clapeyron equation for vapor pressure vs temperature."""
    
    equation_id = "clausius_clapeyron"
    name = "Clausius-Clapeyron Equation"
    category = "Thermodynamics"
    description = "Estimate vapor pressure at T2 given vapor pressure at T1 and heat of vaporization"
    reference = "Standard thermodynamics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Vapor pressure at temperature T1", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₁"),
            EquationParameter("T1", "Reference temperature (absolute)", 
                              "temperature", "degF", ParameterType.INPUT, symbol="T₁"),
            EquationParameter("T2", "Target temperature (absolute)", 
                              "temperature", "degF", ParameterType.INPUT, symbol="T₂"),
            EquationParameter("Hvap", "Heat of vaporization - latent heat", 
                              "energy", "BTU/lb", ParameterType.INPUT, symbol="ΔHᵥₐₚ",
                              tooltip="Water≈970 BTU/lb, varies with temperature"),
            EquationParameter("MW", "Molecular weight of substance", 
                              "dimensionless", "", ParameterType.INPUT, symbol="MW",
                              tooltip="g/mol, e.g., Water=18, Ethanol=46"),
            EquationParameter("P2", "Vapor pressure at temperature T2", 
                              "pressure", "psia", ParameterType.OUTPUT, symbol="P₂"),
        ]
    
    def _calculate(self, P1: pint.Quantity, T1: pint.Quantity, T2: pint.Quantity,
                   Hvap: pint.Quantity, MW: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('Pa').magnitude
        t1 = T1.to('K').magnitude
        t2 = T2.to('K').magnitude
        hvap = Hvap.to('J/kg').magnitude
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        
        R = 8.314  # J/(mol·K)
        hvap_molar = hvap * mw / 1000  # J/mol
        
        ln_p2_p1 = (hvap_molar / R) * (1/t1 - 1/t2)
        p2 = p1 * np.exp(ln_p2_p1)
        
        return {"P2": ureg.Quantity(p2 / 6894.76, "psia")}


class FlashCalculation(BaseEquation):
    """Simple flash drum vapor fraction calculation."""
    
    equation_id = "flash"
    name = "Flash Drum Vapor Fraction"
    category = "Thermodynamics"
    description = "Calculate vapor fraction in a flash drum using Rachford-Rice"
    reference = "Seader & Henley, Separation Process Principles"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("z", "Feed composition (mole fraction of light component)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="z",
                              typical_range=(0, 1)),
            EquationParameter("K1", "K-value of light component at drum conditions", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K₁",
                              tooltip="K = y/x = Psat/P for ideal mixtures"),
            EquationParameter("K2", "K-value of heavy component at drum conditions", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K₂"),
            EquationParameter("V_F", "Vapor fraction - moles vapor per mole feed", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="V/F",
                              typical_range=(0, 1)),
            EquationParameter("y1", "Vapor composition of light component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y₁"),
            EquationParameter("x1", "Liquid composition of light component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="x₁"),
        ]
    
    def _calculate(self, z: pint.Quantity, K1: pint.Quantity, K2: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        z1 = z.magnitude if hasattr(z, 'magnitude') else float(z)
        z2 = 1 - z1
        k1 = K1.magnitude if hasattr(K1, 'magnitude') else float(K1)
        k2 = K2.magnitude if hasattr(K2, 'magnitude') else float(K2)
        
        # Rachford-Rice: solve sum(zi(Ki-1)/(1+V/F(Ki-1))) = 0
        # For binary: V/F = (z1(K1-1) + z2(K2-1)) / ((K1-1)(K2-1))... simplified
        # Use bisection for robustness
        def rr(vf):
            return z1*(k1-1)/(1+vf*(k1-1)) + z2*(k2-1)/(1+vf*(k2-1))
        
        # Bisection
        vf_low, vf_high = 0.0, 1.0
        for _ in range(50):
            vf_mid = (vf_low + vf_high) / 2
            if rr(vf_mid) > 0:
                vf_low = vf_mid
            else:
                vf_high = vf_mid
        
        vf = (vf_low + vf_high) / 2
        vf = max(0, min(1, vf))
        
        x1 = z1 / (1 + vf*(k1-1))
        y1 = k1 * x1
        
        return {
            "V_F": ureg.Quantity(vf, ""),
            "y1": ureg.Quantity(y1, ""),
            "x1": ureg.Quantity(x1, "")
        }


class CompressorWork(BaseEquation):
    """Calculate ideal compressor work for gas compression."""
    
    equation_id = "compressor_work"
    name = "Compressor Work (Isentropic)"
    category = "Thermodynamics"
    description = "Calculate ideal isentropic work for gas compression"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Suction pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Discharge pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₂"),
            EquationParameter("T1", "Suction temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="T₁"),
            EquationParameter("k", "Ratio of specific heats (Cp/Cv) - heat capacity ratio", 
                              "dimensionless", "", ParameterType.INPUT, symbol="k",
                              typical_range=(1.1, 1.7),
                              tooltip="Air=1.4, Methane=1.31, CO2=1.29"),
            EquationParameter("MW", "Molecular weight of gas", "dimensionless", "", 
                              ParameterType.INPUT, symbol="MW"),
            EquationParameter("m", "Mass flow rate", "mass_flow", "lb/hr", 
                              ParameterType.INPUT, symbol="ṁ"),
            EquationParameter("eta", "Isentropic efficiency", "dimensionless", "", 
                              ParameterType.INPUT, symbol="η", typical_range=(0.7, 0.9)),
            EquationParameter("W_ideal", "Ideal isentropic work", "power", "hp", 
                              ParameterType.OUTPUT, symbol="Wᵢₛ"),
            EquationParameter("W_actual", "Actual shaft work (accounting for efficiency)", 
                              "power", "hp", ParameterType.OUTPUT, symbol="Wₐₓ"),
            EquationParameter("T2_ideal", "Ideal discharge temperature", "temperature", "degF", 
                              ParameterType.OUTPUT, symbol="T₂ᵢₛ"),
        ]
    
    def _calculate(self, P1: pint.Quantity, P2: pint.Quantity, T1: pint.Quantity,
                   k: pint.Quantity, MW: pint.Quantity, m: pint.Quantity,
                   eta: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('psia').magnitude
        p2 = P2.to('psia').magnitude
        t1 = T1.to('degR').magnitude
        k_val = k.magnitude if hasattr(k, 'magnitude') else float(k)
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        m_val = m.to('lb/hr').magnitude
        eff = eta.magnitude if hasattr(eta, 'magnitude') else float(eta)
        
        R = 1545 / mw  # ft·lbf/(lb·°R)
        
        # Isentropic work: W = k/(k-1) * R * T1 * [(P2/P1)^((k-1)/k) - 1]
        ratio = (p2/p1)**((k_val-1)/k_val)
        w_spec = k_val/(k_val-1) * R * t1 * (ratio - 1)  # ft·lbf/lb
        
        w_ideal = m_val * w_spec / (3600 * 550)  # hp
        w_actual = w_ideal / eff
        
        # Discharge temperature
        t2_ideal = t1 * ratio - 459.67  # Convert back to °F
        
        return {
            "W_ideal": ureg.Quantity(w_ideal, "hp"),
            "W_actual": ureg.Quantity(w_actual, "hp"),
            "T2_ideal": ureg.Quantity(t2_ideal, "degF")
        }


THERMODYNAMICS_EQUATIONS = {
    'ideal_gas': IdealGas,
    'antoine': AntoineVaporPressure,
    'raoults_law': RaoultsLaw,
    'clausius_clapeyron': ClausiusClapeyron,
    'flash': FlashCalculation,
    'compressor_work': CompressorWork,
}


class VanDerWaals(BaseEquation):
    """Van der Waals equation of state for real gases."""
    
    equation_id = "van_der_waals"
    name = "Van der Waals Equation"
    category = "Thermodynamics"
    description = "(P + a/V²)(V - b) = RT - real gas behavior"
    reference = "Van der Waals, 1873"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("T", "Temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("P", "Pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P"),
            EquationParameter("Tc", "Critical temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="Tᶜ"),
            EquationParameter("Pc", "Critical pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="Pᶜ"),
            EquationParameter("Z", "Compressibility factor", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Z"),
            EquationParameter("V_molar", "Molar volume", "L/mol", "",
                              ParameterType.OUTPUT, symbol="Vₘ"),
        ]
    
    def _calculate(self, T: pint.Quantity, P: pint.Quantity, Tc: pint.Quantity, 
                   Pc: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        t = T.to('K').magnitude
        p = P.to('Pa').magnitude
        tc = Tc.to('K').magnitude
        pc = Pc.to('Pa').magnitude
        R = 8.314
        
        # Calculate a and b from critical properties
        a = 27 * R**2 * tc**2 / (64 * pc)
        b = R * tc / (8 * pc)
        
        # Solve cubic for V (iterative)
        v = R * t / p  # ideal gas initial guess
        for _ in range(20):
            v_new = (R * t / (p + a/v**2)) + b
            if abs(v_new - v) < 1e-10:
                break
            v = v_new
        
        z = p * v / (R * t)
        v_liters = v * 1000  # L/mol
        
        return {
            "Z": ureg.Quantity(z, ""),
            "V_molar": ureg.Quantity(v_liters, "L/mol")
        }


class ActivityCoefficient(BaseEquation):
    """Margules equation for activity coefficients."""
    
    equation_id = "activity_coefficient"
    name = "Margules Activity Coefficient"
    category = "Thermodynamics"
    description = "Calculate activity coefficients for non-ideal liquid mixtures"
    reference = "Smith, Van Ness & Abbott"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x1", "Mole fraction of component 1", "dimensionless", "",
                              ParameterType.INPUT, symbol="x₁", typical_range=(0, 1)),
            EquationParameter("A12", "Margules parameter A12", "dimensionless", "",
                              ParameterType.INPUT, symbol="A₁₂", typical_range=(0, 3)),
            EquationParameter("A21", "Margules parameter A21", "dimensionless", "",
                              ParameterType.INPUT, symbol="A₂₁", typical_range=(0, 3)),
            EquationParameter("gamma1", "Activity coefficient of component 1", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="γ₁"),
            EquationParameter("gamma2", "Activity coefficient of component 2", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="γ₂"),
        ]
    
    def _calculate(self, x1: pint.Quantity, A12: pint.Quantity, A21: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        x = x1.magnitude if hasattr(x1, 'magnitude') else float(x1)
        a12 = A12.magnitude if hasattr(A12, 'magnitude') else float(A12)
        a21 = A21.magnitude if hasattr(A21, 'magnitude') else float(A21)
        
        x2 = 1 - x
        ln_gamma1 = x2**2 * (a12 + 2*(a21 - a12)*x)
        ln_gamma2 = x**2 * (a21 + 2*(a12 - a21)*x2)
        
        return {
            "gamma1": ureg.Quantity(np.exp(ln_gamma1), ""),
            "gamma2": ureg.Quantity(np.exp(ln_gamma2), "")
        }


class JouleThomson(BaseEquation):
    """Joule-Thomson coefficient for throttling calculations."""
    
    equation_id = "joule_thomson"
    name = "Joule-Thomson Effect"
    category = "Thermodynamics"
    description = "Temperature change during isenthalpic throttling"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("mu_JT", "Joule-Thomson coefficient", "degC/bar", "degF/psi",
                              ParameterType.INPUT, symbol="μⱼₜ", typical_range=(-1, 1)),
            EquationParameter("P1", "Inlet pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Outlet pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P₂"),
            EquationParameter("T1", "Inlet temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="T₁"),
            EquationParameter("delta_T", "Temperature change", "delta_degC", "delta_degF",
                              ParameterType.OUTPUT, symbol="ΔT"),
            EquationParameter("T2", "Outlet temperature", "temperature", "degC",
                              ParameterType.OUTPUT, symbol="T₂"),
        ]
    
    def _calculate(self, mu_JT: pint.Quantity, P1: pint.Quantity, P2: pint.Quantity,
                   T1: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        mu = mu_JT.magnitude if hasattr(mu_JT, 'magnitude') else float(mu_JT)
        p1 = P1.to('bar').magnitude
        p2 = P2.to('bar').magnitude
        t1 = T1.to('degC').magnitude
        
        delta_t = mu * (p2 - p1)  # Usually negative (cooling)
        t2 = t1 + delta_t
        
        return {
            "delta_T": ureg.Quantity(delta_t, "delta_degC"),
            "T2": ureg.Quantity(t2, "degC")
        }


class HeatCapacityMixture(BaseEquation):
    """Heat capacity of ideal gas mixture."""
    
    equation_id = "heat_capacity_mixture"
    name = "Mixture Heat Capacity"
    category = "Thermodynamics"
    description = "Calculate Cp of ideal gas mixture from pure component values"
    reference = "Standard thermodynamics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("y1", "Mole fraction of component 1", "dimensionless", "",
                              ParameterType.INPUT, symbol="y₁", typical_range=(0, 1)),
            EquationParameter("Cp1", "Heat capacity of component 1", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.INPUT, symbol="Cₚ₁"),
            EquationParameter("Cp2", "Heat capacity of component 2", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.INPUT, symbol="Cₚ₂"),
            EquationParameter("Cp_mix", "Mixture heat capacity", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.OUTPUT, symbol="Cₚₘᵢₓ"),
        ]
    
    def _calculate(self, y1: pint.Quantity, Cp1: pint.Quantity, Cp2: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        y = y1.magnitude if hasattr(y1, 'magnitude') else float(y1)
        cp1 = Cp1.to('J/(mol*K)').magnitude
        cp2 = Cp2.to('J/(mol*K)').magnitude
        
        cp_mix = y * cp1 + (1-y) * cp2
        
        return {"Cp_mix": ureg.Quantity(cp_mix, "J/(mol*K)")}


# Update registry
THERMODYNAMICS_EQUATIONS = {
    'ideal_gas': IdealGas,
    'antoine': AntoineVaporPressure,
    'raoults_law': RaoultsLaw,
    'clausius_clapeyron': ClausiusClapeyron,
    'flash': FlashCalculation,
    'compressor_work': CompressorWork,
    'van_der_waals': VanDerWaals,
    'activity_coefficient': ActivityCoefficient,
    'joule_thomson': JouleThomson,
    'heat_capacity_mixture': HeatCapacityMixture,
}
