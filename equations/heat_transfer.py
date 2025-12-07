"""
Heat Transfer equations - LMTD, NTU-effectiveness, heat transfer coefficients.
Each equation includes detailed variable descriptions and typical units.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class LMTD(BaseEquation):
    """Log Mean Temperature Difference calculation."""
    
    equation_id = "lmtd"
    name = "Log Mean Temperature Difference"
    category = "Heat Transfer"
    description = "Calculate LMTD for shell-and-tube or double-pipe heat exchanger design"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("T_hot_in", "Hot fluid inlet temperature - entering hot side", 
                              "temperature", "degF", ParameterType.INPUT, symbol="Tₕᵢ"),
            EquationParameter("T_hot_out", "Hot fluid outlet temperature - leaving hot side", 
                              "temperature", "degF", ParameterType.INPUT, symbol="Tₕₒ"),
            EquationParameter("T_cold_in", "Cold fluid inlet temperature - entering cold side", 
                              "temperature", "degF", ParameterType.INPUT, symbol="Tcᵢ"),
            EquationParameter("T_cold_out", "Cold fluid outlet temperature - leaving cold side", 
                              "temperature", "degF", ParameterType.INPUT, symbol="Tcₒ"),
            EquationParameter("LMTD", "Log mean temperature difference - driving force for heat transfer", 
                              "temperature_difference", "delta_degF", ParameterType.OUTPUT, symbol="LMTD",
                              tooltip="Use in Q = U·A·LMTD"),
        ]
    
    def _calculate(self, T_hot_in: pint.Quantity, T_hot_out: pint.Quantity,
                   T_cold_in: pint.Quantity, T_cold_out: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        thi = T_hot_in.to('degF').magnitude
        tho = T_hot_out.to('degF').magnitude
        tci = T_cold_in.to('degF').magnitude
        tco = T_cold_out.to('degF').magnitude
        
        dt1 = thi - tco
        dt2 = tho - tci
        
        if dt1 <= 0 or dt2 <= 0:
            raise ValueError("Invalid temperature approach - check temperature values")
        
        if abs(dt1 - dt2) < 0.001:
            lmtd = dt1
        else:
            lmtd = (dt1 - dt2) / np.log(dt1 / dt2)
        
        return {"LMTD": ureg.Quantity(lmtd, "delta_degF")}


class HeatExchangerArea(BaseEquation):
    """Calculate required heat exchanger area."""
    
    equation_id = "hx_area"
    name = "Heat Exchanger Area"
    category = "Heat Transfer"
    description = "Calculate required heat transfer area using Q = U·A·ΔTlm"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Heat duty - total heat to be transferred", 
                              "power", "BTU/hr", ParameterType.INPUT, symbol="Q",
                              tooltip="From energy balance: Q = m·Cp·ΔT"),
            EquationParameter("U", "Overall heat transfer coefficient", 
                              "heat_transfer_coefficient", "BTU/(hr*ft**2*delta_degF)", 
                              ParameterType.INPUT, symbol="U",
                              tooltip="Typical: 50-150 for liquid-liquid, 5-50 for gas-liquid"),
            EquationParameter("LMTD", "Log mean temperature difference", 
                              "temperature_difference", "delta_degF", ParameterType.INPUT, symbol="ΔTₗₘ"),
            EquationParameter("F", "LMTD correction factor - for multi-pass exchangers", 
                              "dimensionless", "", ParameterType.INPUT, symbol="F",
                              typical_range=(0.75, 1.0), tooltip="1.0 for true counterflow, <1 for multi-pass"),
            EquationParameter("A", "Required heat transfer area", "area", "ft**2", 
                              ParameterType.OUTPUT, symbol="A"),
        ]
    
    def _calculate(self, Q: pint.Quantity, U: pint.Quantity, LMTD: pint.Quantity,
                   F: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('BTU/hr').magnitude
        u = U.to('BTU/(hr*ft**2*delta_degF)').magnitude
        lmtd = LMTD.to('delta_degF').magnitude
        f = F.magnitude if hasattr(F, 'magnitude') else float(F)
        
        area = q / (u * f * lmtd)
        
        return {"A": ureg.Quantity(area, "ft**2")}


class HeatDuty(BaseEquation):
    """Calculate heat duty from flow and temperature change."""
    
    equation_id = "heat_duty"
    name = "Heat Duty Calculation"
    category = "Heat Transfer"
    description = "Calculate heat transfer rate Q = m·Cp·ΔT"
    reference = "Standard thermodynamics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("m", "Mass flow rate of fluid", "mass_flow", "lb/hr", 
                              ParameterType.INPUT, symbol="ṁ"),
            EquationParameter("Cp", "Specific heat capacity - heat needed to raise 1 lb by 1°F", 
                              "specific_heat", "BTU/(lb*delta_degF)", ParameterType.INPUT, symbol="Cₚ",
                              tooltip="Water=1.0, Oil≈0.5, Air≈0.24 BTU/(lb·°F)"),
            EquationParameter("T_in", "Inlet temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="Tᵢₙ"),
            EquationParameter("T_out", "Outlet temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="Tₒᵤₜ"),
            EquationParameter("Q", "Heat duty - rate of heat transfer", "power", "BTU/hr", 
                              ParameterType.OUTPUT, symbol="Q"),
        ]
    
    def _calculate(self, m: pint.Quantity, Cp: pint.Quantity, T_in: pint.Quantity,
                   T_out: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        mass_flow = m.to('lb/hr').magnitude
        cp = Cp.to('BTU/(lb*delta_degF)').magnitude
        dt = abs(T_out.to('degF').magnitude - T_in.to('degF').magnitude)
        
        q = mass_flow * cp * dt
        
        return {"Q": ureg.Quantity(q, "BTU/hr")}


class NTUEffectiveness(BaseEquation):
    """NTU-Effectiveness method for heat exchangers."""
    
    equation_id = "ntu_eff"
    name = "NTU-Effectiveness Method"
    category = "Heat Transfer"
    description = "Calculate heat exchanger effectiveness from NTU (counterflow)"
    reference = "Kays & London, Compact Heat Exchangers"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("NTU", "Number of transfer units - dimensionless size of HX", 
                              "dimensionless", "", ParameterType.INPUT, symbol="NTU", 
                              typical_range=(0.1, 10),
                              tooltip="NTU = U·A / Cmin"),
            EquationParameter("Cr", "Heat capacity ratio - ratio of Cmin to Cmax", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Cᵣ", 
                              typical_range=(0, 1),
                              tooltip="Cr = Cmin/Cmax, where C = m·Cp"),
            EquationParameter("effectiveness", "Heat exchanger effectiveness - actual Q / max Q", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="ε", 
                              typical_range=(0, 1),
                              tooltip="ε = Q_actual / Q_max"),
        ]
    
    def _calculate(self, NTU: pint.Quantity, Cr: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        ntu = NTU.magnitude if hasattr(NTU, 'magnitude') else float(NTU)
        cr = Cr.magnitude if hasattr(Cr, 'magnitude') else float(Cr)
        
        if cr == 0:  # One fluid at constant temperature (condenser/evaporator)
            eff = 1 - np.exp(-ntu)
        elif cr == 1:
            eff = ntu / (1 + ntu)
        else:  # Counterflow
            eff = (1 - np.exp(-ntu*(1-cr))) / (1 - cr*np.exp(-ntu*(1-cr)))
        
        return {"effectiveness": ureg.Quantity(min(eff, 1.0), "")}


class OverallU(BaseEquation):
    """Calculate overall heat transfer coefficient."""
    
    equation_id = "overall_u"
    name = "Overall Heat Transfer Coefficient"
    category = "Heat Transfer"
    description = "Calculate U from individual thermal resistances"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("hi", "Inside film coefficient - convection on tube side", 
                              "heat_transfer_coefficient", "BTU/(hr*ft**2*delta_degF)", 
                              ParameterType.INPUT, symbol="hᵢ",
                              tooltip="From correlations (Dittus-Boelter, Sieder-Tate)"),
            EquationParameter("ho", "Outside film coefficient - convection on shell side", 
                              "heat_transfer_coefficient", "BTU/(hr*ft**2*delta_degF)", 
                              ParameterType.INPUT, symbol="hₒ"),
            EquationParameter("Rfi", "Inside fouling factor - deposit thermal resistance", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Rfᵢ", 
                              required=False,
                              tooltip="TEMA: 0.001-0.003 hr·ft²·°F/BTU typical"),
            EquationParameter("Rfo", "Outside fouling factor", "dimensionless", "", 
                              ParameterType.INPUT, symbol="Rfₒ", required=False),
            EquationParameter("U", "Overall coefficient - based on outside area", 
                              "heat_transfer_coefficient", "BTU/(hr*ft**2*delta_degF)", 
                              ParameterType.OUTPUT, symbol="U"),
        ]
    
    def _calculate(self, hi: pint.Quantity, ho: pint.Quantity,
                   Rfi: pint.Quantity = None, Rfo: pint.Quantity = None,
                   **kwargs) -> Dict[str, pint.Quantity]:
        hi_val = hi.to('BTU/(hr*ft**2*delta_degF)').magnitude
        ho_val = ho.to('BTU/(hr*ft**2*delta_degF)').magnitude
        
        rfi = Rfi.magnitude if Rfi else 0
        rfo = Rfo.magnitude if Rfo else 0
        
        r_total = 1/hi_val + 1/ho_val + rfi + rfo
        u = 1 / r_total
        
        return {"U": ureg.Quantity(u, "BTU/(hr*ft**2*delta_degF)")}


class ConvectionCoefficient(BaseEquation):
    """Estimate convection coefficient using Dittus-Boelter correlation."""
    
    equation_id = "convection_coeff"
    name = "Convection Coefficient (Dittus-Boelter)"
    category = "Heat Transfer"
    description = "Estimate film coefficient for turbulent flow in tubes"
    reference = "Dittus-Boelter correlation (Re > 10000)"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Re", "Reynolds number - must be > 10000 for this correlation", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Re",
                              typical_range=(10000, 1e7)),
            EquationParameter("Pr", "Prandtl number - Cp·μ/k, fluid property", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Pr",
                              typical_range=(0.5, 100),
                              tooltip="Water≈7, Air≈0.7, Oils≈50-1000"),
            EquationParameter("k", "Thermal conductivity of fluid", "thermal_conductivity", 
                              "BTU/(hr*ft*delta_degF)", ParameterType.INPUT, symbol="k"),
            EquationParameter("D", "Tube inside diameter", "length", "in", 
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("heating", "Is fluid being heated? (vs cooled)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Heat?",
                              tooltip="Enter 1 for heating, 0 for cooling"),
            EquationParameter("h", "Convection heat transfer coefficient", 
                              "heat_transfer_coefficient", "BTU/(hr*ft**2*delta_degF)", 
                              ParameterType.OUTPUT, symbol="h"),
        ]
    
    def _calculate(self, Re: pint.Quantity, Pr: pint.Quantity, k: pint.Quantity,
                   D: pint.Quantity, heating: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        re = Re.magnitude if hasattr(Re, 'magnitude') else float(Re)
        pr = Pr.magnitude if hasattr(Pr, 'magnitude') else float(Pr)
        k_val = k.to('BTU/(hr*ft*delta_degF)').magnitude
        d = D.to('ft').magnitude
        is_heating = heating.magnitude if hasattr(heating, 'magnitude') else float(heating)
        
        # Dittus-Boelter: Nu = 0.023 * Re^0.8 * Pr^n (n=0.4 heating, 0.3 cooling)
        n = 0.4 if is_heating > 0.5 else 0.3
        nu = 0.023 * re**0.8 * pr**n
        
        h = nu * k_val / d
        
        return {"h": ureg.Quantity(h, "BTU/(hr*ft**2*delta_degF)")}


HEAT_TRANSFER_EQUATIONS = {
    'lmtd': LMTD,
    'hx_area': HeatExchangerArea,
    'heat_duty': HeatDuty,
    'ntu_eff': NTUEffectiveness,
    'overall_u': OverallU,
    'convection_coeff': ConvectionCoefficient,
}
