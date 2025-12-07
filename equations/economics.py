"""
Economics equations - cost estimation, NPV, payback.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class EquipmentCostScaling(BaseEquation):
    """Scale equipment cost using six-tenths rule."""
    
    equation_id = "cost_scaling"
    name = "Equipment Cost Scaling"
    category = "Economics"
    description = "Scale equipment cost from known size using power law"
    reference = "Peters & Timmerhaus"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("C_known", "Known equipment cost", "dimensionless", "$",
                              ParameterType.INPUT, symbol="C₁"),
            EquationParameter("S_known", "Known equipment size/capacity", "dimensionless", "",
                              ParameterType.INPUT, symbol="S₁"),
            EquationParameter("S_new", "New equipment size/capacity", "dimensionless", "",
                              ParameterType.INPUT, symbol="S₂"),
            EquationParameter("n", "Scaling exponent", "dimensionless", "",
                              ParameterType.INPUT, symbol="n",
                              tooltip="0.6 typical (six-tenths rule), range 0.4-0.8"),
            EquationParameter("CEPCI_known", "Cost index at known cost", "dimensionless", "",
                              ParameterType.INPUT, symbol="CEPCI₁",
                              tooltip="E.g., 2019: 607, 2023: 800"),
            EquationParameter("CEPCI_new", "Current cost index", "dimensionless", "",
                              ParameterType.INPUT, symbol="CEPCI₂"),
            EquationParameter("C_new", "Estimated new equipment cost", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="C₂ ($)"),
        ]
    
    def _calculate(self, C_known: pint.Quantity, S_known: pint.Quantity, S_new: pint.Quantity,
                   n: pint.Quantity, CEPCI_known: pint.Quantity, CEPCI_new: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        c1 = C_known.magnitude if hasattr(C_known, 'magnitude') else float(C_known)
        s1 = S_known.magnitude if hasattr(S_known, 'magnitude') else float(S_known)
        s2 = S_new.magnitude if hasattr(S_new, 'magnitude') else float(S_new)
        exp = n.magnitude if hasattr(n, 'magnitude') else float(n)
        cepci1 = CEPCI_known.magnitude if hasattr(CEPCI_known, 'magnitude') else float(CEPCI_known)
        cepci2 = CEPCI_new.magnitude if hasattr(CEPCI_new, 'magnitude') else float(CEPCI_new)
        
        c2 = c1 * (s2/s1)**exp * (cepci2/cepci1)
        
        return {"C_new": ureg.Quantity(c2, "")}


class NPVCalculation(BaseEquation):
    """Calculate Net Present Value of a project."""
    
    equation_id = "npv"
    name = "Net Present Value (NPV)"
    category = "Economics"
    description = "Calculate NPV from annual cash flows"
    reference = "Standard financial analysis"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("initial_investment", "Initial investment (negative)", "dimensionless", "$",
                              ParameterType.INPUT, symbol="I₀"),
            EquationParameter("annual_cashflow", "Annual net cash flow", "dimensionless", "$",
                              ParameterType.INPUT, symbol="CF"),
            EquationParameter("years", "Project life", "dimensionless", "years",
                              ParameterType.INPUT, symbol="n"),
            EquationParameter("discount_rate", "Discount rate (decimal)", "dimensionless", "",
                              ParameterType.INPUT, symbol="r",
                              tooltip="E.g., 0.10 for 10%"),
            EquationParameter("salvage", "Salvage value at end", "dimensionless", "$",
                              ParameterType.INPUT, symbol="S"),
            EquationParameter("NPV", "Net Present Value", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="NPV ($)"),
            EquationParameter("PI", "Profitability Index", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="PI"),
        ]
    
    def _calculate(self, initial_investment: pint.Quantity, annual_cashflow: pint.Quantity,
                   years: pint.Quantity, discount_rate: pint.Quantity, salvage: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        i0 = initial_investment.magnitude if hasattr(initial_investment, 'magnitude') else float(initial_investment)
        cf = annual_cashflow.magnitude if hasattr(annual_cashflow, 'magnitude') else float(annual_cashflow)
        n = int(years.magnitude if hasattr(years, 'magnitude') else years)
        r = discount_rate.magnitude if hasattr(discount_rate, 'magnitude') else float(discount_rate)
        s = salvage.magnitude if hasattr(salvage, 'magnitude') else float(salvage)
        
        # NPV = -I0 + sum(CF/(1+r)^t) + S/(1+r)^n
        npv = -i0
        for t in range(1, n + 1):
            npv += cf / (1 + r)**t
        npv += s / (1 + r)**n
        
        # Profitability Index
        pi = (npv + i0) / i0 if i0 != 0 else 0
        
        return {
            "NPV": ureg.Quantity(npv, ""),
            "PI": ureg.Quantity(pi, "")
        }


class PaybackPeriod(BaseEquation):
    """Calculate simple and discounted payback period."""
    
    equation_id = "payback"
    name = "Payback Period"
    category = "Economics"
    description = "Calculate time to recover initial investment"
    reference = "Standard financial analysis"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("initial_investment", "Initial investment", "dimensionless", "$",
                              ParameterType.INPUT, symbol="I₀"),
            EquationParameter("annual_cashflow", "Annual net cash flow", "dimensionless", "$",
                              ParameterType.INPUT, symbol="CF"),
            EquationParameter("discount_rate", "Discount rate (for DPP)", "dimensionless", "",
                              ParameterType.INPUT, symbol="r"),
            EquationParameter("simple_payback", "Simple payback period", "time", "years",
                              ParameterType.OUTPUT, symbol="SPP"),
            EquationParameter("discounted_payback", "Discounted payback period", "time", "years",
                              ParameterType.OUTPUT, symbol="DPP"),
        ]
    
    def _calculate(self, initial_investment: pint.Quantity, annual_cashflow: pint.Quantity,
                   discount_rate: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        i0 = initial_investment.magnitude if hasattr(initial_investment, 'magnitude') else float(initial_investment)
        cf = annual_cashflow.magnitude if hasattr(annual_cashflow, 'magnitude') else float(annual_cashflow)
        r = discount_rate.magnitude if hasattr(discount_rate, 'magnitude') else float(discount_rate)
        
        # Simple payback
        spp = i0 / cf if cf > 0 else float('inf')
        
        # Discounted payback
        cumulative = 0
        year = 0
        while cumulative < i0 and year < 100:
            year += 1
            cumulative += cf / (1 + r)**year
        dpp = year if cumulative >= i0 else float('inf')
        
        return {
            "simple_payback": ureg.Quantity(spp, "year"),
            "discounted_payback": ureg.Quantity(dpp, "year")
        }


class OperatingCost(BaseEquation):
    """Calculate annual operating cost for a process unit."""
    
    equation_id = "operating_cost"
    name = "Operating Cost Estimate"
    category = "Economics"
    description = "Estimate annual operating costs from utilities and labor"
    reference = "Peters & Timmerhaus"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("steam_rate", "Steam usage", "mass_flow", "lb/hr",
                              ParameterType.INPUT, symbol="Steam"),
            EquationParameter("steam_cost", "Steam cost ($/1000 lb)", "dimensionless", "",
                              ParameterType.INPUT, symbol="$/klb",
                              tooltip="Typical: $5-15 per 1000 lb"),
            EquationParameter("elec_rate", "Electricity usage", "power", "kW",
                              ParameterType.INPUT, symbol="Power"),
            EquationParameter("elec_cost", "Electricity cost ($/kWh)", "dimensionless", "",
                              ParameterType.INPUT, symbol="$/kWh",
                              tooltip="Typical: $0.05-0.12"),
            EquationParameter("cw_rate", "Cooling water usage", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="CW"),
            EquationParameter("cw_cost", "CW cost ($/1000 gal)", "dimensionless", "",
                              ParameterType.INPUT, symbol="$/kgal",
                              tooltip="Typical: $0.10-0.50"),
            EquationParameter("operators", "Number of operators", "dimensionless", "",
                              ParameterType.INPUT, symbol="Ops"),
            EquationParameter("labor_rate", "Loaded labor rate ($/hr)", "dimensionless", "",
                              ParameterType.INPUT, symbol="$/hr"),
            EquationParameter("hours", "Operating hours per year", "dimensionless", "",
                              ParameterType.INPUT, symbol="hr/yr",
                              tooltip="8760 = 24/7, 8000 = typical"),
            EquationParameter("annual_cost", "Total annual operating cost", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="$/yr"),
        ]
    
    def _calculate(self, steam_rate: pint.Quantity, steam_cost: pint.Quantity,
                   elec_rate: pint.Quantity, elec_cost: pint.Quantity,
                   cw_rate: pint.Quantity, cw_cost: pint.Quantity,
                   operators: pint.Quantity, labor_rate: pint.Quantity,
                   hours: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        
        def get_val(q):
            return q.magnitude if hasattr(q, 'magnitude') else float(q)
        
        steam = steam_rate.to('lb/hr').magnitude
        steam_c = get_val(steam_cost)
        elec = elec_rate.to('kW').magnitude
        elec_c = get_val(elec_cost)
        cw = cw_rate.to('gpm').magnitude
        cw_c = get_val(cw_cost)
        ops = get_val(operators)
        labor = get_val(labor_rate)
        hr = get_val(hours)
        
        steam_annual = steam * hr * steam_c / 1000
        elec_annual = elec * hr * elec_c
        cw_annual = cw * 60 * hr * cw_c / 1000
        labor_annual = ops * labor * hr
        
        total = steam_annual + elec_annual + cw_annual + labor_annual
        
        return {"annual_cost": ureg.Quantity(total, "")}


ECONOMICS_EQUATIONS = {
    'cost_scaling': EquipmentCostScaling,
    'npv': NPVCalculation,
    'payback': PaybackPeriod,
    'operating_cost': OperatingCost,
}
