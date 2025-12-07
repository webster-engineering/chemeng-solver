"""
Pumps and Compressors equations - affinity laws, sizing, performance.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class AffinityLaws(BaseEquation):
    """Pump affinity laws for speed or impeller diameter changes."""
    
    equation_id = "affinity_laws"
    name = "Pump Affinity Laws"
    category = "Pumps"
    description = "Calculate new Q, H, P when changing pump speed or impeller diameter"
    reference = "Hydraulic Institute Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q1", "Original flow rate", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="Q₁"),
            EquationParameter("H1", "Original head", "length", "ft",
                              ParameterType.INPUT, symbol="H₁"),
            EquationParameter("P1", "Original power", "power", "hp",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("N1", "Original speed", "dimensionless", "rpm",
                              ParameterType.INPUT, symbol="N₁"),
            EquationParameter("N2", "New speed", "dimensionless", "rpm",
                              ParameterType.INPUT, symbol="N₂"),
            EquationParameter("Q2", "New flow rate", "volumetric_flow", "gpm",
                              ParameterType.OUTPUT, symbol="Q₂"),
            EquationParameter("H2", "New head", "length", "ft",
                              ParameterType.OUTPUT, symbol="H₂"),
            EquationParameter("P2", "New power", "power", "hp",
                              ParameterType.OUTPUT, symbol="P₂"),
        ]
    
    def _calculate(self, Q1: pint.Quantity, H1: pint.Quantity, P1: pint.Quantity,
                   N1: pint.Quantity, N2: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q1 = Q1.to('gpm').magnitude
        h1 = H1.to('ft').magnitude
        p1 = P1.to('hp').magnitude
        n1 = N1.magnitude if hasattr(N1, 'magnitude') else float(N1)
        n2 = N2.magnitude if hasattr(N2, 'magnitude') else float(N2)
        
        ratio = n2 / n1
        
        q2 = q1 * ratio
        h2 = h1 * ratio**2
        p2 = p1 * ratio**3
        
        return {
            "Q2": ureg.Quantity(q2, "gpm"),
            "H2": ureg.Quantity(h2, "ft"),
            "P2": ureg.Quantity(p2, "hp")
        }


class SpecificSpeed(BaseEquation):
    """Calculate pump specific speed for pump selection."""
    
    equation_id = "specific_speed"
    name = "Pump Specific Speed"
    category = "Pumps"
    description = "Calculate Ns for pump type selection"
    reference = "Hydraulic Institute"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("N", "Pump speed", "dimensionless", "rpm",
                              ParameterType.INPUT, symbol="N"),
            EquationParameter("Q", "Flow rate at BEP", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("H", "Head per stage at BEP", "length", "ft",
                              ParameterType.INPUT, symbol="H"),
            EquationParameter("Ns", "Specific speed", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Ns",
                              tooltip="<1000:radial, 1000-4000:mixed, >4000:axial"),
            EquationParameter("Nss", "Suction specific speed", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Nss",
                              tooltip="<8500 typical, >11000 may cavitate"),
        ]
    
    def _calculate(self, N: pint.Quantity, Q: pint.Quantity, H: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        n = N.magnitude if hasattr(N, 'magnitude') else float(N)
        q = Q.to('gpm').magnitude
        h = H.to('ft').magnitude
        
        ns = n * np.sqrt(q) / h**0.75
        nss = n * np.sqrt(q) / h**0.75  # Same formula, evaluate with NPSH for suction
        
        return {
            "Ns": ureg.Quantity(ns, ""),
            "Nss": ureg.Quantity(nss, "")
        }


class CompressorStaging(BaseEquation):
    """Calculate multi-stage compressor parameters."""
    
    equation_id = "compressor_staging"
    name = "Compressor Staging"
    category = "Pumps"
    description = "Calculate compression ratio and power for multi-stage compression"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Suction pressure", "pressure", "psia",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Final discharge pressure", "pressure", "psia",
                              ParameterType.INPUT, symbol="P₂"),
            EquationParameter("n_stages", "Number of stages", "dimensionless", "",
                              ParameterType.INPUT, symbol="n",
                              tooltip="Use ratio per stage < 4 typically"),
            EquationParameter("ratio_per_stage", "Compression ratio per stage", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="r"),
            EquationParameter("P_inter", "Intermediate pressures", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Pᵢₙₜ"),
            EquationParameter("total_ratio", "Overall compression ratio", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="r_total"),
        ]
    
    def _calculate(self, P1: pint.Quantity, P2: pint.Quantity, n_stages: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('psia').magnitude
        p2 = P2.to('psia').magnitude
        n = int(n_stages.magnitude if hasattr(n_stages, 'magnitude') else n_stages)
        
        total_ratio = p2 / p1
        ratio_per_stage = total_ratio**(1/n)
        
        # First intermediate pressure
        p_inter = p1 * ratio_per_stage
        
        return {
            "ratio_per_stage": ureg.Quantity(ratio_per_stage, ""),
            "P_inter": ureg.Quantity(p_inter, ""),
            "total_ratio": ureg.Quantity(total_ratio, "")
        }


class SurgeMargin(BaseEquation):
    """Calculate compressor surge margin."""
    
    equation_id = "surge_margin"
    name = "Compressor Surge Margin"
    category = "Pumps"
    description = "Calculate distance from surge line"
    reference = "API 617"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q_op", "Operating flow", "volumetric_flow", "acfm",
                              ParameterType.INPUT, symbol="Qop"),
            EquationParameter("Q_surge", "Surge flow at operating head", "volumetric_flow", "acfm",
                              ParameterType.INPUT, symbol="Qsurge"),
            EquationParameter("surge_margin", "Surge margin (%)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="SM%",
                              tooltip="Minimum 10% recommended"),
            EquationParameter("safe", "Is margin acceptable? (1=yes)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Safe"),
        ]
    
    def _calculate(self, Q_op: pint.Quantity, Q_surge: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        q_op = Q_op.magnitude if hasattr(Q_op, 'magnitude') else float(Q_op)
        q_surge = Q_surge.magnitude if hasattr(Q_surge, 'magnitude') else float(Q_surge)
        
        margin = (q_op - q_surge) / q_surge * 100
        safe = 1 if margin >= 10 else 0
        
        return {
            "surge_margin": ureg.Quantity(margin, ""),
            "safe": ureg.Quantity(safe, "")
        }


PUMPS_EQUATIONS = {
    'affinity_laws': AffinityLaws,
    'specific_speed': SpecificSpeed,
    'compressor_staging': CompressorStaging,
    'surge_margin': SurgeMargin,
}
