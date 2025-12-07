"""
Piping equations - pipe sizing, velocity, pressure drop, fittings.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class PipeVelocity(BaseEquation):
    """Calculate fluid velocity in a pipe from flow rate."""
    
    equation_id = "pipe_velocity"
    name = "Pipe Velocity"
    category = "Piping"
    description = "Calculate fluid velocity from volumetric flow and pipe diameter"
    reference = "Standard fluid mechanics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Volumetric flow rate", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("D", "Pipe inside diameter", "length", "in",
                              ParameterType.INPUT, symbol="D",
                              tooltip="Use schedule tables for ID"),
            EquationParameter("V", "Fluid velocity", "velocity", "ft/s",
                              ParameterType.OUTPUT, symbol="V",
                              tooltip="Liquids: 3-10 ft/s, Gases: 30-100 ft/s typical"),
        ]
    
    def _calculate(self, Q: pint.Quantity, D: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('ft**3/s').magnitude
        d = D.to('ft').magnitude
        
        area = np.pi * (d/2)**2
        v = q / area
        
        return {"V": ureg.Quantity(v, "ft/s")}


class PipeScheduleLookup(BaseEquation):
    """Look up pipe dimensions from nominal size and schedule."""
    
    equation_id = "pipe_schedule"
    name = "Pipe Schedule Lookup"
    category = "Piping"
    description = "Get pipe ID, wall thickness, and weight from NPS and schedule"
    reference = "ASME B36.10M"
    
    # Common pipe data: NPS -> {schedule: (OD, wall_thickness)}
    PIPE_DATA = {
        0.5: {40: (0.840, 0.109), 80: (0.840, 0.147)},
        0.75: {40: (1.050, 0.113), 80: (1.050, 0.154)},
        1: {40: (1.315, 0.133), 80: (1.315, 0.179)},
        1.5: {40: (1.900, 0.145), 80: (1.900, 0.200)},
        2: {40: (2.375, 0.154), 80: (2.375, 0.218)},
        3: {40: (3.500, 0.216), 80: (3.500, 0.300)},
        4: {40: (4.500, 0.237), 80: (4.500, 0.337)},
        6: {40: (6.625, 0.280), 80: (6.625, 0.432)},
        8: {40: (8.625, 0.322), 80: (8.625, 0.500)},
        10: {40: (10.750, 0.365), 80: (10.750, 0.500)},
        12: {40: (12.750, 0.406), 80: (12.750, 0.500)},
    }
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("NPS", "Nominal pipe size", "dimensionless", "",
                              ParameterType.INPUT, symbol="NPS",
                              tooltip="Common: 0.5, 1, 2, 3, 4, 6, 8, 10, 12"),
            EquationParameter("schedule", "Pipe schedule", "dimensionless", "",
                              ParameterType.INPUT, symbol="Sch",
                              tooltip="40 (Std), 80 (XS), etc."),
            EquationParameter("OD", "Outside diameter", "length", "in",
                              ParameterType.OUTPUT, symbol="OD"),
            EquationParameter("wall", "Wall thickness", "length", "in",
                              ParameterType.OUTPUT, symbol="t"),
            EquationParameter("ID", "Inside diameter", "length", "in",
                              ParameterType.OUTPUT, symbol="ID"),
            EquationParameter("area", "Flow area", "area", "in**2",
                              ParameterType.OUTPUT, symbol="A"),
        ]
    
    def _calculate(self, NPS: pint.Quantity, schedule: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        nps = NPS.magnitude if hasattr(NPS, 'magnitude') else float(NPS)
        sch = int(schedule.magnitude if hasattr(schedule, 'magnitude') else schedule)
        
        if nps in self.PIPE_DATA and sch in self.PIPE_DATA[nps]:
            od, wall = self.PIPE_DATA[nps][sch]
        else:
            od, wall = 1.0, 0.1  # Default
        
        id_val = od - 2 * wall
        area = np.pi * (id_val/2)**2
        
        return {
            "OD": ureg.Quantity(od, "in"),
            "wall": ureg.Quantity(wall, "in"),
            "ID": ureg.Quantity(id_val, "in"),
            "area": ureg.Quantity(area, "in**2")
        }


class EquivalentLength(BaseEquation):
    """Calculate equivalent length of pipe fittings."""
    
    equation_id = "equivalent_length"
    name = "Fitting Equivalent Length"
    category = "Piping"
    description = "Calculate equivalent pipe length for fittings using L/D method"
    reference = "Crane Technical Paper 410"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("D", "Pipe inside diameter", "length", "in",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("n_elbow90", "Number of 90° elbows", "dimensionless", "",
                              ParameterType.INPUT, symbol="n₉₀", tooltip="L/D = 30"),
            EquationParameter("n_elbow45", "Number of 45° elbows", "dimensionless", "",
                              ParameterType.INPUT, symbol="n₄₅", tooltip="L/D = 16"),
            EquationParameter("n_tee_thru", "Number of tees (through)", "dimensionless", "",
                              ParameterType.INPUT, symbol="nTₜ", tooltip="L/D = 20"),
            EquationParameter("n_tee_branch", "Number of tees (branch)", "dimensionless", "",
                              ParameterType.INPUT, symbol="nTᵦ", tooltip="L/D = 60"),
            EquationParameter("n_gate", "Number of gate valves", "dimensionless", "",
                              ParameterType.INPUT, symbol="nGV", tooltip="L/D = 8"),
            EquationParameter("n_globe", "Number of globe valves", "dimensionless", "",
                              ParameterType.INPUT, symbol="nGlb", tooltip="L/D = 340"),
            EquationParameter("n_check", "Number of check valves", "dimensionless", "",
                              ParameterType.INPUT, symbol="nChk", tooltip="L/D = 100"),
            EquationParameter("L_eq", "Total equivalent length", "length", "ft",
                              ParameterType.OUTPUT, symbol="Lₑq"),
        ]
    
    def _calculate(self, D: pint.Quantity, n_elbow90: pint.Quantity, n_elbow45: pint.Quantity,
                   n_tee_thru: pint.Quantity, n_tee_branch: pint.Quantity,
                   n_gate: pint.Quantity, n_globe: pint.Quantity, n_check: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        d = D.to('ft').magnitude
        
        def get_val(q):
            return q.magnitude if hasattr(q, 'magnitude') else float(q)
        
        # L/D values from Crane 410
        l_d = {
            'elbow90': 30 * get_val(n_elbow90),
            'elbow45': 16 * get_val(n_elbow45),
            'tee_thru': 20 * get_val(n_tee_thru),
            'tee_branch': 60 * get_val(n_tee_branch),
            'gate': 8 * get_val(n_gate),
            'globe': 340 * get_val(n_globe),
            'check': 100 * get_val(n_check),
        }
        
        total_l_d = sum(l_d.values())
        l_eq = total_l_d * d
        
        return {"L_eq": ureg.Quantity(l_eq, "ft")}


class ThermalExpansion(BaseEquation):
    """Calculate pipe thermal expansion."""
    
    equation_id = "thermal_expansion"
    name = "Pipe Thermal Expansion"
    category = "Piping"
    description = "Calculate length change due to temperature change"
    reference = "Standard thermal expansion"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("L", "Pipe length", "length", "ft",
                              ParameterType.INPUT, symbol="L"),
            EquationParameter("T1", "Initial temperature", "temperature", "degF",
                              ParameterType.INPUT, symbol="T₁"),
            EquationParameter("T2", "Final temperature", "temperature", "degF",
                              ParameterType.INPUT, symbol="T₂"),
            EquationParameter("alpha", "Thermal expansion coefficient", "dimensionless", "",
                              ParameterType.INPUT, symbol="α",
                              tooltip="Carbon steel: 6.5e-6 /°F, SS: 9.6e-6 /°F"),
            EquationParameter("dL", "Length change", "length", "in",
                              ParameterType.OUTPUT, symbol="ΔL"),
            EquationParameter("expansion_loop", "Suggested expansion loop leg", "length", "ft",
                              ParameterType.OUTPUT, symbol="Loop",
                              tooltip="Based on 3D expansion loop"),
        ]
    
    def _calculate(self, L: pint.Quantity, T1: pint.Quantity, T2: pint.Quantity,
                   alpha: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        l = L.to('in').magnitude
        t1 = T1.to('degF').magnitude
        t2 = T2.to('degF').magnitude
        a = alpha.magnitude if hasattr(alpha, 'magnitude') else float(alpha)
        
        dt = t2 - t1
        dl = a * l * dt
        
        # Expansion loop sizing (empirical): L = sqrt(3 * D * dL / allowable_stress_factor)
        # Simplified: leg ≈ sqrt(6 * D * dL) assuming 4" pipe
        loop_leg = np.sqrt(6 * 4 * abs(dl)) if dl != 0 else 0
        
        return {
            "dL": ureg.Quantity(dl, "in"),
            "expansion_loop": ureg.Quantity(loop_leg / 12, "ft")
        }


class PipeSupportSpacing(BaseEquation):
    """Calculate recommended pipe support spacing."""
    
    equation_id = "support_spacing"
    name = "Pipe Support Spacing"
    category = "Piping"
    description = "Calculate maximum span between pipe supports"
    reference = "MSS SP-69"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("NPS", "Nominal pipe size", "dimensionless", "",
                              ParameterType.INPUT, symbol="NPS"),
            EquationParameter("material", "Pipe material (1=CS, 2=SS)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Mat",
                              tooltip="1=Carbon Steel, 2=Stainless Steel"),
            EquationParameter("insulated", "Is pipe insulated? (0=no, 1=yes)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Ins"),
            EquationParameter("span_water", "Max span (water filled)", "length", "ft",
                              ParameterType.OUTPUT, symbol="Span_w"),
            EquationParameter("span_gas", "Max span (gas/empty)", "length", "ft",
                              ParameterType.OUTPUT, symbol="Span_g"),
        ]
    
    # Approximate spans from MSS SP-69 (ft)
    SPANS = {
        1: (7, 9), 2: (10, 13), 3: (12, 15), 4: (14, 17),
        6: (17, 21), 8: (19, 24), 10: (22, 26), 12: (23, 30)
    }
    
    def _calculate(self, NPS: pint.Quantity, material: pint.Quantity, insulated: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        nps = int(NPS.magnitude if hasattr(NPS, 'magnitude') else NPS)
        mat = material.magnitude if hasattr(material, 'magnitude') else 1
        ins = insulated.magnitude if hasattr(insulated, 'magnitude') else 0
        
        if nps in self.SPANS:
            span_w, span_g = self.SPANS[nps]
        else:
            span_w, span_g = 10, 15
        
        # Reduce for insulation
        if ins > 0:
            span_w *= 0.9
            span_g *= 0.9
        
        return {
            "span_water": ureg.Quantity(span_w, "ft"),
            "span_gas": ureg.Quantity(span_g, "ft")
        }


PIPING_EQUATIONS = {
    'pipe_velocity': PipeVelocity,
    'pipe_schedule': PipeScheduleLookup,
    'equivalent_length': EquivalentLength,
    'thermal_expansion': ThermalExpansion,
    'support_spacing': PipeSupportSpacing,
}
