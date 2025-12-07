"""
Vessel and Tank equations - volume, design, agitation.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class TankVolumeVertical(BaseEquation):
    """Calculate vertical cylindrical tank volume."""
    
    equation_id = "tank_vertical"
    name = "Vertical Tank Volume"
    category = "Vessels"
    description = "Calculate volume and dimensions of vertical cylindrical tank"
    reference = "Standard geometry"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("D", "Tank diameter", "length", "ft",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("H", "Straight shell height", "length", "ft",
                              ParameterType.INPUT, symbol="H"),
            EquationParameter("head_type", "Head type (1=flat, 2=2:1 ellip, 3=hemi)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Head",
                              tooltip="1=Flat, 2=2:1 Elliptical, 3=Hemispherical"),
            EquationParameter("V_shell", "Shell volume", "volume", "gal",
                              ParameterType.OUTPUT, symbol="Vshell"),
            EquationParameter("V_heads", "Total head volume (both)", "volume", "gal",
                              ParameterType.OUTPUT, symbol="Vheads"),
            EquationParameter("V_total", "Total tank volume", "volume", "gal",
                              ParameterType.OUTPUT, symbol="Vtotal"),
        ]
    
    def _calculate(self, D: pint.Quantity, H: pint.Quantity, head_type: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        d = D.to('ft').magnitude
        h = H.to('ft').magnitude
        ht = int(head_type.magnitude if hasattr(head_type, 'magnitude') else head_type)
        
        r = d / 2
        v_shell = np.pi * r**2 * h  # ft³
        
        if ht == 3:  # Hemispherical
            v_one_head = (2/3) * np.pi * r**3
        elif ht == 2:  # 2:1 Elliptical
            v_one_head = (np.pi * d**3) / 24
        else:  # Flat
            v_one_head = 0
        
        v_heads = 2 * v_one_head
        v_total = v_shell + v_heads
        
        # Convert to gallons (1 ft³ = 7.48052 gal)
        return {
            "V_shell": ureg.Quantity(v_shell * 7.48052, "gal"),
            "V_heads": ureg.Quantity(v_heads * 7.48052, "gal"),
            "V_total": ureg.Quantity(v_total * 7.48052, "gal")
        }


class TankVolumeHorizontal(BaseEquation):
    """Calculate horizontal cylindrical tank partial volume."""
    
    equation_id = "tank_horizontal"
    name = "Horizontal Tank Volume"
    category = "Vessels"
    description = "Calculate volume at given liquid level in horizontal tank"
    reference = "Standard geometry"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("D", "Tank diameter", "length", "ft",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("L", "Tank length (straight)", "length", "ft",
                              ParameterType.INPUT, symbol="L"),
            EquationParameter("h", "Liquid level height", "length", "ft",
                              ParameterType.INPUT, symbol="h"),
            EquationParameter("V_partial", "Volume at level h", "volume", "gal",
                              ParameterType.OUTPUT, symbol="V"),
            EquationParameter("pct_full", "Percent full", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="%Full"),
            EquationParameter("V_total", "Total tank volume", "volume", "gal",
                              ParameterType.OUTPUT, symbol="Vtotal"),
        ]
    
    def _calculate(self, D: pint.Quantity, L: pint.Quantity, h: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        d = D.to('ft').magnitude
        length = L.to('ft').magnitude
        level = h.to('ft').magnitude
        
        r = d / 2
        
        # Partial cylinder volume formula
        if level >= d:
            a_partial = np.pi * r**2
        elif level <= 0:
            a_partial = 0
        else:
            theta = 2 * np.arccos((r - level) / r)
            a_partial = r**2 * (theta - np.sin(theta)) / 2
        
        v_partial = a_partial * length
        v_total = np.pi * r**2 * length
        pct = (v_partial / v_total) * 100 if v_total > 0 else 0
        
        return {
            "V_partial": ureg.Quantity(v_partial * 7.48052, "gal"),
            "pct_full": ureg.Quantity(pct, ""),
            "V_total": ureg.Quantity(v_total * 7.48052, "gal")
        }


class VesselWallThickness(BaseEquation):
    """Calculate pressure vessel wall thickness (ASME Sec VIII Div 1)."""
    
    equation_id = "wall_thickness"
    name = "Vessel Wall Thickness"
    category = "Vessels"
    description = "Calculate minimum wall thickness per ASME Section VIII"
    reference = "ASME BPVC Section VIII Div 1"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P", "Design pressure", "pressure", "psig",
                              ParameterType.INPUT, symbol="P"),
            EquationParameter("D", "Inside diameter", "length", "in",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("S", "Allowable stress", "pressure", "psi",
                              ParameterType.INPUT, symbol="S",
                              tooltip="SA-516-70: 20000 psi, SS304: 20000 psi"),
            EquationParameter("E", "Joint efficiency", "dimensionless", "",
                              ParameterType.INPUT, symbol="E",
                              typical_range=(0.7, 1.0), tooltip="1.0=full RT, 0.85=spot RT, 0.7=no RT"),
            EquationParameter("CA", "Corrosion allowance", "length", "in",
                              ParameterType.INPUT, symbol="CA",
                              tooltip="Typically 0.0625 - 0.125 in"),
            EquationParameter("t_calc", "Calculated thickness", "length", "in",
                              ParameterType.OUTPUT, symbol="t_calc"),
            EquationParameter("t_min", "Minimum thickness (incl CA)", "length", "in",
                              ParameterType.OUTPUT, symbol="t_min"),
        ]
    
    def _calculate(self, P: pint.Quantity, D: pint.Quantity, S: pint.Quantity,
                   E: pint.Quantity, CA: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p = P.to('psi').magnitude
        d = D.to('in').magnitude
        s = S.to('psi').magnitude
        e = E.magnitude if hasattr(E, 'magnitude') else float(E)
        ca = CA.to('in').magnitude
        
        # ASME UG-27: t = P*R / (S*E - 0.6*P)
        r = d / 2
        t_calc = p * r / (s * e - 0.6 * p)
        t_min = t_calc + ca
        
        return {
            "t_calc": ureg.Quantity(t_calc, "in"),
            "t_min": ureg.Quantity(t_min, "in")
        }


class AgitatorPower(BaseEquation):
    """Calculate agitator power requirement."""
    
    equation_id = "agitator_power"
    name = "Agitator Power"
    category = "Vessels"
    description = "Calculate power for mechanical agitator"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("N", "Impeller speed", "dimensionless", "rpm",
                              ParameterType.INPUT, symbol="N"),
            EquationParameter("D", "Impeller diameter", "length", "in",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("rho", "Fluid density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρ"),
            EquationParameter("Np", "Power number (impeller type)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Np",
                              tooltip="Rushton: 5.0, Pitched blade: 1.5, Propeller: 0.35"),
            EquationParameter("P", "Power requirement", "power", "hp",
                              ParameterType.OUTPUT, symbol="P"),
            EquationParameter("Re", "Impeller Reynolds number", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Re"),
        ]
    
    def _calculate(self, N: pint.Quantity, D: pint.Quantity, rho: pint.Quantity,
                   Np: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        n = N.magnitude if hasattr(N, 'magnitude') else float(N)
        d = D.to('ft').magnitude
        rho_val = rho.to('lb/ft**3').magnitude
        np_val = Np.magnitude if hasattr(Np, 'magnitude') else float(Np)
        
        n_rps = n / 60  # Convert to rev/s
        
        # P = Np * rho * N³ * D⁵ / gc
        gc = 32.174
        p_ftlbf_s = np_val * rho_val * n_rps**3 * d**5 / gc
        p_hp = p_ftlbf_s / 550
        
        # Reynolds number (assuming water viscosity ≈ 1 cP)
        mu = 6.72e-4  # lb/(ft·s) for water
        re = rho_val * n_rps * d**2 / mu
        
        return {
            "P": ureg.Quantity(p_hp, "hp"),
            "Re": ureg.Quantity(re, "")
        }


class StokesSettling(BaseEquation):
    """Calculate particle settling velocity using Stokes' law."""
    
    equation_id = "stokes_settling"
    name = "Stokes' Law Settling"
    category = "Vessels"
    description = "Calculate terminal settling velocity for particles"
    reference = "Standard fluid mechanics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("d_p", "Particle diameter", "length", "micron",
                              ParameterType.INPUT, symbol="dₚ"),
            EquationParameter("rho_p", "Particle density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρₚ"),
            EquationParameter("rho_f", "Fluid density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρf"),
            EquationParameter("mu", "Fluid viscosity", "viscosity_dynamic", "cP",
                              ParameterType.INPUT, symbol="μ"),
            EquationParameter("V_t", "Terminal settling velocity", "velocity", "ft/hr",
                              ParameterType.OUTPUT, symbol="Vt"),
            EquationParameter("Re_p", "Particle Reynolds number", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Rep",
                              tooltip="Stokes law valid for Rep < 1"),
        ]
    
    def _calculate(self, d_p: pint.Quantity, rho_p: pint.Quantity, rho_f: pint.Quantity,
                   mu: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        dp = d_p.to('ft').magnitude
        rho_particle = rho_p.to('lb/ft**3').magnitude
        rho_fluid = rho_f.to('lb/ft**3').magnitude
        mu_val = mu.to('lb/(ft*s)').magnitude
        
        g = 32.174  # ft/s²
        
        # Stokes' law: V_t = g * d² * (ρp - ρf) / (18 * μ)
        v_t = g * dp**2 * (rho_particle - rho_fluid) / (18 * mu_val)
        
        # Reynolds number
        re_p = rho_fluid * v_t * dp / mu_val
        
        return {
            "V_t": ureg.Quantity(v_t * 3600, "ft/hr"),
            "Re_p": ureg.Quantity(re_p, "")
        }


VESSELS_EQUATIONS = {
    'tank_vertical': TankVolumeVertical,
    'tank_horizontal': TankVolumeHorizontal,
    'wall_thickness': VesselWallThickness,
    'agitator_power': AgitatorPower,
    'stokes_settling': StokesSettling,
}
