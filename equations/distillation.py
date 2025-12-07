"""
Distillation equations - column design, reflux, stages.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class FenskeMinStages(BaseEquation):
    """Fenske equation for minimum theoretical stages at total reflux."""
    
    equation_id = "fenske"
    name = "Fenske Minimum Stages"
    category = "Distillation"
    description = "Calculate minimum stages at total reflux using Fenske equation"
    reference = "Fenske, 1932"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x_D", "Light key mole fraction in distillate", "dimensionless", "",
                              ParameterType.INPUT, symbol="xD", typical_range=(0.9, 0.999)),
            EquationParameter("x_B", "Light key mole fraction in bottoms", "dimensionless", "",
                              ParameterType.INPUT, symbol="xB", typical_range=(0.001, 0.1)),
            EquationParameter("alpha", "Average relative volatility", "dimensionless", "",
                              ParameterType.INPUT, symbol="α",
                              tooltip="α = K_light / K_heavy"),
            EquationParameter("N_min", "Minimum theoretical stages", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Nmin"),
        ]
    
    def _calculate(self, x_D: pint.Quantity, x_B: pint.Quantity, alpha: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        xd = x_D.magnitude if hasattr(x_D, 'magnitude') else float(x_D)
        xb = x_B.magnitude if hasattr(x_B, 'magnitude') else float(x_B)
        a = alpha.magnitude if hasattr(alpha, 'magnitude') else float(alpha)
        
        n_min = np.log((xd/(1-xd)) * ((1-xb)/xb)) / np.log(a)
        
        return {"N_min": ureg.Quantity(n_min, "")}


class UnderwoodMinReflux(BaseEquation):
    """Underwood equation for minimum reflux ratio."""
    
    equation_id = "underwood"
    name = "Underwood Minimum Reflux"
    category = "Distillation"
    description = "Calculate minimum reflux ratio using Underwood equations"
    reference = "Underwood, 1948"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("alpha", "Relative volatility of light key", "dimensionless", "",
                              ParameterType.INPUT, symbol="α"),
            EquationParameter("z_F", "Light key mole fraction in feed", "dimensionless", "",
                              ParameterType.INPUT, symbol="zF"),
            EquationParameter("x_D", "Light key mole fraction in distillate", "dimensionless", "",
                              ParameterType.INPUT, symbol="xD"),
            EquationParameter("q", "Feed quality (1=sat liquid, 0=sat vapor)", "dimensionless", "",
                              ParameterType.INPUT, symbol="q",
                              tooltip="1=saturated liquid, 0=saturated vapor"),
            EquationParameter("R_min", "Minimum reflux ratio", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Rmin"),
            EquationParameter("R_actual", "Suggested actual reflux (1.2×Rmin)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="R"),
        ]
    
    def _calculate(self, alpha: pint.Quantity, z_F: pint.Quantity, x_D: pint.Quantity,
                   q: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        a = alpha.magnitude if hasattr(alpha, 'magnitude') else float(alpha)
        zf = z_F.magnitude if hasattr(z_F, 'magnitude') else float(z_F)
        xd = x_D.magnitude if hasattr(x_D, 'magnitude') else float(x_D)
        q_val = q.magnitude if hasattr(q, 'magnitude') else float(q)
        
        # Simplified Underwood for binary mixture
        # θ satisfies: α*zf/(α-θ) + (1-zf)*1/(1-θ) = 1-q
        # For saturated liquid feed (q=1): θ is between 1 and α
        theta = (a + 1) / 2  # Initial guess
        
        # R_min from: α*xd/(α-θ) + (1-xd)*1/(1-θ) = R_min + 1
        r_min = a * xd / (a - theta) + (1 - xd) / (1 - theta) - 1
        r_actual = 1.2 * r_min
        
        return {
            "R_min": ureg.Quantity(max(r_min, 0.1), ""),
            "R_actual": ureg.Quantity(max(r_actual, 0.1), "")
        }


class GillilandCorrelation(BaseEquation):
    """Gilliland correlation for actual stages."""
    
    equation_id = "gilliland"
    name = "Gilliland Correlation"
    category = "Distillation"
    description = "Estimate actual stages from N_min and R_min"
    reference = "Gilliland, 1940"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("N_min", "Minimum stages (from Fenske)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Nmin"),
            EquationParameter("R_min", "Minimum reflux ratio", "dimensionless", "",
                              ParameterType.INPUT, symbol="Rmin"),
            EquationParameter("R", "Actual operating reflux ratio", "dimensionless", "",
                              ParameterType.INPUT, symbol="R"),
            EquationParameter("N", "Actual theoretical stages", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="N"),
            EquationParameter("efficiency", "Assumed tray efficiency", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="η",
                              tooltip="Typically 0.5-0.8 for trayed columns"),
            EquationParameter("N_actual", "Actual trays needed", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="N_actual"),
        ]
    
    def _calculate(self, N_min: pint.Quantity, R_min: pint.Quantity, R: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        n_min = N_min.magnitude if hasattr(N_min, 'magnitude') else float(N_min)
        r_min = R_min.magnitude if hasattr(R_min, 'magnitude') else float(R_min)
        r = R.magnitude if hasattr(R, 'magnitude') else float(R)
        
        # X = (R - R_min) / (R + 1)
        x = (r - r_min) / (r + 1)
        
        # Gilliland correlation (Eduljee form):
        # Y = (N - N_min) / (N + 1) = 0.75 * (1 - X^0.5668)
        y = 0.75 * (1 - x**0.5668)
        
        # Solve for N
        n = (y + n_min + y * n_min) / (1 - y)
        
        efficiency = 0.65  # Typical sieve tray efficiency
        n_actual = np.ceil(n / efficiency)
        
        return {
            "N": ureg.Quantity(n, ""),
            "efficiency": ureg.Quantity(efficiency, ""),
            "N_actual": ureg.Quantity(n_actual, "")
        }


class FeedStageLocation(BaseEquation):
    """Kirkbride equation for optimal feed stage location."""
    
    equation_id = "kirkbride"
    name = "Feed Stage Location (Kirkbride)"
    category = "Distillation"
    description = "Estimate optimal feed tray location"
    reference = "Kirkbride, 1944"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("N", "Total theoretical stages", "dimensionless", "",
                              ParameterType.INPUT, symbol="N"),
            EquationParameter("B_D", "Bottoms to distillate molar ratio", "dimensionless", "",
                              ParameterType.INPUT, symbol="B/D"),
            EquationParameter("x_F_HK", "Heavy key mole fraction in feed", "dimensionless", "",
                              ParameterType.INPUT, symbol="xF,HK"),
            EquationParameter("x_F_LK", "Light key mole fraction in feed", "dimensionless", "",
                              ParameterType.INPUT, symbol="xF,LK"),
            EquationParameter("x_D_HK", "Heavy key mole fraction in distillate", "dimensionless", "",
                              ParameterType.INPUT, symbol="xD,HK"),
            EquationParameter("x_B_LK", "Light key mole fraction in bottoms", "dimensionless", "",
                              ParameterType.INPUT, symbol="xB,LK"),
            EquationParameter("N_R", "Rectifying section stages (above feed)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="NR"),
            EquationParameter("N_S", "Stripping section stages (below feed)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="NS"),
            EquationParameter("feed_stage", "Optimal feed stage (from bottom)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Feed"),
        ]
    
    def _calculate(self, N: pint.Quantity, B_D: pint.Quantity, x_F_HK: pint.Quantity,
                   x_F_LK: pint.Quantity, x_D_HK: pint.Quantity, x_B_LK: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        n = N.magnitude if hasattr(N, 'magnitude') else float(N)
        bd = B_D.magnitude if hasattr(B_D, 'magnitude') else float(B_D)
        xfhk = x_F_HK.magnitude if hasattr(x_F_HK, 'magnitude') else float(x_F_HK)
        xflk = x_F_LK.magnitude if hasattr(x_F_LK, 'magnitude') else float(x_F_LK)
        xdhk = x_D_HK.magnitude if hasattr(x_D_HK, 'magnitude') else float(x_D_HK)
        xblk = x_B_LK.magnitude if hasattr(x_B_LK, 'magnitude') else float(x_B_LK)
        
        # Kirkbride equation: N_R/N_S = [(B/D) * (x_F,HK/x_F,LK) * (x_B,LK/x_D,HK)²]^0.206
        ratio = (bd * (xfhk/xflk) * (xblk/xdhk)**2)**0.206
        
        n_s = n / (1 + ratio)
        n_r = n - n_s
        feed_stage = np.round(n_s)
        
        return {
            "N_R": ureg.Quantity(n_r, ""),
            "N_S": ureg.Quantity(n_s, ""),
            "feed_stage": ureg.Quantity(feed_stage, "")
        }


class TrayHydraulics(BaseEquation):
    """Check tray hydraulics - flooding and weeping."""
    
    equation_id = "tray_hydraulics"
    name = "Tray Flooding Check"
    category = "Distillation"
    description = "Calculate flooding velocity and percent flood"
    reference = "Fair correlation"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("V", "Vapor molar flow", "dimensionless", "lbmol/hr",
                              ParameterType.INPUT, symbol="V"),
            EquationParameter("MW_v", "Vapor molecular weight", "dimensionless", "",
                              ParameterType.INPUT, symbol="MW_v"),
            EquationParameter("rho_v", "Vapor density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρᵥ"),
            EquationParameter("rho_l", "Liquid density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρₗ"),
            EquationParameter("sigma", "Surface tension", "dimensionless", "dyne/cm",
                              ParameterType.INPUT, symbol="σ",
                              tooltip="Water ≈ 72, Organics ≈ 20-30"),
            EquationParameter("A_net", "Net tray area", "area", "ft**2",
                              ParameterType.INPUT, symbol="Anet"),
            EquationParameter("U_flood", "Flooding velocity", "velocity", "ft/s",
                              ParameterType.OUTPUT, symbol="Uflood"),
            EquationParameter("U_actual", "Actual vapor velocity", "velocity", "ft/s",
                              ParameterType.OUTPUT, symbol="Uact"),
            EquationParameter("pct_flood", "Percent of flooding", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="%Flood",
                              tooltip="Design: 70-85%"),
        ]
    
    def _calculate(self, V: pint.Quantity, MW_v: pint.Quantity, rho_v: pint.Quantity,
                   rho_l: pint.Quantity, sigma: pint.Quantity, A_net: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        v = V.magnitude if hasattr(V, 'magnitude') else float(V)
        mw = MW_v.magnitude if hasattr(MW_v, 'magnitude') else float(MW_v)
        rv = rho_v.to('lb/ft**3').magnitude
        rl = rho_l.to('lb/ft**3').magnitude
        sig = sigma.magnitude if hasattr(sigma, 'magnitude') else float(sigma)
        a = A_net.to('ft**2').magnitude
        
        # Fair correlation for Csb
        flv = (v * mw / 3600) * np.sqrt(rv / rl) / (rl * a)
        csb = 0.15  # Simplified, normally from chart
        
        # Flooding velocity
        u_flood = csb * (sig/20)**0.2 * np.sqrt((rl - rv) / rv)
        
        # Actual velocity
        v_mass = v * mw / 3600  # lb/s
        v_vol = v_mass / rv  # ft³/s
        u_actual = v_vol / a
        
        pct_flood = (u_actual / u_flood) * 100
        
        return {
            "U_flood": ureg.Quantity(u_flood, "ft/s"),
            "U_actual": ureg.Quantity(u_actual, "ft/s"),
            "pct_flood": ureg.Quantity(pct_flood, "")
        }


DISTILLATION_EQUATIONS = {
    'fenske': FenskeMinStages,
    'underwood': UnderwoodMinReflux,
    'gilliland': GillilandCorrelation,
    'kirkbride': FeedStageLocation,
    'tray_hydraulics': TrayHydraulics,
}
