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


class HTU_NTU(BaseEquation):
    """Height of Transfer Unit method for packed columns."""
    
    equation_id = "htu_ntu"
    name = "HTU-NTU Method"
    category = "Mass Transfer"
    description = "Calculate packed column height using HTU and NTU"
    reference = "Treybal, Mass Transfer Operations"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("y_in", "Inlet gas mole fraction", "dimensionless", "",
                              ParameterType.INPUT, symbol="yᵢₙ", typical_range=(0, 1)),
            EquationParameter("y_out", "Outlet gas mole fraction", "dimensionless", "",
                              ParameterType.INPUT, symbol="yₒᵤₜ", typical_range=(0, 1)),
            EquationParameter("y_eq", "Equilibrium mole fraction", "dimensionless", "",
                              ParameterType.INPUT, symbol="y*", typical_range=(0, 1)),
            EquationParameter("HTU", "Height of transfer unit", "m", "ft",
                              ParameterType.INPUT, symbol="HTU", typical_range=(0.3, 2)),
            EquationParameter("NTU", "Number of transfer units", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="NTU"),
            EquationParameter("Z", "Packing height", "m", "ft",
                              ParameterType.OUTPUT, symbol="Z"),
        ]
    
    def _calculate(self, y_in: pint.Quantity, y_out: pint.Quantity, y_eq: pint.Quantity,
                   HTU: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        yin = y_in.magnitude if hasattr(y_in, 'magnitude') else float(y_in)
        yout = y_out.magnitude if hasattr(y_out, 'magnitude') else float(y_out)
        yeq = y_eq.magnitude if hasattr(y_eq, 'magnitude') else float(y_eq)
        htu = HTU.to('m').magnitude if hasattr(HTU, 'magnitude') else float(HTU)
        
        # Log mean driving force
        ntu = np.log((yin - yeq) / (yout - yeq))
        z = htu * ntu
        
        return {
            "NTU": ureg.Quantity(ntu, ""),
            "Z": ureg.Quantity(z, "m")
        }


class FilmCoefficient(BaseEquation):
    """Mass transfer film coefficient from correlations."""
    
    equation_id = "film_coefficient"
    name = "Film Mass Transfer Coefficient"
    category = "Mass Transfer"
    description = "Calculate gas or liquid film coefficient (kG or kL)"
    reference = "Sherwood, Pigford & Wilke"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Sh", "Sherwood number", "dimensionless", "",
                              ParameterType.INPUT, symbol="Sh", typical_range=(10, 1000)),
            EquationParameter("D_AB", "Diffusivity", "m^2/s", "ft^2/s",
                              ParameterType.INPUT, symbol="Dₐᵦ", typical_range=(1e-10, 1e-4)),
            EquationParameter("L", "Characteristic length", "m", "ft",
                              ParameterType.INPUT, symbol="L", typical_range=(0.001, 1)),
            EquationParameter("k", "Mass transfer coefficient", "m/s", "ft/s",
                              ParameterType.OUTPUT, symbol="k"),
        ]
    
    def _calculate(self, Sh: pint.Quantity, D_AB: pint.Quantity, L: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        sh = Sh.magnitude if hasattr(Sh, 'magnitude') else float(Sh)
        d = D_AB.to('m^2/s').magnitude if hasattr(D_AB, 'magnitude') else float(D_AB)
        l = L.to('m').magnitude if hasattr(L, 'magnitude') else float(L)
        
        k = sh * d / l
        return {"k": ureg.Quantity(k, "m/s")}


class LiquidExtraction(BaseEquation):
    """Liquid-liquid extraction stage calculation."""
    
    equation_id = "liquid_extraction"
    name = "Liquid-Liquid Extraction"
    category = "Mass Transfer"
    description = "Calculate extraction stages using distribution coefficient"
    reference = "Treybal, Liquid Extraction"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x_feed", "Feed solute concentration", "dimensionless", "",
                              ParameterType.INPUT, symbol="xF", typical_range=(0, 0.5)),
            EquationParameter("x_raffinate", "Raffinate solute concentration", "dimensionless", "",
                              ParameterType.INPUT, symbol="xR", typical_range=(0, 0.5)),
            EquationParameter("K_D", "Distribution coefficient (y/x at equilibrium)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Kᴅ", typical_range=(0.5, 10)),
            EquationParameter("S_F", "Solvent to feed ratio", "dimensionless", "",
                              ParameterType.INPUT, symbol="S/F", typical_range=(0.5, 5)),
            EquationParameter("N_stages", "Number of ideal stages", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="N"),
        ]
    
    def _calculate(self, x_feed: pint.Quantity, x_raffinate: pint.Quantity,
                   K_D: pint.Quantity, S_F: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        xf = x_feed.magnitude if hasattr(x_feed, 'magnitude') else float(x_feed)
        xr = x_raffinate.magnitude if hasattr(x_raffinate, 'magnitude') else float(x_raffinate)
        kd = K_D.magnitude if hasattr(K_D, 'magnitude') else float(K_D)
        sf = S_F.magnitude if hasattr(S_F, 'magnitude') else float(S_F)
        
        E = kd * sf  # Extraction factor
        if abs(E - 1) < 0.01:
            n = (xf - xr) / xf
        else:
            n = np.log((xf/xr * (E - 1) + 1) / E) / np.log(E)
        
        return {"N_stages": ureg.Quantity(max(n, 1), "")}


class MembranePermeation(BaseEquation):
    """Membrane gas permeation calculation."""
    
    equation_id = "membrane_permeation"
    name = "Membrane Permeation"
    category = "Mass Transfer"
    description = "Calculate membrane flux and selectivity"
    reference = "Baker, Membrane Technology"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P", "Permeability coefficient", "barrer", "",
                              ParameterType.INPUT, symbol="P", typical_range=(1, 10000)),
            EquationParameter("thickness", "Membrane thickness", "micrometer", "mil",
                              ParameterType.INPUT, symbol="δ", typical_range=(0.1, 100)),
            EquationParameter("delta_p", "Pressure difference", "bar", "psi",
                              ParameterType.INPUT, symbol="ΔP", typical_range=(1, 50)),
            EquationParameter("area", "Membrane area", "m^2", "ft^2",
                              ParameterType.INPUT, symbol="A", typical_range=(1, 1000)),
            EquationParameter("flux", "Gas flux", "cm^3/s", "",
                              ParameterType.OUTPUT, symbol="J"),
        ]
    
    def _calculate(self, P: pint.Quantity, thickness: pint.Quantity, delta_p: pint.Quantity,
                   area: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        perm = P.magnitude if hasattr(P, 'magnitude') else float(P)  # barrer
        t = thickness.magnitude if hasattr(thickness, 'magnitude') else float(thickness)  # um
        dp = delta_p.magnitude if hasattr(delta_p, 'magnitude') else float(delta_p)  # bar
        a = area.magnitude if hasattr(area, 'magnitude') else float(area)  # m^2
        
        # 1 barrer = 10^-10 cm³(STP)·cm/(cm²·s·cmHg)
        # Convert: flux = P * A * dP / thickness
        flux = perm * 1e-10 * (a * 1e4) * (dp * 75.006) / (t * 1e-4)  # cm³/s
        
        return {"flux": ureg.Quantity(flux, "cm^3/s")}


class PackedColumnPressureDrop(BaseEquation):
    """Pressure drop in packed columns (Ergun equation)."""
    
    equation_id = "packed_column_dp"
    name = "Packed Column Pressure Drop"
    category = "Mass Transfer"
    description = "Ergun equation for packed bed pressure drop"
    reference = "Ergun, 1952"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("G", "Gas mass velocity", "kg/m^2/s", "lb/ft^2/s",
                              ParameterType.INPUT, symbol="G", typical_range=(0.1, 5)),
            EquationParameter("dp", "Packing diameter", "m", "in",
                              ParameterType.INPUT, symbol="dₚ", typical_range=(0.01, 0.1)),
            EquationParameter("epsilon", "Void fraction", "dimensionless", "",
                              ParameterType.INPUT, symbol="ε", typical_range=(0.3, 0.95)),
            EquationParameter("mu", "Gas viscosity", "Pa*s", "cP",
                              ParameterType.INPUT, symbol="μ", typical_range=(1e-5, 5e-5)),
            EquationParameter("rho", "Gas density", "kg/m^3", "lb/ft^3",
                              ParameterType.INPUT, symbol="ρ", typical_range=(0.5, 50)),
            EquationParameter("L", "Bed height", "m", "ft",
                              ParameterType.INPUT, symbol="L", typical_range=(1, 20)),
            EquationParameter("delta_P", "Pressure drop", "Pa", "inH2O",
                              ParameterType.OUTPUT, symbol="ΔP"),
        ]
    
    def _calculate(self, G: pint.Quantity, dp: pint.Quantity, epsilon: pint.Quantity,
                   mu: pint.Quantity, rho: pint.Quantity, L: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        g = G.to('kg/m^2/s').magnitude if hasattr(G, 'magnitude') else float(G)
        d = dp.to('m').magnitude if hasattr(dp, 'magnitude') else float(dp)
        e = epsilon.magnitude if hasattr(epsilon, 'magnitude') else float(epsilon)
        visc = mu.to('Pa*s').magnitude if hasattr(mu, 'magnitude') else float(mu)
        dens = rho.to('kg/m^3').magnitude if hasattr(rho, 'magnitude') else float(rho)
        length = L.to('m').magnitude if hasattr(L, 'magnitude') else float(L)
        
        v = g / dens  # superficial velocity
        # Ergun equation
        dp_per_l = (150 * visc * (1-e)**2 / (d**2 * e**3)) * v + \
                   (1.75 * dens * (1-e) / (d * e**3)) * v**2
        delta_p = dp_per_l * length
        
        return {"delta_P": ureg.Quantity(delta_p, "Pa")}


# Update registry
MASS_TRANSFER_EQUATIONS = {
    'mccabe_thiele': McCabeThiele,
    'kremser': Kremser,
    'htu_ntu': HTU_NTU,
    'film_coefficient': FilmCoefficient,
    'liquid_extraction': LiquidExtraction,
    'membrane_permeation': MembranePermeation,
    'packed_column_dp': PackedColumnPressureDrop,
}
