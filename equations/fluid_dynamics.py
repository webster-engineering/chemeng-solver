"""
Fluid Dynamics equations - pressure drop, pipe sizing, pump calculations.
Each equation includes detailed variable descriptions and typical units.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class ReynoldsNumber(BaseEquation):
    """Calculate Reynolds number to determine flow regime."""
    
    equation_id = "reynolds_number"
    name = "Reynolds Number"
    category = "Fluid Dynamics"
    description = "Dimensionless number indicating laminar (<2100) or turbulent (>4000) flow"
    reference = "Fundamentals of Fluid Mechanics"
    
    derivation = """
**Reynolds Number Derivation**

The Reynolds number is a dimensionless quantity that predicts flow patterns. It represents the ratio of inertial forces to viscous forces:

$$Re = \\frac{\\text{Inertial Forces}}{\\text{Viscous Forces}} = \\frac{\\rho V D}{\\mu}$$

**Physical Meaning:**
- **Inertial forces** (ρV²) - tendency of fluid to keep moving
- **Viscous forces** (μV/D) - fluid's resistance to deformation

**Flow Regime Classification:**
- Re < 2100: **Laminar flow** - smooth, orderly layers
- 2100 < Re < 4000: **Transition** - unstable, intermittent
- Re > 4000: **Turbulent flow** - chaotic, mixing eddies

Named after Osborne Reynolds who demonstrated the transition in 1883 using dye injection experiments.
"""

    examples = [
        {
            "title": "Water in 4-inch Pipe",
            "description": "Water at 68°F flowing at 5 ft/s through a 4-inch schedule 40 pipe",
            "inputs": {"rho": 62.4, "V": 5, "D": 4.026, "mu": 1.0},
            "expected": {"Re": 173420},
            "conclusion": "Re >> 4000, so flow is fully turbulent"
        },
        {
            "title": "Oil in Small Tube",
            "description": "Heavy oil (viscosity 500 cP) at 2 ft/s in 1-inch tube",
            "inputs": {"rho": 54, "V": 2, "D": 1, "mu": 500},
            "expected": {"Re": 56},
            "conclusion": "Re << 2100, so flow is laminar"
        }
    ]
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("rho", "Fluid density - mass per unit volume of the fluid", 
                              "density", "lb/ft**3", ParameterType.INPUT, symbol="ρ",
                              tooltip="Water ≈ 62.4 lb/ft³, Air ≈ 0.075 lb/ft³"),
            EquationParameter("V", "Flow velocity - average fluid speed through pipe", 
                              "velocity", "ft/s", ParameterType.INPUT, symbol="V",
                              tooltip="Typical: 3-10 ft/s for liquids, 30-100 ft/s for gases"),
            EquationParameter("D", "Pipe inside diameter - internal diameter of the pipe", 
                              "length", "in", ParameterType.INPUT, symbol="D",
                              tooltip="Use schedule tables for pipe ID"),
            EquationParameter("mu", "Dynamic viscosity - fluid resistance to shear", 
                              "viscosity_dynamic", "cP", ParameterType.INPUT, symbol="μ",
                              tooltip="Water@68°F ≈ 1 cP, Honey ≈ 2000-10000 cP"),
            EquationParameter("Re", "Reynolds number - dimensionless flow regime indicator", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Re",
                              tooltip="<2100: Laminar, 2100-4000: Transition, >4000: Turbulent"),
        ]
    
    def _calculate(self, rho: pint.Quantity, V: pint.Quantity, D: pint.Quantity,
                   mu: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        rho_val = rho.to('kg/m**3').magnitude
        v_val = V.to('m/s').magnitude
        d_val = D.to('m').magnitude
        mu_val = mu.to('Pa*s').magnitude
        
        re = rho_val * v_val * d_val / mu_val
        
        return {"Re": ureg.Quantity(re, "")}


class BernoulliEquation(BaseEquation):
    """Bernoulli equation for ideal fluid flow energy balance."""
    
    equation_id = "bernoulli"
    name = "Bernoulli Equation"
    category = "Fluid Dynamics"
    description = "Energy balance between two points in a flowing fluid (no friction)"
    reference = "Fundamentals of Fluid Mechanics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Pressure at point 1 - upstream pressure", 
                              "pressure", "psi", ParameterType.INPUT, symbol="P₁"),
            EquationParameter("V1", "Velocity at point 1 - upstream fluid velocity", 
                              "velocity", "ft/s", ParameterType.INPUT, symbol="V₁"),
            EquationParameter("z1", "Elevation at point 1 - height above reference datum", 
                              "length", "ft", ParameterType.INPUT, symbol="z₁"),
            EquationParameter("P2", "Pressure at point 2 - downstream pressure", 
                              "pressure", "psi", ParameterType.INPUT, symbol="P₂"),
            EquationParameter("z2", "Elevation at point 2 - height above reference datum", 
                              "length", "ft", ParameterType.INPUT, symbol="z₂"),
            EquationParameter("rho", "Fluid density", "density", "lb/ft**3", 
                              ParameterType.INPUT, symbol="ρ"),
            EquationParameter("V2", "Velocity at point 2 - downstream fluid velocity", 
                              "velocity", "ft/s", ParameterType.OUTPUT, symbol="V₂"),
        ]
    
    def _calculate(self, P1: pint.Quantity, V1: pint.Quantity, z1: pint.Quantity,
                   P2: pint.Quantity, z2: pint.Quantity, rho: pint.Quantity, 
                   **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('Pa').magnitude
        v1 = V1.to('m/s').magnitude
        z1_val = z1.to('m').magnitude
        p2 = P2.to('Pa').magnitude
        z2_val = z2.to('m').magnitude
        rho_val = rho.to('kg/m**3').magnitude
        g = 9.81
        
        # P1/ρg + V1²/2g + z1 = P2/ρg + V2²/2g + z2
        v2_squared = 2 * ((p1 - p2)/(rho_val) + 0.5*v1**2 + g*(z1_val - z2_val))
        v2 = np.sqrt(max(v2_squared, 0))
        
        return {"V2": ureg.Quantity(v2 * 3.28084, "ft/s")}  # Convert m/s to ft/s


class OrificeFlow(BaseEquation):
    """Orifice plate flow measurement calculation."""
    
    equation_id = "orifice_flow"
    name = "Orifice Plate Flow"
    category = "Fluid Dynamics"
    description = "Calculate volumetric flow through an orifice plate"
    reference = "ISA/ASME Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Cd", "Discharge coefficient - accounts for real flow effects", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Cᵈ",
                              typical_range=(0.6, 0.65), tooltip="Typically 0.60-0.65 for sharp-edge orifice"),
            EquationParameter("d", "Orifice bore diameter - diameter of the orifice hole", 
                              "length", "in", ParameterType.INPUT, symbol="d"),
            EquationParameter("dP", "Differential pressure - pressure drop across orifice", 
                              "pressure", "inH2O", ParameterType.INPUT, symbol="ΔP",
                              tooltip="Measured by differential pressure transmitter"),
            EquationParameter("rho", "Fluid density at flowing conditions", 
                              "density", "lb/ft**3", ParameterType.INPUT, symbol="ρ"),
            EquationParameter("Q", "Volumetric flow rate", "volumetric_flow", "gpm", 
                              ParameterType.OUTPUT, symbol="Q"),
        ]
    
    def _calculate(self, Cd: pint.Quantity, d: pint.Quantity, dP: pint.Quantity,
                   rho: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        cd = Cd.magnitude if hasattr(Cd, 'magnitude') else float(Cd)
        d_ft = d.to('ft').magnitude
        dp_psf = dP.to('lbf/ft**2').magnitude
        rho_val = rho.to('lb/ft**3').magnitude
        
        area = np.pi * (d_ft/2)**2
        velocity = np.sqrt(2 * dp_psf * 32.174 / rho_val)
        q_cfs = cd * area * velocity
        q_gpm = q_cfs * 448.831
        
        return {"Q": ureg.Quantity(q_gpm, "gpm")}


class PipeFrictionFactor(BaseEquation):
    """Calculate Darcy friction factor using Swamee-Jain approximation."""
    
    equation_id = "friction_factor"
    name = "Pipe Friction Factor (Swamee-Jain)"
    category = "Fluid Dynamics"
    description = "Calculate Darcy friction factor for turbulent flow (Re > 4000)"
    reference = "Swamee & Jain, 1976"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("epsilon", "Pipe roughness - absolute surface roughness", 
                              "length", "in", ParameterType.INPUT, symbol="ε",
                              tooltip="Steel: 0.0018 in, PVC: 0.00006 in, Cast iron: 0.01 in"),
            EquationParameter("D", "Pipe inside diameter", "length", "in", 
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("Re", "Reynolds number - must be > 4000 for turbulent flow", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Re",
                              typical_range=(4000, 1e8)),
            EquationParameter("f", "Darcy friction factor - dimensionless resistance coefficient", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="f",
                              tooltip="Use in Darcy-Weisbach equation"),
        ]
    
    def _calculate(self, epsilon: pint.Quantity, D: pint.Quantity,
                   Re: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        eps = epsilon.to('in').magnitude
        d = D.to('in').magnitude
        re = Re.magnitude if hasattr(Re, 'magnitude') else float(Re)
        
        relative_roughness = eps / d
        
        # Swamee-Jain equation (explicit approximation of Colebrook-White)
        f = 0.25 / (np.log10(relative_roughness/3.7 + 5.74/re**0.9))**2
        
        return {"f": ureg.Quantity(f, "")}


class DarcyWeisbach(BaseEquation):
    """Darcy-Weisbach equation for pressure drop in pipes."""
    
    equation_id = "darcy_weisbach"
    name = "Darcy-Weisbach Pressure Drop"
    category = "Fluid Dynamics"
    description = "Calculate pressure drop in pipes using friction factor"
    reference = "Perry's Chemical Engineers' Handbook, 8th Ed"
    
    derivation = """
**Darcy-Weisbach Equation Derivation**

The Darcy-Weisbach equation is derived from dimensional analysis and energy balance principles:

$$\\Delta P = f \\cdot \\frac{L}{D} \\cdot \\frac{\\rho V^2}{2}$$

Or expressed as head loss:
$$h_f = f \\cdot \\frac{L}{D} \\cdot \\frac{V^2}{2g}$$

**Origin:**
1. **Force Balance**: Shear stress at pipe wall balanced against pressure difference
2. **Dimensional Analysis**: π-theorem reveals dimensionless groups
3. **The friction factor (f)** encapsulates all Reynolds and roughness effects

**Key Points:**
- Works for ANY Newtonian fluid (unlike Hazen-Williams)
- Requires friction factor from Moody Chart, Colebrook, or Swamee-Jain equation
- For laminar flow: f = 64/Re (exact solution)
- For turbulent flow: f depends on Re and relative roughness ε/D
"""

    examples = [
        {
            "title": "Water through 100 ft of 4-inch Steel Pipe",
            "description": "Water at 5 ft/s through 4\" Sch 40 steel pipe (f = 0.018)",
            "inputs": {"f": 0.018, "L": 100, "D": 4.026, "V": 5, "rho": 62.4},
            "expected": {"dP": 0.46},
            "steps": [
                "Step 1: Convert D to ft: 4.026 in ÷ 12 = 0.3355 ft",
                "Step 2: Calculate L/D = 100 / 0.3355 = 298.1",
                "Step 3: Calculate velocity head: ρV²/2g = 62.4 × 25 / 64.35 = 24.23 lbf/ft²",
                "Step 4: ΔP = f × (L/D) × (ρV²/2g) = 0.018 × 298.1 × 24.23 = 130 lbf/ft²",
                "Step 5: Convert to psi: 130 ÷ 144 = 0.90 psi"
            ],
            "conclusion": "Pressure drop of ~0.9 psi per 100 ft is reasonable for water"
        }
    ]
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("f", "Darcy friction factor - from Moody chart or Colebrook equation", 
                              "dimensionless", "", ParameterType.INPUT, symbol="f", 
                              typical_range=(0.001, 0.1),
                              tooltip="Use Friction Factor equation to calculate"),
            EquationParameter("L", "Pipe length - total length of pipe section", 
                              "length", "ft", ParameterType.INPUT, symbol="L"),
            EquationParameter("D", "Pipe inside diameter - internal diameter from pipe schedule", 
                              "length", "in", ParameterType.INPUT, symbol="D"),
            EquationParameter("V", "Flow velocity - average velocity = Q/A", 
                              "velocity", "ft/s", ParameterType.INPUT, symbol="V",
                              tooltip="Calculate from V = Q/(πD²/4)"),
            EquationParameter("rho", "Fluid density - mass per unit volume", 
                              "density", "lb/ft**3", ParameterType.INPUT, symbol="ρ"),
            EquationParameter("dP", "Pressure drop - pressure loss due to friction", 
                              "pressure", "psi", ParameterType.OUTPUT, symbol="ΔP"),
        ]
    
    def _calculate(self, f: pint.Quantity, L: pint.Quantity, D: pint.Quantity,
                   V: pint.Quantity, rho: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        f_val = f.magnitude if hasattr(f, 'magnitude') else float(f)
        L_val = L.to('ft').magnitude
        D_val = D.to('ft').magnitude
        V_val = V.to('ft/s').magnitude
        rho_val = rho.to('lb/ft**3').magnitude
        
        g = 32.174
        dp_lbf_ft2 = f_val * (L_val / D_val) * (rho_val * V_val**2) / (2 * g)
        dp_psi = dp_lbf_ft2 / 144
        
        return {"dP": ureg.Quantity(dp_psi, "psi")}


class HazenWilliams(BaseEquation):
    """Hazen-Williams equation for water flow pressure loss."""
    
    equation_id = "hazen_williams"
    name = "Hazen-Williams Head Loss"
    category = "Fluid Dynamics"
    description = "Calculate head loss for water flow in pipes (empirical equation)"
    reference = "Fire protection and water distribution standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Volumetric flow rate - volume per unit time", 
                              "volumetric_flow", "gpm", ParameterType.INPUT, symbol="Q"),
            EquationParameter("C", "Hazen-Williams C-factor - pipe roughness coefficient", 
                              "dimensionless", "", ParameterType.INPUT, symbol="C", 
                              typical_range=(60, 150),
                              tooltip="New steel:140, Old steel:100, Cast iron:130, PVC:150"),
            EquationParameter("D", "Pipe inside diameter", "length", "in", 
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("L", "Pipe length", "length", "ft", 
                              ParameterType.INPUT, symbol="L"),
            EquationParameter("hL", "Head loss - energy loss expressed as fluid height", 
                              "length", "ft", ParameterType.OUTPUT, symbol="hₗ",
                              tooltip="Multiply by 0.433×SG to convert to psi"),
        ]
    
    def _calculate(self, Q: pint.Quantity, C: pint.Quantity, D: pint.Quantity,
                   L: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gpm').magnitude
        c = C.magnitude if hasattr(C, 'magnitude') else float(C)
        d = D.to('in').magnitude
        length = L.to('ft').magnitude
        
        hl = 0.002083 * length * (100/c)**1.852 * (q**1.852) / (d**4.8655)
        
        return {"hL": ureg.Quantity(hl, "ft")}


class PumpPower(BaseEquation):
    """Calculate pump hydraulic and brake horsepower."""
    
    equation_id = "pump_power"
    name = "Pump Power Calculation"
    category = "Fluid Dynamics"
    description = "Calculate hydraulic power and motor power required for pumping"
    reference = "Hydraulic Institute Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Volumetric flow rate", "volumetric_flow", "gpm", 
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("H", "Total dynamic head - total head the pump must overcome", 
                              "length", "ft", ParameterType.INPUT, symbol="H",
                              tooltip="From pump head calculation"),
            EquationParameter("SG", "Specific gravity - density relative to water at 60°F", 
                              "dimensionless", "", ParameterType.INPUT, symbol="SG",
                              tooltip="Water=1.0, Gasoline≈0.72, Brine≈1.2"),
            EquationParameter("eta_pump", "Pump efficiency - hydraulic efficiency of pump", 
                              "dimensionless", "", ParameterType.INPUT, symbol="ηₚ",
                              typical_range=(0.5, 0.9), tooltip="Typically 0.6-0.85"),
            EquationParameter("eta_motor", "Motor efficiency", "dimensionless", "", 
                              ParameterType.INPUT, symbol="ηₘ", typical_range=(0.85, 0.97),
                              tooltip="Typically 0.90-0.96"),
            EquationParameter("WHP", "Water (hydraulic) horsepower - ideal power to fluid", 
                              "power", "hp", ParameterType.OUTPUT, symbol="WHP"),
            EquationParameter("BHP", "Brake horsepower - power required at pump shaft", 
                              "power", "hp", ParameterType.OUTPUT, symbol="BHP"),
            EquationParameter("Motor_HP", "Motor horsepower - electrical power input", 
                              "power", "hp", ParameterType.OUTPUT, symbol="Motor HP"),
        ]
    
    def _calculate(self, Q: pint.Quantity, H: pint.Quantity, SG: pint.Quantity,
                   eta_pump: pint.Quantity, eta_motor: pint.Quantity, 
                   **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gpm').magnitude
        h = H.to('ft').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        eff_p = eta_pump.magnitude if hasattr(eta_pump, 'magnitude') else float(eta_pump)
        eff_m = eta_motor.magnitude if hasattr(eta_motor, 'magnitude') else float(eta_motor)
        
        whp = q * h * sg / 3960
        bhp = whp / eff_p
        motor_hp = bhp / eff_m
        
        return {
            "WHP": ureg.Quantity(whp, "hp"),
            "BHP": ureg.Quantity(bhp, "hp"),
            "Motor_HP": ureg.Quantity(motor_hp, "hp")
        }


class PumpHead(BaseEquation):
    """Calculate total dynamic head for pump sizing."""
    
    equation_id = "pump_head"
    name = "Pump Total Dynamic Head"
    category = "Fluid Dynamics"
    description = "Calculate TDH = static head + friction head + pressure head"
    reference = "Hydraulic Institute Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("h_static", "Static head - elevation difference (discharge - suction)", 
                              "length", "ft", ParameterType.INPUT, symbol="hₛ",
                              tooltip="Positive if pumping upward"),
            EquationParameter("h_friction", "Friction head loss - total losses in suction and discharge piping", 
                              "length", "ft", ParameterType.INPUT, symbol="hf",
                              tooltip="From Darcy-Weisbach or Hazen-Williams"),
            EquationParameter("P_discharge", "Discharge pressure - required pressure at discharge point", 
                              "pressure", "psi", ParameterType.INPUT, symbol="Pᵈ"),
            EquationParameter("P_suction", "Suction pressure - pressure at pump suction", 
                              "pressure", "psi", ParameterType.INPUT, symbol="Pₛ"),
            EquationParameter("SG", "Specific gravity of fluid", "dimensionless", "", 
                              ParameterType.INPUT, symbol="SG"),
            EquationParameter("TDH", "Total dynamic head - total head pump must provide", 
                              "length", "ft", ParameterType.OUTPUT, symbol="TDH"),
        ]
    
    def _calculate(self, h_static: pint.Quantity, h_friction: pint.Quantity,
                   P_discharge: pint.Quantity, P_suction: pint.Quantity,
                   SG: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        hs = h_static.to('ft').magnitude
        hf = h_friction.to('ft').magnitude
        pd = P_discharge.to('psi').magnitude
        ps = P_suction.to('psi').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        
        h_pressure = (pd - ps) * 2.31 / sg
        tdh = hs + hf + h_pressure
        
        return {"TDH": ureg.Quantity(tdh, "ft")}


class NPSH(BaseEquation):
    """Net Positive Suction Head available calculation."""
    
    equation_id = "npsh"
    name = "NPSHa (Available)"
    category = "Fluid Dynamics"
    description = "Calculate available NPSH to compare against pump NPSHr"
    reference = "Hydraulic Institute Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P_atm", "Atmospheric pressure - local barometric pressure", 
                              "pressure", "psia", ParameterType.INPUT, symbol="Pₐₜₘ", 
                              typical_range=(10, 16),
                              tooltip="Sea level ≈ 14.7 psia, 5000 ft ≈ 12.2 psia"),
            EquationParameter("P_vapor", "Vapor pressure - saturation pressure of liquid at temperature", 
                              "pressure", "psia", ParameterType.INPUT, symbol="Pᵥ",
                              tooltip="Water@68°F=0.34 psia, @150°F=3.7 psia, @212°F=14.7 psia"),
            EquationParameter("h_s", "Static suction head - height of liquid above pump centerline", 
                              "length", "ft", ParameterType.INPUT, symbol="hₛ",
                              tooltip="Positive if liquid is above pump, negative if below"),
            EquationParameter("h_f", "Suction friction losses - head loss in suction piping", 
                              "length", "ft", ParameterType.INPUT, symbol="hf"),
            EquationParameter("SG", "Specific gravity of fluid", "dimensionless", "", 
                              ParameterType.INPUT, symbol="SG"),
            EquationParameter("NPSHa", "Available NPSH - must exceed pump NPSHr by margin", 
                              "length", "ft", ParameterType.OUTPUT, symbol="NPSHₐ",
                              tooltip="Should be > NPSHr + 2-3 ft margin"),
        ]
    
    def _calculate(self, P_atm: pint.Quantity, P_vapor: pint.Quantity,
                   h_s: pint.Quantity, h_f: pint.Quantity,
                   SG: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        patm = P_atm.to('psia').magnitude
        pvap = P_vapor.to('psia').magnitude
        hs = h_s.to('ft').magnitude
        hf = h_f.to('ft').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        
        h_atm = patm * 2.31 / sg
        h_vap = pvap * 2.31 / sg
        npsha = h_atm + hs - hf - h_vap
        
        return {"NPSHa": ureg.Quantity(npsha, "ft")}


FLUID_DYNAMICS_EQUATIONS = {
    'reynolds_number': ReynoldsNumber,
    'bernoulli': BernoulliEquation,
    'orifice_flow': OrificeFlow,
    'friction_factor': PipeFrictionFactor,
    'darcy_weisbach': DarcyWeisbach,
    'hazen_williams': HazenWilliams,
    'pump_power': PumpPower,
    'pump_head': PumpHead,
    'npsh': NPSH,
}
