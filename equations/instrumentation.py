"""
Instrumentation equations - orifice, transmitter, control valve sizing.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class OrificeDesign(BaseEquation):
    """Design orifice plate for flow measurement."""
    
    equation_id = "orifice_design"
    name = "Orifice Plate Sizing"
    category = "Instrumentation"
    description = "Size orifice bore for given flow and differential pressure"
    reference = "ISO 5167 / AGA 3"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Design flow rate", "volumetric_flow", "gpm",
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("D", "Pipe inside diameter", "length", "in",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("dP_design", "Design differential pressure", "pressure", "inH2O",
                              ParameterType.INPUT, symbol="ΔP",
                              tooltip="Typically 100-200 inH2O at max flow"),
            EquationParameter("rho", "Fluid density", "density", "lb/ft**3",
                              ParameterType.INPUT, symbol="ρ"),
            EquationParameter("d", "Orifice bore diameter", "length", "in",
                              ParameterType.OUTPUT, symbol="d"),
            EquationParameter("beta", "Beta ratio (d/D)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="β",
                              tooltip="Should be 0.2-0.75 for accuracy"),
            EquationParameter("perm_loss", "Permanent pressure loss (%)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Loss%"),
        ]
    
    def _calculate(self, Q: pint.Quantity, D: pint.Quantity, dP_design: pint.Quantity,
                   rho: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('ft**3/s').magnitude
        d_pipe = D.to('ft').magnitude
        dp = dP_design.to('lbf/ft**2').magnitude
        rho_val = rho.to('lb/ft**3').magnitude
        
        # Iterative solution for orifice diameter
        # Q = Cd * A * sqrt(2 * dP / rho), assume Cd = 0.61
        cd = 0.61
        a_needed = q / (cd * np.sqrt(2 * dp * 32.174 / rho_val))
        d_orifice = np.sqrt(4 * a_needed / np.pi)
        
        beta = d_orifice / d_pipe
        
        # Permanent pressure loss (approximate)
        perm_loss = (1 - beta**2) * 100
        
        return {
            "d": ureg.Quantity(d_orifice * 12, "in"),
            "beta": ureg.Quantity(beta, ""),
            "perm_loss": ureg.Quantity(perm_loss, "")
        }


class ThermowellWakeFrequency(BaseEquation):
    """Calculate thermowell wake frequency for vibration analysis."""
    
    equation_id = "thermowell_wake"
    name = "Thermowell Wake Frequency"
    category = "Instrumentation"
    description = "Check thermowell for flow-induced vibration (ASME PTC 19.3)"
    reference = "ASME PTC 19.3 TW-2016"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("V", "Fluid velocity", "velocity", "ft/s",
                              ParameterType.INPUT, symbol="V"),
            EquationParameter("d_tip", "Thermowell tip diameter", "length", "in",
                              ParameterType.INPUT, symbol="d",
                              tooltip="Typically 0.26-1.0 inches"),
            EquationParameter("L", "Unsupported length", "length", "in",
                              ParameterType.INPUT, symbol="L"),
            EquationParameter("f_wake", "Wake (Strouhal) frequency", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="fₛ (Hz)"),
            EquationParameter("f_natural", "Estimated natural frequency", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="fₙ (Hz)",
                              tooltip="Using simplified beam formula"),
            EquationParameter("freq_ratio", "Frequency ratio (should be <0.8)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="fₛ/fₙ"),
        ]
    
    def _calculate(self, V: pint.Quantity, d_tip: pint.Quantity, L: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        v = V.to('ft/s').magnitude
        d = d_tip.to('ft').magnitude
        length = L.to('ft').magnitude
        
        # Strouhal frequency: f = St * V / d, St ≈ 0.22
        st = 0.22
        f_wake = st * v / d
        
        # Simplified natural frequency estimation (cantilever beam)
        # f_n ≈ 3.52 / (2π) * sqrt(E*I / (ρ*A*L^4))
        # Simplified: f_n ≈ 1000 / L² for steel thermowell
        f_natural = 1000 / (length * 12)**1.5 * 100  # Rough estimate
        
        freq_ratio = f_wake / f_natural if f_natural > 0 else 0
        
        return {
            "f_wake": ureg.Quantity(f_wake, ""),
            "f_natural": ureg.Quantity(f_natural, ""),
            "freq_ratio": ureg.Quantity(freq_ratio, "")
        }


class TransmitterSpan(BaseEquation):
    """Calculate transmitter span and calibration range."""
    
    equation_id = "transmitter_span"
    name = "Transmitter Span Calculation"
    category = "Instrumentation"
    description = "Calculate LRV, URV, and span for 4-20mA transmitter"
    reference = "Standard instrumentation"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("min_process", "Minimum expected process value", "dimensionless", "",
                              ParameterType.INPUT, symbol="PV_min"),
            EquationParameter("max_process", "Maximum expected process value", "dimensionless", "",
                              ParameterType.INPUT, symbol="PV_max"),
            EquationParameter("margin", "Margin above/below range (%)", "dimensionless", "",
                              ParameterType.INPUT, symbol="Margin",
                              typical_range=(5, 25), tooltip="Typically 10-20%"),
            EquationParameter("LRV", "Lower range value (4 mA)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="LRV"),
            EquationParameter("URV", "Upper range value (20 mA)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="URV"),
            EquationParameter("span", "Calibrated span (URV - LRV)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Span"),
            EquationParameter("turndown", "Turndown ratio", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="TD"),
        ]
    
    def _calculate(self, min_process: pint.Quantity, max_process: pint.Quantity,
                   margin: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        pv_min = min_process.magnitude if hasattr(min_process, 'magnitude') else float(min_process)
        pv_max = max_process.magnitude if hasattr(max_process, 'magnitude') else float(max_process)
        m = margin.magnitude if hasattr(margin, 'magnitude') else float(margin)
        
        range_val = pv_max - pv_min
        margin_val = range_val * m / 100
        
        lrv = pv_min - margin_val
        urv = pv_max + margin_val
        span = urv - lrv
        td = urv / lrv if lrv != 0 else float('inf')
        
        return {
            "LRV": ureg.Quantity(lrv, ""),
            "URV": ureg.Quantity(urv, ""),
            "span": ureg.Quantity(span, ""),
            "turndown": ureg.Quantity(td, "")
        }


class ControlValveNoise(BaseEquation):
    """Estimate control valve aerodynamic noise."""
    
    equation_id = "valve_noise"
    name = "Control Valve Noise Prediction"
    category = "Instrumentation"
    description = "Estimate valve noise level for gas service"
    reference = "IEC 60534-8-3"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Inlet pressure", "pressure", "psia",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Outlet pressure", "pressure", "psia",
                              ParameterType.INPUT, symbol="P₂"),
            EquationParameter("W", "Mass flow rate", "mass_flow", "lb/hr",
                              ParameterType.INPUT, symbol="W"),
            EquationParameter("D", "Downstream pipe diameter", "length", "in",
                              ParameterType.INPUT, symbol="D"),
            EquationParameter("dB_A", "Estimated noise level at 1m (dBA)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="dBA"),
            EquationParameter("attenuation", "Pipe wall attenuation (dB)", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Atten"),
        ]
    
    def _calculate(self, P1: pint.Quantity, P2: pint.Quantity, W: pint.Quantity,
                   D: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('psia').magnitude
        p2 = P2.to('psia').magnitude
        w = W.to('lb/hr').magnitude
        d = D.to('in').magnitude
        
        # Simplified IEC 60534-8-3 correlation
        dp = p1 - p2
        pressure_ratio = p2 / p1
        
        # Basic noise estimate (simplified)
        if pressure_ratio > 0.5:  # Subcritical
            spl = 50 + 30 * np.log10(w/1000) + 40 * np.log10(dp/10)
        else:  # Critical/choked
            spl = 80 + 30 * np.log10(w/1000) + 20 * np.log10(p1/100)
        
        # Pipe wall attenuation (approximate)
        attenuation = 10 + 20 * np.log10(d/4)
        
        return {
            "dB_A": ureg.Quantity(min(spl, 120), ""),
            "attenuation": ureg.Quantity(attenuation, "")
        }


INSTRUMENTATION_EQUATIONS = {
    'orifice_design': OrificeDesign,
    'thermowell_wake': ThermowellWakeFrequency,
    'transmitter_span': TransmitterSpan,
    'valve_noise': ControlValveNoise,
}


class RTDTemperature(BaseEquation):
    """RTD temperature from resistance measurement."""
    
    equation_id = "rtd_temperature"
    name = "RTD Temperature Calculation"
    category = "Instrumentation"
    description = "Calculate temperature from RTD resistance (Pt100/Pt1000)"
    reference = "IEC 60751"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("R", "Measured resistance", "ohm", "",
                              ParameterType.INPUT, symbol="R"),
            EquationParameter("R0", "Reference resistance at 0°C (100 for Pt100)", "ohm", "",
                              ParameterType.INPUT, symbol="R₀"),
            EquationParameter("alpha", "Temperature coefficient (0.00385 typical)", "1/degC", "",
                              ParameterType.INPUT, symbol="α", typical_range=(0.003, 0.004)),
            EquationParameter("T", "Temperature", "temperature", "degC",
                              ParameterType.OUTPUT, symbol="T"),
        ]
    
    def _calculate(self, R: pint.Quantity, R0: pint.Quantity, alpha: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        r = R.magnitude if hasattr(R, 'magnitude') else float(R)
        r0 = R0.magnitude if hasattr(R0, 'magnitude') else float(R0)
        a = alpha.magnitude if hasattr(alpha, 'magnitude') else float(alpha)
        
        # Linear approximation: R = R0*(1 + alpha*T)
        t = (r/r0 - 1) / a
        
        return {"T": ureg.Quantity(t, "degC")}


class SignalConversion(BaseEquation):
    """Convert between 4-20mA, 0-10V, and engineering units."""
    
    equation_id = "signal_conversion"
    name = "Signal Conversion (4-20mA)"
    category = "Instrumentation"
    description = "Convert between mA signal and engineering units"
    reference = "Standard instrumentation"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("mA", "Current signal", "mA", "",
                              ParameterType.INPUT, symbol="mA", typical_range=(4, 20)),
            EquationParameter("LRV", "Lower range value (at 4mA)", "dimensionless", "",
                              ParameterType.INPUT, symbol="LRV"),
            EquationParameter("URV", "Upper range value (at 20mA)", "dimensionless", "",
                              ParameterType.INPUT, symbol="URV"),
            EquationParameter("PV", "Process value in engineering units", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="PV"),
            EquationParameter("percent", "Percent of span", "dimensionless", "%",
                              ParameterType.OUTPUT, symbol="%"),
        ]
    
    def _calculate(self, mA: pint.Quantity, LRV: pint.Quantity, URV: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        signal = mA.magnitude if hasattr(mA, 'magnitude') else float(mA)
        lrv = LRV.magnitude if hasattr(LRV, 'magnitude') else float(LRV)
        urv = URV.magnitude if hasattr(URV, 'magnitude') else float(URV)
        
        # 4mA = 0%, 20mA = 100%
        pct = (signal - 4) / 16 * 100
        pv = lrv + (urv - lrv) * pct / 100
        
        return {
            "PV": ureg.Quantity(pv, ""),
            "percent": ureg.Quantity(pct, "")
        }


class ControlValveCv(BaseEquation):
    """Control valve Cv sizing for liquid service."""
    
    equation_id = "valve_cv"
    name = "Control Valve Cv (Liquid)"
    category = "Instrumentation"
    description = "Size control valve Cv for liquid service"
    reference = "ISA 75.01 / IEC 60534"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Volumetric flow rate", "gal/min", "m^3/hr",
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("dP", "Pressure drop across valve", "psi", "bar",
                              ParameterType.INPUT, symbol="ΔP"),
            EquationParameter("SG", "Specific gravity (water = 1)", "dimensionless", "",
                              ParameterType.INPUT, symbol="SG", typical_range=(0.5, 2)),
            EquationParameter("Cv", "Valve flow coefficient", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Cv"),
        ]
    
    def _calculate(self, Q: pint.Quantity, dP: pint.Quantity, SG: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gal/min').magnitude
        dp = dP.to('psi').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        
        # ISA equation: Cv = Q * sqrt(SG/dP)
        cv = q * np.sqrt(sg / dp)
        
        return {"Cv": ureg.Quantity(cv, "")}


class SensorResponseTime(BaseEquation):
    """First-order sensor response time calculation."""
    
    equation_id = "sensor_response"
    name = "Sensor Response Time"
    category = "Instrumentation"
    description = "Calculate sensor response for step change (first-order)"
    reference = "Standard control theory"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("tau", "Time constant", "s", "",
                              ParameterType.INPUT, symbol="τ"),
            EquationParameter("t", "Elapsed time", "s", "",
                              ParameterType.INPUT, symbol="t"),
            EquationParameter("PV_final", "Final steady-state value", "dimensionless", "",
                              ParameterType.INPUT, symbol="PVf"),
            EquationParameter("PV_initial", "Initial value", "dimensionless", "",
                              ParameterType.INPUT, symbol="PVi"),
            EquationParameter("PV", "Process value at time t", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="PV(t)"),
            EquationParameter("percent_response", "Percent of final value", "dimensionless", "%",
                              ParameterType.OUTPUT, symbol="%"),
        ]
    
    def _calculate(self, tau: pint.Quantity, t: pint.Quantity, PV_final: pint.Quantity,
                   PV_initial: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        tc = tau.to('s').magnitude
        time = t.to('s').magnitude
        pvf = PV_final.magnitude if hasattr(PV_final, 'magnitude') else float(PV_final)
        pvi = PV_initial.magnitude if hasattr(PV_initial, 'magnitude') else float(PV_initial)
        
        # First-order response: PV(t) = PVi + (PVf - PVi) * (1 - exp(-t/tau))
        pv = pvi + (pvf - pvi) * (1 - np.exp(-time / tc))
        pct = (1 - np.exp(-time / tc)) * 100
        
        return {
            "PV": ureg.Quantity(pv, ""),
            "percent_response": ureg.Quantity(pct, "")
        }


# Update registry
INSTRUMENTATION_EQUATIONS = {
    'orifice_design': OrificeDesign,
    'thermowell_wake': ThermowellWakeFrequency,
    'transmitter_span': TransmitterSpan,
    'valve_noise': ControlValveNoise,
    'rtd_temperature': RTDTemperature,
    'signal_conversion': SignalConversion,
    'valve_cv': ControlValveCv,
    'sensor_response': SensorResponseTime,
}
