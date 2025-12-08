"""
Basic Math equations - fundamental calculations for engineering.
Includes conversions, geometry, statistics, and common formulas.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg


class QuadraticFormula(BaseEquation):
    """Solve quadratic equation ax² + bx + c = 0."""
    
    equation_id = "quadratic"
    name = "Quadratic Formula"
    category = "Basic Math"
    description = "Solve ax² + bx + c = 0 for x"
    reference = "Standard algebra"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("a", "Coefficient of x² term", 
                              "dimensionless", "", ParameterType.INPUT, symbol="a"),
            EquationParameter("b", "Coefficient of x term", 
                              "dimensionless", "", ParameterType.INPUT, symbol="b"),
            EquationParameter("c", "Constant term", 
                              "dimensionless", "", ParameterType.INPUT, symbol="c"),
            EquationParameter("x1", "First root", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="x₁"),
            EquationParameter("x2", "Second root", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="x₂"),
            EquationParameter("discriminant", "Discriminant (b²-4ac) - indicates nature of roots", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Δ",
                              tooltip=">0: 2 real roots, =0: 1 real root, <0: complex roots"),
        ]
    
    def _calculate(self, a: pint.Quantity, b: pint.Quantity, c: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        a_val = a.magnitude if hasattr(a, 'magnitude') else float(a)
        b_val = b.magnitude if hasattr(b, 'magnitude') else float(b)
        c_val = c.magnitude if hasattr(c, 'magnitude') else float(c)
        
        discriminant = b_val**2 - 4*a_val*c_val
        
        if discriminant >= 0:
            x1 = (-b_val + np.sqrt(discriminant)) / (2*a_val)
            x2 = (-b_val - np.sqrt(discriminant)) / (2*a_val)
        else:
            x1 = float('nan')
            x2 = float('nan')
        
        return {
            "x1": ureg.Quantity(x1, ""),
            "x2": ureg.Quantity(x2, ""),
            "discriminant": ureg.Quantity(discriminant, "")
        }


class LinearInterpolation(BaseEquation):
    """Linear interpolation between two points."""
    
    equation_id = "linear_interp"
    name = "Linear Interpolation"
    category = "Basic Math"
    description = "Find y at x between two known points (x1,y1) and (x2,y2)"
    reference = "Standard mathematics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x1", "First x value (lower bound)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="x₁"),
            EquationParameter("y1", "First y value (at x1)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="y₁"),
            EquationParameter("x2", "Second x value (upper bound)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="x₂"),
            EquationParameter("y2", "Second y value (at x2)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="y₂"),
            EquationParameter("x", "Target x value to interpolate", 
                              "dimensionless", "", ParameterType.INPUT, symbol="x"),
            EquationParameter("y", "Interpolated y value", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y"),
            EquationParameter("slope", "Line slope (Δy/Δx)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="m"),
        ]
    
    def _calculate(self, x1: pint.Quantity, y1: pint.Quantity, x2: pint.Quantity,
                   y2: pint.Quantity, x: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        x1_val = x1.magnitude if hasattr(x1, 'magnitude') else float(x1)
        y1_val = y1.magnitude if hasattr(y1, 'magnitude') else float(y1)
        x2_val = x2.magnitude if hasattr(x2, 'magnitude') else float(x2)
        y2_val = y2.magnitude if hasattr(y2, 'magnitude') else float(y2)
        x_val = x.magnitude if hasattr(x, 'magnitude') else float(x)
        
        slope = (y2_val - y1_val) / (x2_val - x1_val)
        y = y1_val + slope * (x_val - x1_val)
        
        return {
            "y": ureg.Quantity(y, ""),
            "slope": ureg.Quantity(slope, "")
        }


class PercentChange(BaseEquation):
    """Calculate percent change between two values."""
    
    equation_id = "percent_change"
    name = "Percent Change"
    category = "Basic Math"
    description = "Calculate % change from old value to new value"
    reference = "Standard mathematics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("old_value", "Original/initial value", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Old"),
            EquationParameter("new_value", "New/final value", 
                              "dimensionless", "", ParameterType.INPUT, symbol="New"),
            EquationParameter("percent_change", "Percent change (can be negative)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Δ%"),
            EquationParameter("ratio", "Ratio of new to old", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="New/Old"),
        ]
    
    def _calculate(self, old_value: pint.Quantity, new_value: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        old = old_value.magnitude if hasattr(old_value, 'magnitude') else float(old_value)
        new = new_value.magnitude if hasattr(new_value, 'magnitude') else float(new_value)
        
        pct_change = (new - old) / old * 100
        ratio = new / old
        
        return {
            "percent_change": ureg.Quantity(pct_change, ""),
            "ratio": ureg.Quantity(ratio, "")
        }


class MeanAndStdDev(BaseEquation):
    """Calculate mean and standard deviation from up to 10 values."""
    
    equation_id = "statistics"
    name = "Mean and Standard Deviation"
    category = "Basic Math"
    description = "Calculate mean, std dev, min, max from values (enter 0s to skip unused inputs)"
    reference = "Standard statistics"
    
    def get_parameters(self) -> List[EquationParameter]:
        params = []
        for i in range(1, 11):
            params.append(EquationParameter(
                f"x{i}", f"Value {i} (enter 0 to skip)", 
                "dimensionless", "", ParameterType.INPUT, symbol=f"x{i}",
                required=(i <= 2)
            ))
        params.extend([
            EquationParameter("mean", "Arithmetic mean", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="μ"),
            EquationParameter("std_dev", "Sample standard deviation", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="σ"),
            EquationParameter("min_val", "Minimum value", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Min"),
            EquationParameter("max_val", "Maximum value", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Max"),
            EquationParameter("range_val", "Range (max - min)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Range"),
            EquationParameter("n", "Count of non-zero values", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="n"),
        ])
        return params
    
    def _calculate(self, **kwargs) -> Dict[str, pint.Quantity]:
        values = []
        for i in range(1, 11):
            key = f"x{i}"
            if key in kwargs and kwargs[key] is not None:
                val = kwargs[key]
                v = val.magnitude if hasattr(val, 'magnitude') else float(val)
                if v != 0:
                    values.append(v)
        
        if len(values) < 2:
            values = [0, 0]
        
        mean = np.mean(values)
        std = np.std(values, ddof=1)
        min_v = np.min(values)
        max_v = np.max(values)
        
        return {
            "mean": ureg.Quantity(mean, ""),
            "std_dev": ureg.Quantity(std, ""),
            "min_val": ureg.Quantity(min_v, ""),
            "max_val": ureg.Quantity(max_v, ""),
            "range_val": ureg.Quantity(max_v - min_v, ""),
            "n": ureg.Quantity(len(values), "")
        }


class CircleCalculations(BaseEquation):
    """Calculate circle properties from radius or diameter."""
    
    equation_id = "circle"
    name = "Circle Properties"
    category = "Basic Math"
    description = "Calculate area, circumference from radius"
    reference = "Geometry"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("r", "Radius of circle", 
                              "length", "in", ParameterType.INPUT, symbol="r"),
            EquationParameter("area", "Area = πr²", 
                              "area", "in**2", ParameterType.OUTPUT, symbol="A"),
            EquationParameter("circumference", "Circumference = 2πr", 
                              "length", "in", ParameterType.OUTPUT, symbol="C"),
            EquationParameter("diameter", "Diameter = 2r", 
                              "length", "in", ParameterType.OUTPUT, symbol="D"),
        ]
    
    def _calculate(self, r: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        r_val = r.to('in').magnitude
        
        area = np.pi * r_val**2
        circumference = 2 * np.pi * r_val
        diameter = 2 * r_val
        
        return {
            "area": ureg.Quantity(area, "in**2"),
            "circumference": ureg.Quantity(circumference, "in"),
            "diameter": ureg.Quantity(diameter, "in")
        }


class CylinderVolume(BaseEquation):
    """Calculate cylinder volume and surface area."""
    
    equation_id = "cylinder"
    name = "Cylinder Volume"
    category = "Basic Math"
    description = "Calculate volume and surface area of cylinder"
    reference = "Geometry"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("r", "Radius of cylinder", 
                              "length", "in", ParameterType.INPUT, symbol="r"),
            EquationParameter("h", "Height of cylinder", 
                              "length", "in", ParameterType.INPUT, symbol="h"),
            EquationParameter("volume", "Volume = πr²h", 
                              "volume", "in**3", ParameterType.OUTPUT, symbol="V"),
            EquationParameter("volume_gal", "Volume in gallons", 
                              "volume", "gal", ParameterType.OUTPUT, symbol="V (gal)"),
            EquationParameter("surface_area", "Total surface area", 
                              "area", "in**2", ParameterType.OUTPUT, symbol="A"),
            EquationParameter("lateral_area", "Lateral (side) surface area", 
                              "area", "in**2", ParameterType.OUTPUT, symbol="Aₗₐₜ"),
        ]
    
    def _calculate(self, r: pint.Quantity, h: pint.Quantity, 
                   **kwargs) -> Dict[str, pint.Quantity]:
        r_val = r.to('in').magnitude
        h_val = h.to('in').magnitude
        
        volume = np.pi * r_val**2 * h_val
        volume_gal = volume / 231  # 231 in³ per gallon
        lateral = 2 * np.pi * r_val * h_val
        surface = lateral + 2 * np.pi * r_val**2
        
        return {
            "volume": ureg.Quantity(volume, "in**3"),
            "volume_gal": ureg.Quantity(volume_gal, "gal"),
            "surface_area": ureg.Quantity(surface, "in**2"),
            "lateral_area": ureg.Quantity(lateral, "in**2")
        }




class PythagoreanTheorem(BaseEquation):
    """Calculate hypotenuse or leg of right triangle."""
    
    equation_id = "pythagorean"
    name = "Pythagorean Theorem"
    category = "Basic Math"
    description = "a² + b² = c² - solve for any side of right triangle"
    reference = "Geometry"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("a", "First leg (side a)", 
                              "length", "ft", ParameterType.INPUT, symbol="a"),
            EquationParameter("b", "Second leg (side b)", 
                              "length", "ft", ParameterType.INPUT, symbol="b"),
            EquationParameter("c", "Hypotenuse = √(a² + b²)", 
                              "length", "ft", ParameterType.OUTPUT, symbol="c"),
            EquationParameter("angle_A", "Angle opposite side a", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="∠A (deg)"),
            EquationParameter("angle_B", "Angle opposite side b", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="∠B (deg)"),
        ]
    
    def _calculate(self, a: pint.Quantity, b: pint.Quantity, 
                   **kwargs) -> Dict[str, pint.Quantity]:
        a_val = a.to('ft').magnitude
        b_val = b.to('ft').magnitude
        
        c = np.sqrt(a_val**2 + b_val**2)
        angle_a = np.degrees(np.arctan(a_val / b_val)) if b_val != 0 else 90
        angle_b = 90 - angle_a
        
        return {
            "c": ureg.Quantity(c, "ft"),
            "angle_A": ureg.Quantity(angle_a, ""),
            "angle_B": ureg.Quantity(angle_b, "")
        }


class ExponentialDecay(BaseEquation):
    """Calculate exponential decay or growth."""
    
    equation_id = "exponential"
    name = "Exponential Decay/Growth"
    category = "Basic Math"
    description = "Calculate y = y₀ × e^(kt) for decay (k<0) or growth (k>0)"
    reference = "Standard mathematics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("y0", "Initial value at t=0", 
                              "dimensionless", "", ParameterType.INPUT, symbol="y₀"),
            EquationParameter("k", "Rate constant (positive=growth, negative=decay)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="k",
                              tooltip="Also related to half-life: k = -ln(2)/t_half"),
            EquationParameter("t", "Time elapsed", 
                              "time", "min", ParameterType.INPUT, symbol="t"),
            EquationParameter("y", "Value at time t", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y"),
            EquationParameter("half_life", "Half-life (if k<0)", 
                              "time", "min", ParameterType.OUTPUT, symbol="t½"),
            EquationParameter("doubling_time", "Doubling time (if k>0)", 
                              "time", "min", ParameterType.OUTPUT, symbol="t₂ₓ"),
        ]
    
    def _calculate(self, y0: pint.Quantity, k: pint.Quantity, t: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        y0_val = y0.magnitude if hasattr(y0, 'magnitude') else float(y0)
        k_val = k.magnitude if hasattr(k, 'magnitude') else float(k)
        t_val = t.to('min').magnitude
        
        y = y0_val * np.exp(k_val * t_val)
        
        if k_val < 0:
            half_life = -np.log(2) / k_val
            doubling = float('inf')
        elif k_val > 0:
            half_life = float('inf')
            doubling = np.log(2) / k_val
        else:
            half_life = float('inf')
            doubling = float('inf')
        
        return {
            "y": ureg.Quantity(y, ""),
            "half_life": ureg.Quantity(half_life, "min"),
            "doubling_time": ureg.Quantity(doubling, "min")
        }


BASIC_MATH_EQUATIONS = {
    'quadratic': QuadraticFormula,
    'linear_interp': LinearInterpolation,
    'percent_change': PercentChange,
    'statistics': MeanAndStdDev,
    'circle': CircleCalculations,
    'cylinder': CylinderVolume,
    'pythagorean': PythagoreanTheorem,
    'exponential': ExponentialDecay,
}
