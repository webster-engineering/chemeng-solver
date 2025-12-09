"""
Heat Transfer equations - LMTD, NTU-effectiveness, heat transfer coefficients.
Each equation includes detailed variable descriptions and typical units.
"""

import numpy as np
from typing import Dict, List
import pint

from .base import BaseEquation, EquationParameter, ParameterType
from core.units import ureg
from core.learning import (
    LearningContent, QuizQuestion, QuestionType, Difficulty,
    WorkedExample, CalculationStep
)


# ===============================
# LMTD LEARNING CONTENT
# ===============================
LMTD_LEARNING = LearningContent(
    background_theory="""
The **Log Mean Temperature Difference (LMTD)** is the true driving force for heat transfer 
in heat exchangers. It accounts for the changing temperature difference along the exchanger.

**The Equation:** LMTD = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)

**For Counterflow:**
- ΔT₁ = T_hot,in - T_cold,out (hot end)
- ΔT₂ = T_hot,out - T_cold,in (cold end)

**Why Log Mean?**
Temperature difference varies exponentially along the exchanger length. The logarithmic 
mean correctly averages this for heat transfer calculations.

**The Design Equation:** Q = U × A × F × LMTD
Where F is the correction factor for multi-pass arrangements (F=1 for true counterflow).
""",
    key_concepts=[
        "LMTD is the effective average temperature difference",
        "Always larger at one end (approach temperature matters!)",
        "Counterflow gives highest LMTD for given terminal temperatures",
        "If ΔT₁ ≈ ΔT₂, then LMTD ≈ ΔT (arithmetic mean works)"
    ],
    real_world_applications=[
        "Shell and tube heat exchanger design",
        "Plate heat exchanger sizing",
        "Condenser and reboiler design",
        "Process heater calculations"
    ],
    common_mistakes=[
        "Using wrong flow arrangement (co-current vs counter-current)",
        "Temperature crossover (physically impossible - need larger area)",
        "Forgetting F factor for shell-side multi-pass",
        "Using LMTD when outlet temperatures are unknown (use NTU instead)"
    ],
    pro_tips=[
        "If designing, counterflow ALWAYS gives smaller area than co-current",
        "Target minimum approach temp of 10-20°F for process exchangers",
        "For close approach temps, check if NTU method is easier"
    ],
    quiz_questions=[
        QuizQuestion(
            id="lmtd_q1",
            question="If ΔT₁ = 50°F and ΔT₂ = 50°F (parallel approaches), what is LMTD?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["50°F", "25°F", "0°F", "Undefined"],
            correct_answer="50°F",
            explanation="When ΔT₁ = ΔT₂, the limit of the LMTD formula equals the arithmetic mean, which is simply ΔT.",
            difficulty=Difficulty.BEGINNER,
            points=10
        ),
        QuizQuestion(
            id="lmtd_q2",
            question="Hot in=200°F, Hot out=120°F, Cold in=80°F, Cold out=160°F (counterflow). What is LMTD?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="40",
            explanation="ΔT₁=200-160=40°F, ΔT₂=120-80=40°F. Since equal, LMTD=40°F.",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Calculate ΔT at each end first",
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="LMTD for Oil-Water Heat Exchanger",
        scenario="Size a counterflow HX to cool oil from 180°F to 120°F using cooling water from 75°F to 95°F.",
        given_values={
            "T_hot_in": "180°F (oil inlet)",
            "T_hot_out": "120°F (oil outlet)",
            "T_cold_in": "75°F (CW inlet)",
            "T_cold_out": "95°F (CW outlet)"
        },
        find=["LMTD for design"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate Temperature Approaches",
                description="Find ΔT at each end (counterflow arrangement)",
                formula="ΔT₁ = T_hot,in - T_cold,out, ΔT₂ = T_hot,out - T_cold,in",
                substitution="ΔT₁ = 180 - 95 = 85°F, ΔT₂ = 120 - 75 = 45°F",
                result="ΔT₁ = 85°F, ΔT₂ = 45°F"
            ),
            CalculationStep(
                step_number=2,
                title="Apply LMTD Formula",
                description="Calculate log mean",
                formula="LMTD = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)",
                substitution="LMTD = (85 - 45) / ln(85/45)",
                computation="LMTD = 40 / ln(1.889) = 40 / 0.636 = 62.9°F",
                result="LMTD = 62.9°F"
            )
        ],
        final_answer="LMTD = 62.9°F",
        real_world_context="Use this in Q = U·A·LMTD to find required exchanger area."
    ),
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=12,
    prerequisites=[],
    related_equations=["hx_area", "heat_duty", "ntu_eff"],
    references=[
        "Perry's Chemical Engineers' Handbook, 9th Ed., Section 11",
        "Kern, D.Q. 'Process Heat Transfer', McGraw-Hill",
        "TEMA (Tubular Exchanger Manufacturers Association) Standards"
    ],
    derivation_summary="LMTD comes from integrating the heat transfer rate along the exchanger length, assuming constant U. The log mean properly averages the varying temperature difference.",
    limitations_assumptions=[
        "Assumes constant overall heat transfer coefficient U",
        "Valid for single-pass or true counterflow/cocurrent",
        "Multi-pass requires correction factor F",
        "Temperature cross situations may require special analysis"
    ],
    diagram_type="heat_exchanger",
    diagram_labels={
        "shell": "Shell Side (Hot)",
        "tubes": "Tube Side (Cold)",
        "inlet_hot": "Hot Fluid In",
        "outlet_cold": "Cold Fluid Out"
    }
)


# ===============================
# HEAT DUTY LEARNING CONTENT
# ===============================
HEAT_DUTY_LEARNING = LearningContent(
    background_theory="""
**Heat Duty (Q)** is theA rate of heat transfer in a process. It's the starting point for 
all heat exchanger and heater/cooler designs.

**The Equation:** Q = ṁ × Cp × ΔT

Where:
- **Q** = Heat duty (BTU/hr, kW)
- **ṁ** = Mass flow rate (lb/hr, kg/s)
- **Cp** = Specific heat capacity (BTU/lb·°F, kJ/kg·K)
- **ΔT** = Temperature change (T_out - T_in)

**Energy Balance:**
For heat exchangers without heat loss:
Q_hot = Q_cold (heat lost by hot fluid = heat gained by cold fluid)

This is the fundamental energy balance that links both sides of the exchanger.
""",
    key_concepts=[
        "Heat duty is the total rate of energy transferred",
        "Q = ṁCpΔT is the sensible heat equation",
        "For phase change, add latent heat: Q = ṁ × h_fg",
        "Energy balance: Q_hot = Q_cold for adiabatic operation"
    ],
    real_world_applications=[
        "Sizing heaters and coolers",
        "Determining utility requirements (steam, cooling water)",
        "Energy balance on distillation columns",
        "Calculating condenser and reboiler duties"
    ],
    common_mistakes=[
        "Using Cp at wrong temperature (it varies!)",
        "Forgetting latent heat during phase change",
        "Sign errors in ΔT (heating vs cooling)",
        "Mixing unit systems (SI vs Imperial)"
    ],
    quiz_questions=[
        QuizQuestion(
            id="hd_q1",
            question="10,000 lb/hr of water cools from 150°F to 100°F. What's Q?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="500000",
            explanation="Q = ṁCpΔT = 10000 × 1.0 × 50 = 500,000 BTU/hr",
            difficulty=Difficulty.BEGINNER,
            hint="Water Cp = 1.0 BTU/(lb·°F)",
            points=10
        )
    ],
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=10,
    prerequisites=[],
    related_equations=["lmtd", "hx_area"],
    references=[
        "Perry's Chemical Engineers' Handbook, 9th Ed., Section 11",
        "McCabe, Smith & Harriott, Unit Operations of Chemical Engineering",
        "Incropera & DeWitt, Fundamentals of Heat and Mass Transfer"
    ],
    derivation_summary="Derived from the First Law of Thermodynamics applied to a flowing fluid. The rate of enthalpy change equals the heat transfer rate: Q = ṁΔH = ṁCpΔT for sensible heat.",
    limitations_assumptions=[
        "Assumes constant Cp over temperature range",
        "Does not include phase change - add latent heat separately",
        "Assumes no heat loss to surroundings",
        "For non-ideal mixing, use enthalpy values directly"
    ]
)


# ===============================
# NTU-EFFECTIVENESS LEARNING CONTENT
# ===============================
NTU_EFF_LEARNING = LearningContent(
    background_theory="""
The **NTU-Effectiveness (ε-NTU) method** is an alternative to LMTD for heat exchanger analysis.
It's especially useful when **outlet temperatures are unknown**.

**The Number of Transfer Units (NTU):** NTU = UA / C_min

Where C = ṁ·Cp (heat capacity rate)

**Effectiveness (ε):** ε = Q_actual / Q_max = Q / (C_min × (T_h,in - T_c,in))

**For Counterflow:**
ε = [1 - exp(-NTU(1-Cr))] / [1 - Cr·exp(-NTU(1-Cr))]  (Cr ≠ 1)

Where Cr = C_min / C_max (heat capacity ratio)

**Special Cases:**
- Cr = 0 (condenser/evaporator): ε = 1 - exp(-NTU)
- Cr = 1: ε = NTU / (1 + NTU)
""",
    key_concepts=[
        "NTU is a dimensionless measure of exchanger size",
        "ε tells you what fraction of maximum possible heat is transferred",
        "C_min determines maximum possible heat transfer",
        "Use ε-NTU when outlet temps are unknown (rating problem)"
    ],
    real_world_applications=[
        "Rating existing heat exchangers",
        "Performance prediction after fouling",
        "Optimizing exchanger networks",
        "Pinch analysis in process integration"
    ],
    common_mistakes=[
        "Confusing C_min and C_max assignment",
        "Using wrong ε-NTU formula for the exchanger type",
        "Forgetting Cr = C_min/C_max, not C_max/C_min"
    ],
    quiz_questions=[
        QuizQuestion(
            id="ntu_q1",
            question="What does effectiveness ε = 0.7 mean?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["70% of maximum possible heat is transferred", "70% efficiency", "70% of area is used", "UA = 0.7"],
            correct_answer="70% of maximum possible heat is transferred",
            explanation="ε = Q_actual/Q_max, so 70% of the thermodynamically maximum heat transfer is achieved.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["lmtd", "heat_duty"],
    related_equations=["lmtd", "hx_area", "overall_u"],
    references=[
        "Kays, W.M. & London, A.L. 'Compact Heat Exchangers', McGraw-Hill",
        "Incropera & DeWitt, Fundamentals of Heat and Mass Transfer",
        "Shah, R.K. & Sekulic, D.P. 'Fundamentals of Heat Exchanger Design'"
    ],
    derivation_summary="NTU is defined as UA/C_min. Effectiveness equations are derived by integrating the heat transfer equations along the exchanger. Different configurations yield different ε-NTU relationships.",
    limitations_assumptions=[
        "Assumes constant U along the exchanger",
        "Different formulas apply to each exchanger configuration",
        "Phase change (Cr→0) uses simplified exponential formula",
        "For rating problems (known UA), this method is preferred over LMTD"
    ]
)


class LMTD(BaseEquation):
    """Log Mean Temperature Difference calculation."""
    
    equation_id = "lmtd"
    name = "Log Mean Temperature Difference"
    category = "Heat Transfer"
    description = "Calculate LMTD for shell-and-tube or double-pipe heat exchanger design"
    reference = "Perry's Chemical Engineers' Handbook"
    learning_content = LMTD_LEARNING
    
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
    learning_content = HEAT_DUTY_LEARNING
    
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
    learning_content = NTU_EFF_LEARNING
    
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


# ===============================
# CONVECTION COEFFICIENT LEARNING
# ===============================
CONVECTION_LEARNING = LearningContent(
    background_theory="""
The **Dittus-Boelter correlation** estimates the convection heat transfer 
coefficient for turbulent flow inside tubes.

**The Equation:** Nu = 0.023 × Re^0.8 × Pr^n

Where:
- Nu = Nusselt number = h×D/k
- Re = Reynolds number (must be > 10,000)
- Pr = Prandtl number = Cp×μ/k
- n = 0.4 for heating, 0.3 for cooling

**Rearranged for h:** h = 0.023 × (k/D) × Re^0.8 × Pr^n

**Valid Range:**
- Re > 10,000 (turbulent flow)
- 0.6 < Pr < 160
- L/D > 10 (fully developed)
""",
    key_concepts=[
        "Nu = hD/k is dimensionless heat transfer",
        "Re^0.8 shows turbulence enhances transfer",
        "n = 0.4 (heating) or 0.3 (cooling)",
        "Only valid for turbulent flow in tubes"
    ],
    common_mistakes=[
        "Using for laminar flow (Re < 2300)",
        "Forgetting heating vs cooling exponent",
        "Not using bulk fluid properties"
    ],
    quiz_questions=[
        QuizQuestion(
            id="conv_q1",
            question="In Dittus-Boelter, what is the Prandtl exponent for heating?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["0.4", "0.3", "0.8", "1.0"],
            correct_answer="0.4",
            explanation="For heating, use Pr^0.4. For cooling, use Pr^0.3.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=["reynolds_number"],
    related_equations=["overall_u", "ntu_eff"],
    references=[
        "Dittus, F.W. & Boelter, L.M.K., University of California Publications in Engineering, 1930",
        "Incropera & DeWitt, Fundamentals of Heat and Mass Transfer",
        "Perry's Chemical Engineers' Handbook, 9th Ed., Section 5"
    ],
    derivation_summary="Empirical correlation from 1930 based on experimental data for turbulent flow in smooth tubes. The 0.8 exponent on Re reflects the enhancing effect of turbulent mixing on heat transfer.",
    limitations_assumptions=[
        "Valid only for turbulent flow (Re > 10,000)",
        "Limited Prandtl number range: 0.6 < Pr < 160",
        "Requires fully developed flow (L/D > 10)",
        "Smooth tubes only - modified correlations exist for rough surfaces"
    ]
)


class ConvectionCoefficient(BaseEquation):
    """Estimate convection coefficient using Dittus-Boelter correlation."""
    
    equation_id = "convection_coeff"
    name = "Convection Coefficient (Dittus-Boelter)"
    category = "Heat Transfer"
    description = "Estimate film coefficient for turbulent flow in tubes"
    reference = "Dittus-Boelter correlation (Re > 10000)"
    learning_content = CONVECTION_LEARNING
    
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
