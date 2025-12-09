"""
Distillation equations - column design, reflux, stages.
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
# FENSKE LEARNING CONTENT
# ===============================
FENSKE_LEARNING = LearningContent(
    background_theory="""
The **Fenske equation** calculates the minimum number of theoretical stages required for a 
binary distillation at **total reflux** (all condensate returned, no product withdrawal).

**The Equation:**
N_min = ln[(xD/(1-xD)) × ((1-xB)/xB)] / ln(α)

Or equivalently:
N_min = ln[(xD/xB) × ((1-xB)/(1-xD))] / ln(α)

Where:
- **N_min** = Minimum theoretical stages (including reboiler)
- **xD** = Light key mole fraction in distillate
- **xB** = Light key mole fraction in bottoms
- **α** = Average relative volatility = K_light / K_heavy

**Why Total Reflux?**
At total reflux, operating lines merge with the y=x diagonal on a McCabe-Thiele diagram,
giving the largest possible steps and therefore the fewest stages. This represents the
theoretical minimum - real columns need more stages.

**Key Insight:** Fenske shows that higher α (easier separation) means fewer stages needed.
""",
    key_concepts=[
        "N_min is achieved only at total reflux (infinite energy, no product)",
        "Higher relative volatility = fewer stages needed",
        "Include reboiler as a theoretical stage",
        "For multicomponent: use light and heavy key components"
    ],
    real_world_applications=[
        "Quick feasibility check for new column designs",
        "Benchmarking column performance during troubleshooting",
        "First step in shortcut column design (Fenske-Underwood-Gilliland)",
        "Estimating minimum equipment size for cost studies"
    ],
    industry_examples=[
        "Propane/butane splitter minimum stages estimation",
        "Crude unit atmospheric column scoping",
        "Benzene/toluene separation preliminary design"
    ],
    common_mistakes=[
        "Forgetting that N_min includes the reboiler as a stage",
        "Using arithmetic average for relative volatility instead of geometric mean",
        "Expecting to operate at N_min (impossible without total reflux)",
        "Confusing mole fractions of light key vs total composition"
    ],
    pro_tips=[
        "Use geometric mean for alpha: alpha_avg = sqrt(alpha_top * alpha_bottom)",
        "If alpha varies significantly, column may need more detailed simulation",
        "Rule of thumb: actual stages ~ 2.0 to 2.5 times N_min at typical reflux",
        "For multicomponent: use only the two key components in Fenske"
    ],
    variable_sources={
        "xD": "Distillate product specification. From customer requirements, downstream unit needs, or sales grade specifications. Typical high-purity products: 0.95-0.999.",
        "xB": "Bottoms product specification. Set by product value recovery, environmental limits, or downstream process needs. Lower xB = purer bottoms but more stages.",
        "alpha": "Calculate from VLE data or simulation at average column conditions. Use alpha = K_LK/K_HK. Geometric mean of top and bottom is typical: alpha = sqrt(alpha_top * alpha_bottom)."
    },
    quiz_questions=[
        QuizQuestion(
            id="fen_q1",
            question="What operating condition gives the minimum number of stages in the Fenske equation?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Total reflux", "Minimum reflux", "Zero reflux", "80% flood"],
            correct_answer="Total reflux",
            explanation="At total reflux, all condensate is returned to the column (no distillate product). This gives maximum driving force and minimum stages.",
            difficulty=Difficulty.BEGINNER,
            hint="What happens when R approaches infinity?",
            points=10
        ),
        QuizQuestion(
            id="fen_q2",
            question="If alpha = 2.5, xD = 0.98, xB = 0.02, calculate N_min.",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="8.5",
            explanation="N_min = ln[(0.98/0.02)(0.98/0.02)] / ln(2.5) = ln[(49)(49)] / ln(2.5) = ln(2401) / 0.916 = 7.78/0.916 = 8.5 stages",
            difficulty=Difficulty.INTERMEDIATE,
            hint="N_min = ln[(xD/(1-xD))((1-xB)/xB)] / ln(alpha)",
            points=15
        ),
        QuizQuestion(
            id="fen_q3",
            question="The Fenske equation applies to which mixture type?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Binary or pseudo-binary (with key components)", "Only pure binary", "Only azeotropic", "Only ideal mixtures"],
            correct_answer="Binary or pseudo-binary (with key components)",
            explanation="Fenske works for binary mixtures and can be applied to multicomponent mixtures by focusing on the light key and heavy key components.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="Depropanizer Minimum Stages",
        scenario="A depropanizer separates propane (LK) from butane (HK). The distillate should be 98% propane, and bottoms should contain only 2% propane. The average relative volatility is 2.5.",
        given_values={
            "xD": "0.98 (98% propane in distillate)",
            "xB": "0.02 (2% propane in bottoms)",
            "alpha": "2.5"
        },
        find=["Minimum theoretical stages at total reflux"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate separation factor",
                description="Compute the combined distillate/bottoms purity term.",
                formula="Separation = (xD/(1-xD)) × ((1-xB)/xB)",
                substitution="Separation = (0.98/0.02) × (0.98/0.02)",
                computation="Separation = 49 × 49 = 2401",
                result="Separation factor = 2401"
            ),
            CalculationStep(
                step_number=2,
                title="Apply Fenske equation",
                description="Take natural log and divide by ln(alpha).",
                formula="N_min = ln(Separation) / ln(alpha)",
                substitution="N_min = ln(2401) / ln(2.5)",
                computation="N_min = 7.78 / 0.916 = 8.49",
                result="N_min = 8.5 theoretical stages"
            ),
            CalculationStep(
                step_number=3,
                title="Interpret result",
                description="Understand what this means for design.",
                result="8.5 stages minimum (including reboiler)",
                notes="At real operating reflux, expect 15-20 theoretical stages. With 60% tray efficiency, actual trays would be 25-35."
            )
        ],
        final_answer="Minimum 8.5 theoretical stages at total reflux for 98% propane purity",
        real_world_context="A real depropanizer typically operates at R/R_min = 1.2-1.3, requiring about 20 theoretical stages or 30+ actual trays."
    ),
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=12,
    prerequisites=["raoults_law"],
    related_equations=["underwood", "gilliland", "mccabe_thiele"],
    references=[
        "Fenske, M.R. 'Fractionation of Straight-Run Pennsylvania Gasoline', Ind. Eng. Chem., 1932",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley, 4th Ed.",
        "Perry's Chemical Engineers' Handbook, Section 13"
    ],
    derivation_summary="Fenske derived the equation by applying material balances at total reflux, where the operating line coincides with y=x. This leads to a geometric series of composition ratios that simplifies to the logarithmic form.",
    limitations_assumptions=[
        "Applies only at total reflux (no products withdrawn)",
        "Assumes constant relative volatility throughout column",
        "Binary or pseudo-binary (key component) mixtures",
        "Theoretical stages - real trays need efficiency correction"
    ]
)


class FenskeMinStages(BaseEquation):
    """Fenske equation for minimum theoretical stages at total reflux."""
    
    equation_id = "fenske"
    name = "Fenske Minimum Stages"
    category = "Distillation"
    description = "Calculate minimum stages at total reflux using Fenske equation"
    reference = "Fenske, 1932"
    learning_content = FENSKE_LEARNING
    
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


# ===============================
# UNDERWOOD LEARNING CONTENT
# ===============================
UNDERWOOD_LEARNING = LearningContent(
    background_theory="""
The **Underwood equations** calculate the minimum reflux ratio (R_min) for a distillation 
column. At minimum reflux, an infinite number of stages would be required.

**The Two Underwood Equations:**

1. **First equation** (find theta):
   Sum over all components: (alpha_i * z_i) / (alpha_i - theta) = 1 - q

2. **Second equation** (find R_min):
   Sum over all components: (alpha_i * x_D,i) / (alpha_i - theta) = R_min + 1

Where:
- **theta** = Underwood root (lies between alpha values of light and heavy keys)
- **alpha_i** = Relative volatility of component i (to heavy key)
- **z_i** = Mole fraction of component i in feed
- **x_D,i** = Mole fraction of component i in distillate
- **q** = Feed thermal condition (1 = saturated liquid, 0 = saturated vapor)
- **R_min** = Minimum reflux ratio = L_min / D

**Feed Quality (q):**
- q = 1: Saturated liquid (all feed enters as liquid)
- q = 0: Saturated vapor (all feed enters as vapor)
- q > 1: Subcooled liquid
- 0 < q < 1: Two-phase feed
""",
    key_concepts=[
        "R_min corresponds to infinite stages (pinch point)",
        "Theta root must be found iteratively from first equation",
        "Feed quality q significantly affects R_min",
        "Saturated vapor feed (q=0) gives lower R_min than liquid feed"
    ],
    real_world_applications=[
        "Column design optimization (R vs stages tradeoff)",
        "Feed conditioning decisions (preheat to vaporize?)",
        "Revamp studies for existing columns",
        "Energy cost estimates for early-stage projects"
    ],
    industry_examples=[
        "Crude unit atmospheric column reflux sizing",
        "NGL fractionation train design",
        "Aromatics recovery column optimization"
    ],
    common_mistakes=[
        "Using wrong alpha basis (should be relative to heavy key)",
        "Forgetting that minimum reflux gives infinite stages",
        "Not iterating to find correct theta value",
        "Ignoring feed quality effects on R_min"
    ],
    pro_tips=[
        "Design at R = 1.1 to 1.5 times R_min for practical operation",
        "Higher R = fewer stages but more energy (condenser/reboiler duty)",
        "For binary: simplified formula R_min = (1/(alpha-1)) * (xD/zF - alpha*(1-xD)/(1-zF))",
        "If R_min < 0, the separation is trivially easy"
    ],
    variable_sources={
        "alpha": "Calculate from VLE K-values at average column conditions: alpha_i = K_i / K_HK. Use simulation or experimental data for accuracy.",
        "z_F": "Feed composition from laboratory analysis or process sampling. For design, use typical or worst-case values.",
        "x_D": "Distillate specification from product requirements. Sum of all x_D,i must equal 1.0.",
        "q": "Feed thermal condition. Determine from feed temperature vs bubble point. q = (H_V - H_F)/(H_V - H_L) where H is enthalpy."
    },
    quiz_questions=[
        QuizQuestion(
            id="und_q1",
            question="At minimum reflux, how many theoretical stages are required?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Infinite stages", "Minimum stages", "Zero stages", "One stage"],
            correct_answer="Infinite stages",
            explanation="Minimum reflux is the limiting case where a pinch occurs and infinite stages would be needed. It's the lower bound on reflux ratio.",
            difficulty=Difficulty.BEGINNER,
            hint="Think about what 'minimum' means in terms of stages",
            points=10
        ),
        QuizQuestion(
            id="und_q2",
            question="How does feed quality q affect minimum reflux?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Vapor feed (q=0) gives lower R_min than liquid (q=1)", 
                     "Liquid feed always gives lower R_min",
                     "Feed quality has no effect on R_min",
                     "Only subcooled feed affects R_min"],
            correct_answer="Vapor feed (q=0) gives lower R_min than liquid (q=1)",
            explanation="Vapor feed provides internal reflux directly, reducing external reflux needs. Liquid feed must be vaporized in the stripping section first.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["fenske", "raoults_law"],
    related_equations=["fenske", "gilliland", "mccabe_thiele"],
    references=[
        "Underwood, A.J.V. 'Fractional Distillation of Multicomponent Mixtures', Chem. Eng. Prog., 1948",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley, 4th Ed.",
        "Perry's Chemical Engineers' Handbook, Section 13"
    ],
    derivation_summary="Underwood derived the equations from material and energy balances at the pinch point (where operating and equilibrium lines touch). The theta root represents the composition at this pinch.",
    limitations_assumptions=[
        "Constant molar overflow (CMO) assumption",
        "Constant relative volatility",
        "Sharp split between key components (distributing components need modification)",
        "Single feed location"
    ]
)


class UnderwoodMinReflux(BaseEquation):
    """Underwood equation for minimum reflux ratio."""
    
    equation_id = "underwood"
    name = "Underwood Minimum Reflux"
    category = "Distillation"
    description = "Calculate minimum reflux ratio using Underwood equations"
    reference = "Underwood, 1948"
    learning_content = UNDERWOOD_LEARNING
    
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


# ===============================
# GILLILAND LEARNING CONTENT
# ===============================
GILLILAND_LEARNING = LearningContent(
    background_theory="""
The **Gilliland correlation** bridges the gap between Fenske (minimum stages) and Underwood 
(minimum reflux) by relating actual stages to actual reflux ratio.

**The Correlation:**
Y = (N - N_min) / (N + 1)  vs  X = (R - R_min) / (R + 1)

**Eduljee Fit (most common):**
Y = 0.75 × (1 - X^0.5668)

**The Physical Meaning:**
- At X = 0 (R = R_min): Y approaches 1, meaning N approaches infinity ✓
- At X = 1 (R = infinity): Y = 0, meaning N = N_min ✓

**Solving for Actual Stages:**
Given N_min (from Fenske), R_min (from Underwood), and chosen R:
1. Calculate X = (R - R_min) / (R + 1)
2. Calculate Y from correlation
3. Solve: N = (N_min + Y) / (1 - Y)

**Design Practice:**
- Choose R = 1.1 to 1.5 × R_min (typically 1.2-1.3)
- Lower R means more stages but less energy
- Higher R means fewer stages but more energy
""",
    key_concepts=[
        "Gilliland correlates ACTUAL stages vs ACTUAL reflux",
        "Connects the two extremes: total reflux and minimum reflux",
        "Empirical correlation - valid for most systems",
        "Actual trays = Theoretical stages / Tray efficiency"
    ],
    real_world_applications=[
        "Final step in Fenske-Underwood-Gilliland shortcut method",
        "Quick sizing for cost estimates and feasibility",
        "Trade-off analysis between stages and energy",
        "Revamp analysis - can column hit new spec?"
    ],
    industry_examples=[
        "Depropanizer column optimization study",
        "Debottlenecking analysis for crude unit",
        "Green field LNG fractionation preliminary design"
    ],
    common_mistakes=[
        "Using R < R_min (impossible operation)",
        "Forgetting to convert theoretical stages to actual trays",
        "Applying outside typical R/R_min range (1.05-2.0)",
        "Mixing up X and Y definitions"
    ],
    pro_tips=[
        "Economic optimum is typically R/R_min = 1.2 to 1.3",
        "At R/R_min = 1.2, expect N/N_min = 2.0 to 2.5",
        "Tray efficiency typically 50-80% for sieve/valve trays",
        "Add extra trays for design margin (10-20%)"
    ],
    variable_sources={
        "N_min": "Calculate from Fenske equation using alpha and xD, xB specifications. This is the minimum stages at total reflux.",
        "R_min": "Calculate from Underwood equations using alpha, feed composition, and feed quality. This is minimum reflux for infinite stages.",
        "R": "Design choice. Start at 1.2 x R_min and optimize based on energy vs capital cost tradeoff. Higher R = more operating cost but less capital."
    },
    quiz_questions=[
        QuizQuestion(
            id="gil_q1",
            question="At what R/R_min ratio do most columns operate?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["1.2 to 1.3", "0.8 to 0.9", "2.0 to 3.0", "Exactly 1.0"],
            correct_answer="1.2 to 1.3",
            explanation="R/R_min = 1.2-1.3 typically gives the best economic balance between capital cost (fewer stages at higher R) and operating cost (less energy at lower R).",
            difficulty=Difficulty.BEGINNER,
            points=10
        ),
        QuizQuestion(
            id="gil_q2",
            question="If N_min = 10, R_min = 2.0, and R = 3.0, what is X?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="0.25",
            explanation="X = (R - R_min)/(R + 1) = (3.0 - 2.0)/(3.0 + 1) = 1.0/4.0 = 0.25",
            difficulty=Difficulty.INTERMEDIATE,
            hint="X = (R - R_min)/(R + 1)",
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=["fenske", "underwood"],
    related_equations=["fenske", "underwood", "kirkbride"],
    references=[
        "Gilliland, E.R. 'Multicomponent Rectification', Ind. Eng. Chem., 1940",
        "Eduljee, H.E. 'Equations Replace Gilliland Plot', Hydrocarbon Processing, 1975",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley"
    ],
    derivation_summary="Gilliland plotted experimental data for actual N vs actual R. Eduljee later fit this data to the commonly used equation Y = 0.75(1 - X^0.5668).",
    limitations_assumptions=[
        "Empirical correlation (not rigorous)",
        "Best for R/R_min between 1.05 and 2.0",
        "Theoretical stages - multiply by efficiency for actual trays",
        "Assumes feed at optimal location"
    ]
)


class GillilandCorrelation(BaseEquation):
    """Gilliland correlation for actual stages."""
    
    equation_id = "gilliland"
    name = "Gilliland Correlation"
    category = "Distillation"
    description = "Estimate actual stages from N_min and R_min"
    reference = "Gilliland, 1940"
    learning_content = GILLILAND_LEARNING
    
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


# ============== Kirkbride Feed Stage Learning Content ==============
KIRKBRIDE_LEARNING = LearningContent(
    background_theory="""
## Kirkbride Equation for Feed Stage Location

The **Kirkbride equation** provides a quick estimate of the optimal feed stage location 
in a distillation column. Proper feed location is critical for column efficiency - 
feeding at the wrong stage wastes energy and reduces separation capability.

### The Equation

The ratio of rectifying to stripping stages is:

$$\\frac{N_R}{N_S} = \\left[ \\frac{B}{D} \\cdot \\frac{x_{F,HK}}{x_{F,LK}} \\cdot \\left(\\frac{x_{B,LK}}{x_{D,HK}}\\right)^2 \\right]^{0.206}$$

### Physical Significance

- **Matching compositions**: The optimal feed stage is where the internal column 
  liquid composition most closely matches the feed composition
- **Minimum remixing**: Proper feed location minimizes energy waste from remixing
- **Key components**: The equation uses heavy key (HK) and light key (LK) distributions
""",
    key_concepts=[
        "Feed should enter where internal compositions match feed composition",
        "Rectifying section is above the feed, stripping section is below",
        "The 0.206 exponent is empirically derived from industrial data",
        "B/D ratio affects the section split - higher bottoms flow pushes feed up",
        "Key component purity affects the balance between sections"
    ],
    real_world_applications=[
        "Designing new distillation column feed locations",
        "Troubleshooting existing columns with suboptimal performance",
        "Debottlenecking studies when changing feed compositions",
        "Rating existing column designs for different operating conditions"
    ],
    industry_examples=[
        "Crude distillation: multiple feed and side draw locations",
        "Ethylene splitter optimization for varying feed compositions",
        "Aromatics separation with BTX fractionators",
        "Depropanizer feed point adjustment for seasonal propane specs"
    ],
    common_mistakes=[
        "Forgetting the 0.206 exponent is on the entire bracketed term",
        "Confusing N_R (rectifying) with upper section numbering conventions",
        "Not accounting for feed quality (q) effects on optimal location",
        "Using molar compositions when weight fractions are given"
    ],
    pro_tips=[
        "Start with Kirkbride for initial estimate, refine with simulation",
        "Feed location sensitivity is highest for sharp separations (high purity)",
        "Two-stage adjustment: count from bottom for reboiler, from top for condenser",
        "For multicomponent systems, use the key components that define the split"
    ],
    variable_sources={
        "N": "Total theoretical stages from Fenske-Gilliland shortcut sequence. This is the sum of N_R + N_S. Typically 10-100+ for most separations.",
        "B_D": "Bottoms-to-distillate molar flow ratio calculated from overall mass balance: B/D = (xD - xF)/(xF - xB). Range 0.1-10 depending on cut point.",
        "x_F_HK": "Mole fraction of heavy key component in the feed. Obtained from feed composition analysis (GC, lab sampling). Typical LK+HK sum is 0.90+ for key-dominant feeds.",
        "x_F_LK": "Mole fraction of light key component in the feed. From same feed analysis as x_F_HK. Often paired with HK as the dominant components.",
        "x_D_HK": "Heavy key mole fraction in distillate spec. Set by product specification - typically 0.001-0.05 for overhead HK impurity limit.",
        "x_B_LK": "Light key mole fraction in bottoms spec. Set by bottoms product specification - typically 0.001-0.05 for LK loss limit."
    },
    quiz_questions=[
        QuizQuestion(
            id="kirk_q1",
            question="If B/D = 1, x_F,HK/x_F,LK = 1, and x_B,LK/x_D,HK = 1, what is N_R/N_S?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["0.5", "1.0", "2.0", "Undefined"],
            correct_answer="1.0",
            explanation="With all ratios = 1, (1 * 1 * 1)^0.206 = 1.0, so N_R = N_S (feed at midpoint).",
            difficulty=Difficulty.BEGINNER,
            hint="All bracketed terms equal 1",
            points=10
        ),
        QuizQuestion(
            id="kirk_q2",
            question="A column has 30 total stages and Kirkbride gives N_R/N_S = 0.5. What is the optimal feed stage from the bottom?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="20",
            explanation="N_R/N_S = 0.5 means N_R = 0.5*N_S. Since N_R + N_S = 30: 0.5*N_S + N_S = 30, so N_S = 20. Feed stage from bottom = N_S = 20.",
            difficulty=Difficulty.INTERMEDIATE,
            hint="N_S is the stripping section (below feed)",
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=10,
    prerequisites=["fenske", "gilliland"],
    related_equations=["fenske", "underwood", "gilliland"],
    references=[
        "Kirkbride, C.G. 'Process Design Procedure for Multicomponent Fractionators', Petroleum Refiner, 1944",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley",
        "Ludwig, 'Applied Process Design for Chemical and Petrochemical Plants'"
    ],
    derivation_summary="Kirkbride empirically correlated rectifying/stripping stage ratios with flow ratios and key component distributions from industrial column data.",
    limitations_assumptions=[
        "Empirical correlation - use simulation for rigorous design",
        "Assumes binary-like behavior with key components",
        "Does not account for feed thermal condition (q)",
        "Best accuracy for columns with 10-80 stages"
    ]
)

# ============== Tray Hydraulics Learning Content ==============
TRAY_HYDRAULICS_LEARNING = LearningContent(
    background_theory="""
## Tray Hydraulics and Flooding

**Tray hydraulics** determines the operational limits of distillation columns. 
The most critical limit is **flooding**, which occurs when vapor velocity is so high 
that liquid cannot flow down the column properly.

### Flooding Mechanism

At flooding:
1. Liquid accumulates on trays (liquid holdup increases dramatically)
2. Pressure drop spikes
3. Separation efficiency collapses
4. Column may become unstable

### Fair Correlation

The flooding velocity is calculated from:

$$U_{flood} = C_{sb} \\left(\\frac{\\sigma}{20}\\right)^{0.2} \\sqrt{\\frac{\\rho_L - \\rho_V}{\\rho_V}}$$

Where $C_{sb}$ is the Souders-Brown capacity factor (from correlations or charts).

### Design Approach

Columns typically operate at **70-85% of flooding** to provide operating margin 
while maintaining good efficiency.
""",
    key_concepts=[
        "Flooding is the upper vapor velocity limit for stable operation",
        "Weeping is the lower liquid flow limit where liquid falls through holes",
        "Surface tension affects bubble/froth behavior (sigma^0.2 term)",
        "Density difference drives liquid-vapor separation on trays",
        "Tray spacing affects maximum capacity - taller spacing allows higher velocities"
    ],
    real_world_applications=[
        "Rating existing columns for capacity increases",
        "Designing new columns for required throughput",
        "Troubleshooting flooding or weeping problems",
        "Debottlenecking studies for capacity expansion"
    ],
    industry_examples=[
        "Crude unit atmospheric column capacity checks",
        "FCC main fractionator flood assessment",
        "Amine regenerator tray hydraulic evaluation",
        "Vacuum column flood checks at low absolute pressures"
    ],
    common_mistakes=[
        "Using volumetric vapor flow instead of velocity",
        "Forgetting to convert surface tension units (dyne/cm is standard)",
        "Not accounting for foaming tendency which reduces effective capacity",
        "Ignoring tray geometry effects on Csb"
    ],
    pro_tips=[
        "70-75% flood for clean services, 65-70% for fouling or foaming services",
        "Check flood at both ends of column - vapor rate varies with reflux",
        "Low pressure systems have higher vapor volumes - check flooding carefully",
        "Structured packing has different flooding correlations than trays"
    ],
    variable_sources={
        "V": "Vapor molar flow rate from column simulation or heat/mass balance. Typically from process simulator or hand calculation. Units often kmol/hr or lbmol/hr.",
        "MW_v": "Vapor molecular weight calculated as mole-fraction weighted average of component MWs. For mixtures: sum(yi * MWi). Typical range 20-150 for organic systems.",
        "rho_v": "Vapor density from ideal gas law or equation of state at column T and P. PV = nRT rearranged gives rho_v = P*MW/(R*T). Typical 0.5-5 kg/m3 for atmospheric columns.",
        "rho_l": "Liquid density from pure component data or mixing rules at column temperature. Typically 500-900 kg/m3 for organics, 1000 kg/m3 for water.",
        "sigma": "Surface tension from pure component data or estimation (Parachor method). Water is 72 dyne/cm, most organics 15-35 dyne/cm. Use dyne/cm units.",
        "A_net": "Net tray area = total tower cross-section minus downcomer area. A_net = (pi*D^2/4) - A_downcomer. Typically 85-92% of total area."
    },
    quiz_questions=[
        QuizQuestion(
            id="tray_q1",
            question="If vapor density doubles (rho_v * 2) while liquid density stays constant, what happens to flooding velocity?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Increases by sqrt(2)", "Decreases by sqrt(2)", "Stays the same", "Doubles"],
            correct_answer="Decreases by sqrt(2)",
            explanation="U_flood is proportional to sqrt((rho_L - rho_V)/rho_V). Higher rho_V decreases numerator and increases denominator, lowering U_flood.",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Consider both numerator and denominator effects",
            points=15
        ),
        QuizQuestion(
            id="tray_q2",
            question="A column operates at 2.5 ft/s actual velocity with flooding velocity of 3.5 ft/s. What is the percent flood?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="71.4",
            explanation="%Flood = (U_actual / U_flood) * 100 = (2.5 / 3.5) * 100 = 71.4%",
            difficulty=Difficulty.BEGINNER,
            hint="Percent flood = (actual/flooding) * 100",
            points=10
        )
    ],
    difficulty=Difficulty.ADVANCED,
    estimated_time_minutes=15,
    prerequisites=["gilliland"],
    related_equations=["packed_column_dp"],
    references=[
        "Fair, J.R. 'How to Design Packed Columns', Chem. Eng., 1961",
        "Kister, H.Z. 'Distillation Design', McGraw-Hill",
        "Perry's Chemical Engineers' Handbook, Section 14"
    ],
    derivation_summary="The Fair/Souders-Brown correlation is empirically derived from flooding data on commercial trays. Csb values are tabulated or correlated with flow parameter FLV.",
    limitations_assumptions=[
        "Fair correlation - other correlations exist (Glitsch, Koch, etc.)",
        "Assumes non-foaming system (apply derating for foaming)",
        "Does not include weeping check (lower operating limit)",
        "Csb varies with tray spacing - correlation assumes 24-inch spacing"
    ]
)


class FeedStageLocation(BaseEquation):
    """Kirkbride equation for optimal feed stage location."""
    
    equation_id = "kirkbride"
    name = "Feed Stage Location (Kirkbride)"
    category = "Distillation"
    description = "Estimate optimal feed tray location"
    reference = "Kirkbride, 1944"
    learning_content = KIRKBRIDE_LEARNING
    
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
    learning_content = TRAY_HYDRAULICS_LEARNING
    
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
