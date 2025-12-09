"""
Mass Transfer equations - absorption, distillation, extraction.
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
# MCCABE-THIELE LEARNING CONTENT
# ===============================
MCCABE_THIELE_LEARNING = LearningContent(
    background_theory="""
The **McCabe-Thiele method** is a graphical technique for determining the number of theoretical 
stages required in a binary distillation column. It's one of the most important methods taught 
in chemical engineering for understanding separation processes.

**Key Concept:** The method plots the vapor composition (y) vs liquid composition (x) and 
steps off between the equilibrium curve and operating lines to count stages.

**The Three Lines:**
1. **Equilibrium Curve:** y = αx / [1 + (α-1)x] where α = relative volatility
2. **Rectifying Operating Line:** y = (R/(R+1))x + xD/(R+1)
3. **Stripping Operating Line:** Connects feed to bottoms composition

**Important Variables:**
- **xD** = Distillate composition (light key mole fraction)
- **xB** = Bottoms composition (light key mole fraction)
- **xF** = Feed composition (light key mole fraction)
- **R** = Reflux ratio (L/D = liquid returned / distillate withdrawn)
- **α** = Relative volatility (K_light / K_heavy)

The method assumes constant molar overflow (CMO), meaning liquid and vapor rates remain 
constant in each section of the column.
""",
    key_concepts=[
        "Equilibrium curve from VLE data or relative volatility",
        "Rectifying section operates above feed",
        "Stripping section operates below feed",
        "More stages needed at lower reflux ratio",
        "Minimum reflux = infinite stages, Total reflux = minimum stages"
    ],
    real_world_applications=[
        "Crude oil distillation column design",
        "Ethanol-water separation in biofuel production",
        "BTX (benzene-toluene-xylene) separation in refineries",
        "Solvent recovery in pharmaceutical manufacturing"
    ],
    industry_examples=[
        "Sizing a methanol-water column in a chemical plant",
        "Designing depropanizer columns in natural gas processing",
        "Optimizing reflux ratio for energy vs capital tradeoff"
    ],
    common_mistakes=[
        "Using mole fraction > 1 (check: xD + xB compositions must be valid)",
        "Forgetting that α may vary with composition (average α)",
        "Not checking if R > R_min (cannot operate below minimum reflux)",
        "Confusing theoretical stages with actual trays (need efficiency)"
    ],
    pro_tips=[
        "Plot the diagram to visualize - even approximate sketches help",
        "α closer to 1.0 = harder separation = more stages needed",
        "Typical rule: Use R = 1.2 to 1.5 × R_min for economic optimum",
        "Add a partial reboiler and condenser as stages (+1 to +2)"
    ],
    variable_sources={
        "xD": "Product specification from downstream process requirements. Set by customer spec, sales grade, or next unit feed requirement. Typical values 0.90-0.99 for high-purity products.",
        "xB": "Bottoms spec from product purity or waste disposal limits. Environmental permits, economic value recovery, or residue quality requirements set this value.",
        "xF": "Feed analysis from laboratory GC/LC analysis or online process analyzer. Sample upstream of column and verify with multiple measurements.",
        "R": "Design choice - start at 1.2-1.5 times minimum reflux. Higher R = fewer stages but more energy; lower R = more stages but less energy. Economic optimization required.",
        "α": "Calculate from VLE data: α = K_light/K_heavy = (y/x)/(y'/x') at column conditions. Use average of top and bottom, or geometric mean. Check experiment or simulation."
    },
    quiz_questions=[
        QuizQuestion(
            id="mt_q1",
            question="What happens to the number of stages if you increase the reflux ratio?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Number of stages decreases", "Number of stages increases", 
                     "Stages stay the same", "Column diameter increases"],
            correct_answer="Number of stages decreases",
            explanation="Higher reflux brings operating line closer to equilibrium curve, allowing larger steps. At total reflux (R=∞), you get minimum stages.",
            difficulty=Difficulty.BEGINNER,
            hint="Think about what happens at total reflux vs minimum reflux",
            points=10
        ),
        QuizQuestion(
            id="mt_q2",
            question="If α = 2.0, xD = 0.95, xB = 0.05, what is the minimum number of stages (Fenske)?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="8.5",
            explanation="N_min = ln[(0.95/0.05)(0.95/0.05)] / ln(2.0) = ln(361) / ln(2) = 5.89/0.693 ≈ 8.5 stages",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Use Fenske equation: N_min = ln[(xD/(1-xD))((1-xB)/xB)] / ln(α)",
            points=15
        ),
        QuizQuestion(
            id="mt_q3",
            question="The McCabe-Thiele method assumes constant molar overflow. This is valid when:",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Heats of vaporization of components are similar", 
                     "The feed is always saturated liquid",
                     "The column operates at high pressure",
                     "The relative volatility is very high"],
            correct_answer="Heats of vaporization of components are similar",
            explanation="Constant molar overflow (CMO) assumes equal molar heats of vaporization so that when one mole condenses, one mole evaporates.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="Ethanol-Water Distillation Column",
        scenario="Design a distillation column to produce 95% ethanol from a 40% feed. The bottoms should contain no more than 5% ethanol. Use a reflux ratio of 2.0.",
        given_values={
            "xD": "0.95 (95% ethanol in distillate)",
            "xB": "0.05 (5% ethanol in bottoms)",
            "xF": "0.40 (40% ethanol in feed)",
            "R": "2.0 (reflux ratio L/D)",
            "α": "2.5 (average relative volatility)"
        },
        find=["Minimum stages (Fenske)", "Minimum reflux (Underwood)", "Actual stages (Gilliland)"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate Minimum Stages (Fenske Equation)",
                description="At total reflux, use the Fenske equation for minimum stages.",
                formula="N_min = ln[(xD/(1-xD)) × ((1-xB)/xB)] / ln(α)",
                substitution="N_min = ln[(0.95/0.05) × (0.95/0.05)] / ln(2.5)",
                computation="N_min = ln(361) / ln(2.5) = 5.89 / 0.916 = 6.43",
                result="N_min ≈ 6.4 theoretical stages",
                notes="This is the absolute minimum - only achievable at infinite reflux"
            ),
            CalculationStep(
                step_number=2,
                title="Estimate Minimum Reflux (Simplified Underwood)",
                description="Calculate the minimum reflux ratio below which separation is impossible.",
                formula="R_min ≈ (1/(α-1)) × (xD/xF - α(1-xD)/(1-xF))",
                substitution="R_min = (1/1.5) × (0.95/0.40 - 2.5×0.05/0.60)",
                computation="R_min = 0.667 × (2.375 - 0.208) = 0.667 × 2.167 = 1.44",
                result="R_min ≈ 1.44"
            ),
            CalculationStep(
                step_number=3,
                title="Apply Gilliland Correlation",
                description="Use the Gilliland correlation to find actual stages at R = 2.0.",
                formula="X = (R - R_min)/(R + 1), then Y = 1 - exp[(1+54.4X)/(11+117.2X) × (X-1)/√X]",
                substitution="X = (2.0 - 1.44)/(2.0 + 1) = 0.56/3.0 = 0.187",
                computation="Y = 1 - exp[...] ≈ 0.45 (from Gilliland curve)",
                result="Y ≈ 0.45"
            ),
            CalculationStep(
                step_number=4,
                title="Calculate Actual Stages",
                description="Convert Y value to actual theoretical stages.",
                formula="N = (N_min + Y)/(1 - Y)",
                substitution="N = (6.4 + 0.45)/(1 - 0.45)",
                computation="N = 6.85/0.55 = 12.5",
                result="N ≈ 12-13 theoretical stages",
                notes="Add 1 for reboiler, so ~13-14 stages including reboiler"
            )
        ],
        final_answer="Design requires approximately 12-13 theoretical stages (plus reboiler) at R = 2.0",
        real_world_context="Ethanol-water actually has a maximum-boiling azeotrope at 95.6% ethanol, limiting distillation. This example uses a simplified constant α."
    ),
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=20,
    prerequisites=["raoults_law", "antoine"],
    related_equations=["kremser", "fenske", "underwood"],
    references=[
        "McCabe, W.L. & Thiele, E.W. 'Graphical Design of Fractionating Columns', Ind. Eng. Chem., 1925",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley, 4th Ed.",
        "Perry's Chemical Engineers' Handbook, 9th Ed., Section 13"
    ],
    derivation_summary="McCabe-Thiele combines material balance equations with VLE equilibrium. The operating lines come from component balances around each section. Stepping between equilibrium and operating lines graphically solves the stage-by-stage calculation.",
    limitations_assumptions=[
        "Assumes constant molar overflow (CMO) - valid if ΔHvap similar",
        "Binary mixtures only - multicomponent requires other methods",
        "Theoretical stages - actual trays need efficiency correction (typically 50-80%)",
        "May fail near azeotropes where equilibrium curve crosses diagonal"
    ]
)


class McCabeThiele(BaseEquation):
    """McCabe-Thiele method for distillation stages."""
    
    equation_id = "mccabe_thiele"
    name = "McCabe-Thiele Stages"
    category = "Mass Transfer"
    description = "Estimate theoretical stages for binary distillation"
    reference = "Perry's Chemical Engineers' Handbook"
    learning_content = MCCABE_THIELE_LEARNING
    
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


# ===============================
# KREMSER LEARNING CONTENT
# ===============================
KREMSER_LEARNING = LearningContent(
    background_theory="""
The **Kremser equation** provides an analytical solution for calculating the number of 
theoretical stages required in absorption or stripping columns. Unlike the graphical 
McCabe-Thiele method, Kremser gives a direct algebraic formula.

**The Equation:**
N = ln[(y_in - y*)/(y_out - y*) × (1 - 1/A) + 1/A] / ln(A)

**Key Parameter - Absorption Factor (A):**
A = L / (m × G)

Where:
- **L** = Liquid molar flow rate
- **G** = Gas molar flow rate
- **m** = Henry's law constant slope (y = mx)
- **y_in** = Gas composition entering the column
- **y_out** = Gas composition leaving (target)
- **y*** = Equilibrium gas composition with entering liquid

**Physical Meaning:**
- **A > 1**: Absorption favored (more liquid capacity than needed)
- **A = 1**: Operating line parallel to equilibrium - infinite stages
- **A < 1**: Absorption limited - may not achieve target

The Kremser equation is the analytical equivalent of stepping off stages graphically.
""",
    key_concepts=[
        "Absorption factor A = L/(mG) controls effectiveness",
        "A > 1 needed for effective absorption",
        "Higher A → fewer stages but more liquid required",
        "Stripping uses S = 1/A = mG/L instead"
    ],
    real_world_applications=[
        "Natural gas sweetening (H2S, CO2 removal)",
        "Air pollution control (SO2 scrubbing)",
        "Solvent recovery from vent gases",
        "Stripping dissolved gases from water treatment"
    ],
    industry_examples=[
        "Amine absorber for CO2 capture in LNG plants",
        "TEG dehydration contactors for natural gas",
        "NH3 scrubbing in fertilizer plants"
    ],
    common_mistakes=[
        "Confusing A (absorption) with S (stripping) - they are reciprocals",
        "Using A < 1 and expecting good separation",
        "Forgetting to convert m to consistent units",
        "Not checking equilibrium data temperature/pressure basis"
    ],
    pro_tips=[
        "For absorption: target A = 1.2 to 2.0 for economic design",
        "For stripping: use S = mG/L > 1",
        "At A = 1 exactly, use simplified formula: N = (y_in - y_out)/(y_out - y*)",
        "Check Henry's law constant units carefully (pressure-based vs mole fraction)"
    ],
    variable_sources={
        "A": "Calculate A = L/(m*G) from measured liquid rate L, gas rate G, and equilibrium slope m from VLE data. Process control historian or DCS screens show flow rates.",
        "y_in": "Process feed gas composition from online analyzer (GC, IR) or laboratory grab sample. Check that analyzer calibration is current.",
        "y_out": "Target spec set by environmental permit, downstream process requirement, or product quality spec. Typically given as ppm or mole fraction.",
        "y*": "Equilibrium composition with entering liquid. For pure/regenerated solvent, y* typically equals 0. Check Henry's law: y* = H*x_liquid."
    },
    quiz_questions=[
        QuizQuestion(
            id="kr_q1",
            question="If the absorption factor A = 2.0, what does this mean physically?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Liquid can absorb twice the gas present", 
                     "Gas flow is twice liquid flow",
                     "Two stages are needed",
                     "Equilibrium is shifted toward liquid"],
            correct_answer="Liquid can absorb twice the gas present",
            explanation="A = L/(mG) > 1 means the liquid's capacity to dissolve solute exceeds what's in the gas. Higher A = easier absorption.",
            difficulty=Difficulty.BEGINNER,
            hint="Think about capacity - liquid vs gas",
            points=10
        ),
        QuizQuestion(
            id="kr_q2",
            question="What happens to the number of stages when A approaches 1?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Stages approach infinity", "Stages approach zero", 
                     "Stages equal to 1", "Stages become negative"],
            correct_answer="Stages approach infinity",
            explanation="At A = 1, the operating line is parallel to equilibrium, so you can never reach the target composition - infinite stages.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        ),
        QuizQuestion(
            id="kr_q3",
            question="For a stripper, you would use S = mG/L. If S = 0.5, is stripping effective?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["No, S should be > 1", "Yes, lower is better", 
                     "Depends on temperature", "Only at high pressure"],
            correct_answer="No, S should be > 1",
            explanation="For stripping, S > 1 is required just as A > 1 is needed for absorption. S < 1 means gas cannot carry away enough solute.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="CO2 Absorption in Amine Unit",
        scenario="An amine absorber removes CO2 from natural gas. The feed contains 5% CO2, and the spec is 50 ppm (0.005%). Pure amine enters (y* = 0). If A = 1.5, how many stages are needed?",
        given_values={
            "y_in": "0.05 (5% CO2)",
            "y_out": "0.00005 (50 ppm)",
            "y*": "0 (pure amine, no CO2)",
            "A": "1.5"
        },
        find=["Number of theoretical stages N"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate the separation ratio",
                description="Determine how much CO2 must be removed.",
                formula="φ = (y_in - y*)/(y_out - y*)",
                substitution="φ = (0.05 - 0)/(0.00005 - 0)",
                computation="φ = 0.05/0.00005 = 1000",
                result="φ = 1000 (need 1000x reduction)",
                notes="This is a deep removal - typical for gas treating"
            ),
            CalculationStep(
                step_number=2,
                title="Apply Kremser equation",
                description="Use the analytical formula for stages.",
                formula="N = ln[(φ - 1/A)/(1 - 1/A)] / ln(A)",
                substitution="N = ln[(1000 - 1/1.5)/(1 - 1/1.5)] / ln(1.5)",
                computation="N = ln[(1000 - 0.667)/(0.333)] / ln(1.5) = ln(2999) / 0.405",
                result="N = 8.0 / 0.405 = 19.8 stages"
            ),
            CalculationStep(
                step_number=3,
                title="Round and interpret",
                description="Report practical stage count.",
                computation="Round up to 20 theoretical stages",
                result="N = 20 theoretical stages",
                notes="Actual trays = N / efficiency (typically 20-30% for amine)"
            )
        ],
        final_answer="20 theoretical stages required for 5% → 50 ppm CO2 removal at A = 1.5",
        real_world_context="Real amine contactors typically have 15-25 actual trays with 20-30% tray efficiency, so this translates to many more actual trays."
    ),
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["raoults_law"],
    related_equations=["mccabe_thiele", "htu_ntu"],
    references=[
        "Kremser, A. 'Theoretical Analysis of Absorption Process', Nat. Petrol. News, 1930",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley, 4th Ed.",
        "Kohl, A. & Nielsen, R. 'Gas Purification', Gulf Publishing, 5th Ed."
    ],
    derivation_summary="Kremser derived the analytical solution by assuming linear equilibrium (y = mx) and constant molar flows. The result comes from solving the geometric series of stage-by-stage material balances.",
    limitations_assumptions=[
        "Assumes linear equilibrium relationship (y = mx)",
        "Constant molar overflow (liquid and gas rates constant)",
        "Single transferring component (other species inert)",
        "Isothermal operation (constant m)"
    ]
)


class Kremser(BaseEquation):
    """Kremser equation for absorption/stripping."""
    
    equation_id = "kremser"
    name = "Kremser Equation"
    category = "Mass Transfer"
    description = "Calculate stages for absorption or stripping"
    reference = "Seader & Henley, Separation Process Principles"
    learning_content = KREMSER_LEARNING
    
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


# ===============================
# HTU-NTU LEARNING CONTENT
# ===============================
HTU_NTU_LEARNING = LearningContent(
    background_theory="""
The **HTU-NTU method** is the standard approach for designing packed columns. Unlike tray 
columns with discrete stages, packed columns provide continuous contact between phases.

**Fundamental Relationship:**
Z = HTU × NTU

Where:
- **Z** = Required packing height (m or ft)
- **HTU** = Height of a Transfer Unit (characteristic of packing/system)
- **NTU** = Number of Transfer Units (determined by separation required)

**Calculating NTU (for dilute systems):**
NTU = ln[(y_in - y*)/(y_out - y*)]

Or using log-mean driving force:
NTU = (y_in - y_out) / ΔyLM

**HTU depends on:**
- Packing type and size
- Gas and liquid flow rates
- Physical properties (viscosity, diffusivity)
- Operating conditions

HTU is analogous to HETP (Height Equivalent to Theoretical Plate) but based on 
transfer unit concept rather than equilibrium stage.
""",
    key_concepts=[
        "Z = HTU × NTU separates packing characteristics from separation duty",
        "HTU determined by mass transfer coefficients and flow rates",
        "NTU determined by inlet/outlet concentrations and equilibrium",
        "Lower HTU = more efficient packing = shorter column"
    ],
    real_world_applications=[
        "Acid gas scrubbers (H2S, SO2, HCl removal)",
        "Packed distillation columns",
        "Cooling towers",
        "VOC emission control wet scrubbers"
    ],
    industry_examples=[
        "FGD (flue gas desulfurization) scrubbers in power plants",
        "Ammonia strippers in fertilizer production",
        "Structured packing in vacuum distillation"
    ],
    common_mistakes=[
        "Confusing HTU with HETP (different concepts)",
        "Not checking that driving force (y - y*) stays positive",
        "Using correlation-based HTU outside its valid range",
        "Forgetting to account for liquid-phase resistance"
    ],
    pro_tips=[
        "Typical HTU values: 0.3-1.0 m for structured packing, 0.5-1.5 m for random",
        "HTU_OG = HTU_G + (mG/L) × HTU_L for overall gas-phase",
        "If y* ≈ 0 (irreversible absorption), NTU = ln(y_in/y_out)",
        "Always add 20-50% safety factor to calculated height"
    ],
    variable_sources={
        "HTU": "From vendor correlations (Onda, Bravo-Fair), pilot plant tests, or literature for your packing type. Vendors often provide HETP which approximates HTU. Typical range 0.3-1.5 m.",
        "y_in": "Feed gas composition measured by online analyzer or laboratory sample. For combustion gases, use flue gas analyzer. Verify with mass balance.",
        "y_out": "Target outlet spec set by environmental permit, downstream process requirement, or product quality. Convert ppm to mole fraction: y = ppm / 1,000,000.",
        "y*": "Equilibrium composition with liquid leaving column bottom. For chemical absorption (reaction), y* is often near zero. For physical absorption, calculate from Henry's law."
    },
    quiz_questions=[
        QuizQuestion(
            id="htu_q1",
            question="If HTU = 0.5 m and NTU = 4, what packing height is needed?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="2",
            explanation="Z = HTU × NTU = 0.5 × 4 = 2.0 m",
            difficulty=Difficulty.BEGINNER,
            hint="Just multiply!",
            points=10
        ),
        QuizQuestion(
            id="htu_q2",
            question="Which packing type typically has lower HTU?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Structured packing", "Random packing (Raschig rings)", 
                     "They are always equal", "Depends only on column diameter"],
            correct_answer="Structured packing",
            explanation="Structured packing provides more surface area per volume and better liquid distribution, giving lower HTU (more efficient mass transfer).",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        ),
        QuizQuestion(
            id="htu_q3",
            question="For a 99% removal of a pollutant (y_in = 0.01, y_out = 0.0001, y* = 0), what is NTU?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="4.6",
            explanation="NTU = ln[(0.01 - 0)/(0.0001 - 0)] = ln(100) = 4.605 ≈ 4.6",
            difficulty=Difficulty.INTERMEDIATE,
            hint="NTU = ln(y_in/y_out) when y* = 0",
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="SO2 Scrubber Design",
        scenario="Design a packed scrubber to remove 95% of SO2 from flue gas. Inlet: 2000 ppm SO2. The scrubbing liquid (NaOH) is in large excess so y* ≈ 0. HTU for the system is 0.6 m.",
        given_values={
            "y_in": "2000 ppm = 0.002",
            "y_out": "100 ppm = 0.0001 (95% removal)",
            "y*": "0 (excess alkali)",
            "HTU": "0.6 m"
        },
        find=["NTU required", "Packing height Z"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate NTU",
                description="For y* = 0, the NTU formula simplifies.",
                formula="NTU = ln(y_in / y_out)",
                substitution="NTU = ln(0.002 / 0.0001)",
                computation="NTU = ln(20) = 3.0",
                result="NTU = 3.0 transfer units"
            ),
            CalculationStep(
                step_number=2,
                title="Calculate packing height",
                description="Apply the fundamental HTU-NTU relationship.",
                formula="Z = HTU × NTU",
                substitution="Z = 0.6 m × 3.0",
                computation="Z = 1.8 m",
                result="Z = 1.8 m required packing height"
            ),
            CalculationStep(
                step_number=3,
                title="Add safety factor",
                description="Apply engineering margin for uncertainties.",
                computation="Z_design = 1.8 × 1.3 = 2.34 m",
                result="Z_design ≈ 2.4 m (with 30% margin)",
                notes="Round to nearest standard packing bed depth"
            )
        ],
        final_answer="Design packing height = 2.4 m (including safety factor)",
        real_world_context="FGD scrubbers in power plants use similar calculations. Actual designs may have multiple spray zones and include mist eliminators."
    ),
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["kremser"],
    related_equations=["kremser", "film_coefficient", "packed_column_dp"],
    references=[
        "Treybal, R.E. 'Mass Transfer Operations', McGraw-Hill, 3rd Ed.",
        "Onda, Takeuchi & Okumoto, J. Chem. Eng. Japan, 1968",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley"
    ],
    derivation_summary="NTU is defined as the integral of dy/(y-y*) from outlet to inlet. For linear equilibrium and operating lines, this integral gives the logarithmic formula. HTU = G/(kya*a) relates mass transfer fundamentals to the empirical transfer unit concept.",
    limitations_assumptions=[
        "Assumes dilute concentrations (constant molar fluxes)",
        "Linear equilibrium relationship",
        "Plug flow of gas (no axial mixing)",
        "HTU is constant through the column"
    ]
)


class HTU_NTU(BaseEquation):
    """Height of Transfer Unit method for packed columns."""
    
    equation_id = "htu_ntu"
    name = "HTU-NTU Method"
    category = "Mass Transfer"
    description = "Calculate packed column height using HTU and NTU"
    reference = "Treybal, Mass Transfer Operations"
    learning_content = HTU_NTU_LEARNING
    
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


# ===============================
# FILM COEFFICIENT LEARNING CONTENT
# ===============================
FILM_COEFFICIENT_LEARNING = LearningContent(
    background_theory="""
The **mass transfer coefficient (k)** quantifies the rate of mass transfer between phases.
It's analogous to the heat transfer coefficient (h) in convection problems.

**Fundamental Definition:**
k = Sh × D_AB / L

Where:
- **k** = Mass transfer coefficient (m/s or mol/m²·s·Pa)
- **Sh** = Sherwood number (dimensionless mass transfer analog of Nusselt)
- **D_AB** = Molecular diffusivity of A in B (m²/s)
- **L** = Characteristic length (diameter, film thickness, etc.)

**Sherwood Number Correlations:**
- Laminar flow over flat plate: Sh = 0.664 × Re^0.5 × Sc^0.33
- Turbulent pipe flow: Sh = 0.023 × Re^0.8 × Sc^0.33
- Packed beds (Onda): Sh = 0.0051 × (Re_L)^0.7 × Sc_L^0.5 × (a_t × d_p)^(-2)

**Units Matter!**
- k_c = molar basis (mol/m²·s·(mol/m³)) = m/s
- k_G = pressure basis (mol/m²·s·Pa)
- k_y = mole fraction basis (mol/m²·s)
""",
    key_concepts=[
        "Sh = kL/D is the mass transfer analog of Nu = hL/k_thermal",
        "Sherwood depends on Re (flow) and Sc (fluid properties)",
        "Gas-phase coefficients typically larger than liquid-phase",
        "Overall coefficient combines gas and liquid film resistances"
    ],
    real_world_applications=[
        "Packed column design (kG, kL values)",
        "Membrane separation flux estimation",
        "Dissolution rate of solids",
        "Evaporation from liquid surfaces"
    ],
    common_mistakes=[
        "Using wrong characteristic length in Sh correlation",
        "Mixing up mass and molar diffusivity",
        "Applying correlation outside valid Re range",
        "Forgetting to convert units between different k bases"
    ],
    pro_tips=[
        "Sc = μ/(ρ·D) is typically 0.5-2 for gases, 100-1000 for liquids",
        "For gases in air at 25°C: D_AB ≈ 0.1-0.5 cm²/s",
        "For organics in water: D_AB ≈ 0.5-2 × 10⁻⁵ cm²/s",
        "1/k_overall = 1/k_G + m/k_L (series resistance)"
    ],
    variable_sources={
        "Sh": "Calculate from appropriate correlation for your geometry. For pipe flow: Sh = 0.023*Re^0.8*Sc^0.33. For packed beds: use Onda or similar. Check Re and Sc ranges.",
        "D_AB": "Diffusivity from Wilke-Chang (liquids) or Chapman-Enskog/FSG (gases) correlations. Perry's Handbook Table 5-16/5-17 has values. Units are m^2/s.",
        "L": "Characteristic length for your system: pipe diameter, packing size, particle diameter, or boundary layer thickness. Must match the correlation being used."
    },
    quiz_questions=[
        QuizQuestion(
            id="fc_q1",
            question="If Sh = 100, D = 2×10⁻⁵ m²/s, and L = 0.02 m, what is k?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="0.1",
            explanation="k = Sh × D / L = 100 × 2×10⁻⁵ / 0.02 = 0.1 m/s",
            difficulty=Difficulty.BEGINNER,
            hint="Substitute into k = Sh × D / L",
            points=10
        ),
        QuizQuestion(
            id="fc_q2",
            question="The Schmidt number (Sc = μ/ρD) is the mass transfer analog of which number?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Prandtl number", "Reynolds number", "Nusselt number", "Grashof number"],
            correct_answer="Prandtl number",
            explanation="Sc (momentum/mass diffusivity) is analogous to Pr (momentum/thermal diffusivity). Both appear in boundary layer correlations.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=["reynolds"],
    related_equations=["htu_ntu", "convection"],
    references=[
        "Sherwood, Pigford & Wilke, 'Mass Transfer', McGraw-Hill",
        "Bird, Stewart & Lightfoot, 'Transport Phenomena', Wiley",
        "Perry's Chemical Engineers' Handbook, Section 5"
    ],
    derivation_summary="The Sherwood number is defined by analogy to the Nusselt number. k = Sh × D/L comes from dimensional analysis of the convective-diffusive mass transport equation.",
    limitations_assumptions=[
        "Correlation-specific (check Re, Sc ranges)",
        "Assumes steady-state mass transfer",
        "Single transferring component",
        "Properties evaluated at film conditions"
    ]
)


class FilmCoefficient(BaseEquation):
    """Mass transfer film coefficient from correlations."""
    
    equation_id = "film_coefficient"
    name = "Film Mass Transfer Coefficient"
    category = "Mass Transfer"
    description = "Calculate gas or liquid film coefficient (kG or kL)"
    reference = "Sherwood, Pigford & Wilke"
    learning_content = FILM_COEFFICIENT_LEARNING
    
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


# ===============================
# LIQUID EXTRACTION LEARNING CONTENT
# ===============================
LIQUID_EXTRACTION_LEARNING = LearningContent(
    background_theory="""
**Liquid-liquid extraction** (also called solvent extraction) separates components based on 
their different solubilities in two immiscible liquid phases. It's widely used when 
distillation is impractical (e.g., heat-sensitive compounds, azeotropes).

**Key Concept - Distribution Coefficient:**
K_D = y/x = (concentration in extract)/(concentration in raffinate)

**Extraction Factor:**
E = K_D × (S/F)

Where:
- **K_D** = Distribution coefficient (partition ratio)
- **S** = Solvent flow rate
- **F** = Feed flow rate
- **x_F** = Solute concentration in feed
- **x_R** = Solute concentration in raffinate (what remains)

**The Kremser analog for extraction:**
N = ln[(x_F/x_R)(E-1) + 1] / ln(E)

**Physical Meaning:**
- **E > 1**: Extraction favored - solvent has excess capacity
- **E = 1**: Infinite stages required
- **E < 1**: Cannot achieve desired separation

Common extraction systems: water-organic, acid-base, metal chelation.
""",
    key_concepts=[
        "Distribution coefficient K_D from partition equilibrium",
        "Extraction factor E = K_D × (S/F) determines feasibility",
        "Counter-current flow gives best separation",
        "Choose solvent with high K_D for target compound"
    ],
    real_world_applications=[
        "Pharmaceutical purification (antibiotic extraction)",
        "Petrochemical aromatics extraction (sulfolane process)",
        "Metal recovery (copper from leach solutions)",
        "Food processing (decaffeination, flavor extraction)",
        "Nuclear fuel reprocessing (PUREX process)"
    ],
    industry_examples=[
        "Acetic acid recovery from aqueous waste using ethyl acetate",
        "Phenol extraction from wastewater using diisopropyl ether",
        "Rare earth separation using organic extractants"
    ],
    common_mistakes=[
        "Forgetting that K_D varies with concentration (use average)",
        "Not accounting for solvent in raffinate (mutual solubility)",
        "Ignoring temperature effects on K_D",
        "Using E < 1 and expecting good extraction"
    ],
    pro_tips=[
        "Higher S/F ratio → better extraction but more solvent cost",
        "Target E = 1.5 to 3.0 for economic multi-stage design",
        "Check mutual solubility - some 'immiscible' solvents dissolve significantly",
        "Consider raffinate recycling if valuable solute remains"
    ],
    variable_sources={
        "x_F": "Feed analysis from laboratory GC/HPLC sample or online analyzer. For process design, use typical operating data or pilot plant measurements.",
        "x_R": "Target raffinate concentration from environmental permit, product spec, or economic recovery analysis. Lower targets require more stages or solvent.",
        "K_D": "Distribution coefficient from shake-flask experiments at your temperature. Check literature for similar systems. K_D varies with concentration - use average.",
        "S/F": "Design variable - start at 1.0-2.0. Higher S/F gives fewer stages but more solvent cost. Optimize with economic analysis of solvent vs equipment cost."
    },
    quiz_questions=[
        QuizQuestion(
            id="lx_q1",
            question="If K_D = 3 and S/F = 2, what is the extraction factor E?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="6",
            explanation="E = K_D × (S/F) = 3 × 2 = 6. This is a favorable extraction (E >> 1).",
            difficulty=Difficulty.BEGINNER,
            hint="Just multiply K_D by S/F",
            points=10
        ),
        QuizQuestion(
            id="lx_q2",
            question="A solute has K_D = 0.5 in water/toluene. Which phase does it prefer?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["It prefers water (raffinate)", "It prefers toluene (extract)", 
                     "Equal distribution", "Cannot determine"],
            correct_answer="It prefers water (raffinate)",
            explanation="K_D = y/x < 1 means more solute stays in the raffinate (feed phase). The solute prefers water.",
            difficulty=Difficulty.BEGINNER,
            points=10
        ),
        QuizQuestion(
            id="lx_q3",
            question="Why might liquid extraction be preferred over distillation?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Heat-sensitive compounds", "Low relative volatility or azeotrope", 
                     "Very dilute streams", "All of the above"],
            correct_answer="All of the above",
            explanation="Extraction is mild (room temperature) and works on solubility differences, not volatility. It's great for thermally labile compounds, azeotropic mixtures, and recovering organics from dilute aqueous streams.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="Phenol Recovery from Wastewater",
        scenario="A wastewater stream contains 5% phenol and must be reduced to 0.1% before discharge. Using diisopropyl ether (DIPE) as solvent, K_D = 40. What solvent ratio and stages are needed?",
        given_values={
            "x_F": "0.05 (5% phenol in feed)",
            "x_R": "0.001 (0.1% target in raffinate)",
            "K_D": "40",
            "Target": "Minimize solvent use"
        },
        find=["Minimum S/F ratio", "Number of stages at S/F = 1.0"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate separation required",
                description="Determine the extraction ratio needed.",
                formula="Separation = x_F / x_R",
                substitution="Separation = 0.05 / 0.001",
                computation="Separation = 50× reduction required",
                result="Need 50× reduction in phenol concentration"
            ),
            CalculationStep(
                step_number=2,
                title="Calculate extraction factor at S/F = 1",
                description="E = K_D × (S/F)",
                substitution="E = 40 × 1.0",
                result="E = 40",
                notes="Very favorable - K_D is high for phenol/DIPE"
            ),
            CalculationStep(
                step_number=3,
                title="Apply Kremser equation for extraction",
                description="Calculate stages needed.",
                formula="N = ln[(x_F/x_R)(E-1) + 1] / ln(E)",
                substitution="N = ln[(50)(39) + 1] / ln(40)",
                computation="N = ln(1951) / ln(40) = 7.58 / 3.69",
                result="N ≈ 2.1 stages"
            ),
            CalculationStep(
                step_number=4,
                title="Practical design",
                description="Round up and add margin.",
                computation="Use 3 theoretical stages (mixer-settlers or column)",
                result="N = 3 stages at S/F = 1.0",
                notes="Very few stages because K_D is so high!"
            )
        ],
        final_answer="3 stages at S/F = 1.0 (equal solvent and feed flow) achieves 50× phenol reduction",
        real_world_context="Phenol is easily extracted because of its high K_D with DIPE. The solvent is then sent to a stripper to regenerate it and recover phenol."
    ),
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["kremser"],
    related_equations=["kremser", "mccabe_thiele"],
    references=[
        "Treybal, R.E. 'Liquid Extraction', McGraw-Hill, 2nd Ed.",
        "Lo, Baird & Hanson, 'Handbook of Solvent Extraction', Wiley",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley, 4th Ed."
    ],
    derivation_summary="The equation derives from material balance around each stage assuming counter-current flow and constant K_D. The geometric series solution is analogous to the Kremser equation for absorption.",
    limitations_assumptions=[
        "Assumes constant distribution coefficient (dilute solutions)",
        "Counter-current flow pattern",
        "Immiscible solvents (mutual solubility corrections may be needed)",
        "Theoretical stages - real efficiency typically 80-95%"
    ]
)


class LiquidExtraction(BaseEquation):
    """Liquid-liquid extraction stage calculation."""
    
    equation_id = "liquid_extraction"
    name = "Liquid-Liquid Extraction"
    category = "Mass Transfer"
    description = "Calculate extraction stages using distribution coefficient"
    reference = "Treybal, Liquid Extraction"
    learning_content = LIQUID_EXTRACTION_LEARNING
    
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


# ===============================
# MEMBRANE PERMEATION LEARNING CONTENT
# ===============================
MEMBRANE_PERMEATION_LEARNING = LearningContent(
    background_theory="""
**Membrane gas separation** uses selective permeation through polymer or inorganic membranes
to separate gas mixtures. The driving force is the partial pressure difference across the membrane.

**Permeation Equation:**
J = P × A × ΔP / δ

Where:
- **J** = Molar flux (mol/s or cm³(STP)/s)
- **P** = Permeability coefficient (Barrer)
- **A** = Membrane area (m² or cm²)
- **ΔP** = Partial pressure difference (bar or cmHg)
- **δ** = Membrane thickness (μm)

**The Barrer Unit:**
1 Barrer = 10⁻¹⁰ cm³(STP)·cm/(cm²·s·cmHg)

**Selectivity (α):**
α = P_A / P_B = (D_A/D_B) × (S_A/S_B)

Where D = diffusivity and S = solubility in the membrane.

**Tradeoff:** High permeability membranes often have low selectivity (Robeson upper bound).
""",
    key_concepts=[
        "Permeability = Diffusivity × Solubility (solution-diffusion model)",
        "Selectivity determines separation quality",
        "Thinner membranes = higher flux but mechanical limits",
        "Pressure ratio affects stage cut and purity"
    ],
    real_world_applications=[
        "Nitrogen generation from air",
        "Hydrogen recovery in refineries",
        "Natural gas CO2 removal",
        "Biogas upgrading (CO2/CH4 separation)"
    ],
    common_mistakes=[
        "Ignoring pressure ratio effects on separation",
        "Using Barrer values at wrong temperature",
        "Forgetting that flux is per unit area",
        "Not accounting for concentration polarization"
    ],
    pro_tips=[
        "Typical industrial permeances: 100-1000 GPU (1 GPU = 10⁻⁶ cm³/cm²·s·cmHg)",
        "For asymmetric membranes, use permeance (P/δ) not permeability",
        "Multi-stage systems needed for high-purity products",
        "Check Robeson plot for realistic P vs α combinations"
    ],
    variable_sources={
        "P": "Permeability from membrane vendor data sheets or literature (Robeson 2008, Baker textbook). Units in Barrer. Check temperature - P increases with T.",
        "delta": "Membrane selective layer thickness from vendor specification. For asymmetric membranes, this is the thin skin layer (0.1-1 um), not total thickness.",
        "delta_p": "Partial pressure difference across membrane. Calculate from feed/permeate pressures and compositions. Higher ratio gives better separation but costs more."
    },
    quiz_questions=[
        QuizQuestion(
            id="mp_q1",
            question="If membrane selectivity α = Pₐ/Pᵦ = 20, and CO2 is 'A', what does this mean?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["CO2 permeates 20x faster than the other gas", 
                     "The membrane is 20 μm thick",
                     "20 bar pressure is needed",
                     "Recovery is 20%"],
            correct_answer="CO2 permeates 20x faster than the other gas",
            explanation="Selectivity is the ratio of permeabilities. α = 20 means CO2 moves through the membrane 20 times faster than the slower component.",
            difficulty=Difficulty.BEGINNER,
            points=10
        ),
        QuizQuestion(
            id="mp_q2",
            question="What is the 'Robeson upper bound' in membrane science?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["The maximum achievable combination of permeability and selectivity", 
                     "The maximum operating pressure for polymeric membranes",
                     "The upper limit on membrane area",
                     "The thickest practical membrane"],
            correct_answer="The maximum achievable combination of permeability and selectivity",
            explanation="Robeson showed there's a tradeoff: high permeability membranes have low selectivity and vice versa. The 'upper bound' line represents the best achievable combinations.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=[],
    related_equations=["film_coefficient"],
    references=[
        "Baker, R.W. 'Membrane Technology and Applications', Wiley, 3rd Ed.",
        "Robeson, L.M. 'The upper bound revisited', J. Membrane Sci., 2008",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley"
    ],
    derivation_summary="The solution-diffusion model: gas dissolves in membrane surface (S), diffuses through (D), and desorbs. Permeability P = D × S. Flux is proportional to P and ΔP, inversely proportional to thickness.",
    limitations_assumptions=[
        "Solution-diffusion mechanism (not Knudsen flow)",
        "Steady-state permeation",
        "Ideal mixing on both sides",
        "No concentration polarization"
    ]
)


class MembranePermeation(BaseEquation):
    """Membrane gas permeation calculation."""
    
    equation_id = "membrane_permeation"
    name = "Membrane Permeation"
    category = "Mass Transfer"
    description = "Calculate membrane flux and selectivity"
    reference = "Baker, Membrane Technology"
    learning_content = MEMBRANE_PERMEATION_LEARNING
    
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


# ===============================
# PACKED COLUMN ΔP LEARNING CONTENT
# ===============================
PACKED_COLUMN_DP_LEARNING = LearningContent(
    background_theory="""
The **Ergun equation** calculates pressure drop through packed beds. It applies to packed
columns, fixed-bed reactors, and filtration systems.

**Ergun Equation:**
ΔP/L = 150 × μ(1-ε)²/dp²ε³ × v + 1.75 × ρ(1-ε)/dp·ε³ × v²

Where:
- **ΔP** = Pressure drop (Pa)
- **L** = Bed height (m)
- **μ** = Fluid viscosity (Pa·s)
- **ε** = Void fraction (typically 0.35-0.45 for random packing)
- **dp** = Particle/packing diameter (m)
- **v** = Superficial velocity (m/s)
- **ρ** = Fluid density (kg/m³)

**The Two Terms:**
1. **First term (150...)**: Viscous losses (dominates at low Re)
2. **Second term (1.75...)**: Inertial losses (dominates at high Re)

**Flooding:** At high gas rates, the column "floods" — liquid cannot drain and pressure drop
spikes. Design typically targets 70-80% of flood velocity.
""",
    key_concepts=[
        "ΔP increases with (1-ε)² term - void fraction is critical",
        "Smaller packing = higher ΔP but better mass transfer",
        "At flooding, ΔP increases dramatically",
        "Low ε (tight packing) → much higher ΔP"
    ],
    real_world_applications=[
        "Packed column design (absorption, distillation)",
        "Fixed-bed reactor pressure drop",
        "Filter bed sizing",
        "Fluidization calculations"
    ],
    common_mistakes=[
        "Using wrong void fraction (varies with column/packing diameter ratio)",
        "Ignoring liquid holdup contribution to ΔP",
        "Not checking flooding limits",
        "Using superficial velocity instead of interstitial"
    ],
    pro_tips=[
        "Typical void fractions: random packing 0.35-0.45, structured 0.90-0.97",
        "For two-phase flow, use Generalized Pressure Drop Correlation (GPDC)",
        "Design at 70-80% of flood velocity for safe operation",
        "ΔP increases roughly with the square of gas velocity"
    ],
    variable_sources={
        "epsilon": "Void fraction from packing vendor data. Measure by water displacement for installed bed. Typical random packing 0.35-0.45, structured 0.90-0.97.",
        "dp": "Packing nominal diameter from vendor catalog. For non-spherical particles, use equivalent diameter (6*Vp/Ap). Units must match other lengths.",
        "mu": "Gas dynamic viscosity from property correlations (Sutherland, Chapman-Enskog). Perry's Table 2-305 has values. Check temperature dependence.",
        "rho": "Gas density from ideal gas law or equation of state at operating T, P. For mixtures, use molar average molecular weight."
    },
    quiz_questions=[
        QuizQuestion(
            id="dp_q1",
            question="If you halve the packing diameter dp, what happens to pressure drop (approximately)?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["ΔP roughly quadruples (×4)", "ΔP doubles (×2)", 
                     "ΔP halves (×0.5)", "ΔP stays the same"],
            correct_answer="ΔP roughly quadruples (×4)",
            explanation="ΔP scales with 1/dp² in the viscous term and 1/dp in the inertial term. For small particles (viscous-dominated), halving dp quadruples ΔP.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        ),
        QuizQuestion(
            id="dp_q2",
            question="At what point does 'flooding' occur in a packed column?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["When gas velocity prevents liquid drainage", 
                     "When liquid flow is too low",
                     "When temperature exceeds boiling point",
                     "When pressure reaches 1 atm"],
            correct_answer="When gas velocity prevents liquid drainage",
            explanation="Flooding occurs when upward gas flow is so high that liquid cannot drain down. The liquid accumulates, void fraction drops, and ΔP spikes.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=["reynolds"],
    related_equations=["htu_ntu", "darcy_weisbach"],
    references=[
        "Ergun, S. 'Fluid Flow Through Packed Columns', Chem. Eng. Prog., 1952",
        "Perry's Chemical Engineers' Handbook, Section 6",
        "Seader, Henley & Roper, 'Separation Process Principles', Wiley"
    ],
    derivation_summary="Ergun combined the Blake-Kozeny equation (viscous, low Re) with the Burke-Plummer equation (inertial, high Re) into a single correlation valid across all flow regimes.",
    limitations_assumptions=[
        "Single-phase flow (correlations exist for two-phase)",
        "Uniform packing (no channeling)",
        "Spherical or near-spherical particles",
        "Incompressible fluid"
    ]
)


class PackedColumnPressureDrop(BaseEquation):
    """Pressure drop in packed columns (Ergun equation)."""
    
    equation_id = "packed_column_dp"
    name = "Packed Column Pressure Drop"
    category = "Mass Transfer"
    description = "Ergun equation for packed bed pressure drop"
    reference = "Ergun, 1952"
    learning_content = PACKED_COLUMN_DP_LEARNING
    
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
