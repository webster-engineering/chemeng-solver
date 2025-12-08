"""
Thermodynamics equations - ideal gas, flash, VLE calculations.
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
# IDEAL GAS LEARNING CONTENT
# ===============================
IDEAL_GAS_LEARNING = LearningContent(
    background_theory="""
The Ideal Gas Law (PV = nRT) is one of the most fundamental equations in thermodynamics 
and chemical engineering. It relates four key properties of a gas:

**The Equation:** PV = nRT

Where:
- **P** = Absolute pressure (not gauge!)
- **V** = Volume occupied by the gas
- **n** = Amount of substance (moles)
- **R** = Universal gas constant (8.314 J/mol·K)
- **T** = Absolute temperature (must be Kelvin or Rankine!)

**Key Assumptions:**
1. Gas molecules have negligible volume compared to the container
2. No intermolecular forces between molecules
3. Collisions are perfectly elastic
4. Random molecular motion

The ideal gas law works best at **low pressures** and **high temperatures** where 
real gases behave most ideally.
""",
    key_concepts=[
        "Absolute vs gauge pressure",
        "Absolute temperature scales (Kelvin, Rankine)",
        "Molar quantities and molecular weight",
        "The universal gas constant R has different values in different unit systems"
    ],
    real_world_applications=[
        "Sizing gas storage tanks and vessels",
        "Calculating gas flow rates in pipelines", 
        "HVAC system design and air handling",
        "Pneumatic system design",
        "Chemical reactor sizing for gas-phase reactions",
        "Estimating gas volumes at different conditions"
    ],
    industry_examples=[
        "Compressed air systems in manufacturing plants",
        "Natural gas processing and transmission",
        "Pharmaceutical lyophilization (freeze-drying)",
        "Semiconductor manufacturing clean rooms",
        "Automotive airbag inflation calculations"
    ],
    common_mistakes=[
        "Using gauge pressure instead of absolute pressure",
        "Using Celsius or Fahrenheit instead of absolute temperature (K or °R)",
        "Forgetting to convert mass to moles using molecular weight",
        "Using inconsistent units for the gas constant R",
        "Applying to high-pressure systems where real gas effects matter"
    ],
    pro_tips=[
        "Always add atmospheric pressure to gauge readings: P_abs = P_gauge + P_atm",
        "Convert °C to K by adding 273.15; convert °F to °R by adding 459.67",
        "At conditions far from ideal (high P, low T), use Van der Waals or other real gas equations",
        "R = 8.314 J/(mol·K) = 10.73 (psia·ft³)/(lbmol·°R) = 0.08206 (L·atm)/(mol·K)"
    ],
    variable_sources={
        "n": "Calculate from mass and molecular weight: n = mass / MW. Look up MW in chemical tables or material data sheets.",
        "T": "Measure with thermometer or thermocouple. Process control systems display this. Must convert to absolute scale!",
        "P": "Read from pressure gauge and ADD atmospheric pressure (~14.7 psia or 101.3 kPa). Check if gauge or absolute.",
        "V": "For unknowns: this is what you're solving for. For verification: measure tank dimensions or use vessel specs."
    },
    quiz_questions=[
        QuizQuestion(
            id="ig_q1",
            question="What type of pressure must be used in the ideal gas law?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Gauge pressure", "Absolute pressure", "Differential pressure", "Vacuum pressure"],
            correct_answer="Absolute pressure",
            explanation="The ideal gas law requires absolute pressure. Gauge pressure reads zero at atmospheric conditions, but there's still ~14.7 psi of absolute pressure present. Always use P_abs = P_gauge + P_atm.",
            difficulty=Difficulty.BEGINNER,
            hint="Think about what happens at 'zero' gauge pressure - is there really no pressure?",
            points=10
        ),
        QuizQuestion(
            id="ig_q2", 
            question="A gas at 25°C needs to be used in the ideal gas law. What temperature should you use?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="298.15",
            explanation="Convert Celsius to Kelvin: K = °C + 273.15 = 25 + 273.15 = 298.15 K. The ideal gas law requires absolute temperature!",
            difficulty=Difficulty.BEGINNER,
            hint="Add 273.15 to convert from Celsius to Kelvin",
            points=10
        ),
        QuizQuestion(
            id="ig_q3",
            question="If you double the absolute temperature of an ideal gas at constant pressure and moles, what happens to the volume?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Volume doubles", "Volume halves", "Volume stays the same", "Volume quadruples"],
            correct_answer="Volume doubles",
            explanation="From PV = nRT, if P and n are constant: V ∝ T. Doubling T doubles V. This is Charles's Law.",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Rearrange PV = nRT to solve for V in terms of T",
            points=15
        ),
        QuizQuestion(
            id="ig_q4",
            question="The ideal gas law works best under which conditions?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["High pressure, low temperature", "Low pressure, high temperature", 
                     "High pressure, high temperature", "Low pressure, low temperature"],
            correct_answer="Low pressure, high temperature",
            explanation="At low pressure, gas molecules are far apart (negligible volume). At high temperature, kinetic energy dominates over intermolecular forces. These conditions match the ideal gas assumptions.",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Think about when gas molecules are most spread out and moving fastest",
            points=15
        ),
        QuizQuestion(
            id="ig_q5",
            question="What is the value of R in SI units (J/mol·K)?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="8.314",
            explanation="R = 8.314 J/(mol·K) is the universal gas constant. Other common values: 0.08206 L·atm/(mol·K), 10.73 psia·ft³/(lbmol·°R).",
            difficulty=Difficulty.BEGINNER,
            hint="This is a fundamental constant worth memorizing",
            points=10
        )
    ],
    worked_example=WorkedExample(
        title="Compressed Air Tank Sizing",
        scenario="A chemical plant needs to store compressed air for pneumatic instruments. The tank will hold 5 kg of air at 25°C and 150 psig. What tank volume is required?",
        given_values={
            "mass": "5 kg",
            "T": "25°C = 298.15 K",
            "P_gauge": "150 psig",
            "P_atm": "14.7 psia",
            "MW_air": "28.97 g/mol"
        },
        find=["Tank volume V in cubic feet"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Convert gauge to absolute pressure",
                description="The pressure gauge reads 150 psig. We must add atmospheric pressure to get absolute pressure.",
                formula="P_abs = P_gauge + P_atm",
                substitution="P_abs = 150 + 14.7",
                computation="P_abs = 164.7 psia",
                result="P = 164.7 psia",
                notes="Always check if pressure is gauge (psig) or absolute (psia)!"
            ),
            CalculationStep(
                step_number=2,
                title="Calculate moles of air",
                description="Convert mass to moles using the molecular weight of air.",
                formula="n = mass / MW",
                substitution="n = 5000 g / 28.97 g/mol",
                computation="n = 172.6 mol",
                result="n = 172.6 mol",
                notes="Air MW ≈ 29 g/mol (78% N₂ at 28 + 21% O₂ at 32)"
            ),
            CalculationStep(
                step_number=3,
                title="Verify temperature is absolute",
                description="Convert Celsius to Kelvin for use in the ideal gas law.",
                formula="T(K) = T(°C) + 273.15",
                substitution="T = 25 + 273.15",
                computation="T = 298.15 K",
                result="T = 298.15 K",
                notes="Never use °C or °F directly in gas law equations!"
            ),
            CalculationStep(
                step_number=4,
                title="Apply the ideal gas law",
                description="Rearrange PV = nRT to solve for volume V.",
                formula="V = nRT / P",
                substitution="V = (172.6 mol)(8.314 J/mol·K)(298.15 K) / (164.7 psia)",
                computation="Converting: 164.7 psia = 1,135,500 Pa. V = (172.6)(8.314)(298.15) / 1,135,500 = 0.377 m³",
                result="V = 0.377 m³ = 13.3 ft³",
                notes="1 m³ = 35.31 ft³"
            )
        ],
        final_answer="V ≈ 13.3 ft³ (0.377 m³)",
        real_world_context="In practice, you'd specify a standard tank size (e.g., 15 ft³) with some margin for safety and operational flexibility."
    ),
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=15,
    prerequisites=[],
    related_equations=["van_der_waals", "compressor_work"],
    diagram_type="gas_tank",
    diagram_labels={
        "tank": "Gas Storage Vessel",
        "pressure_gauge": "Pressure Gauge (P)",
        "thermometer": "Temperature Sensor (T)",
        "volume": "Internal Volume (V)",
        "molecules": "Gas Molecules (n moles)"
    }
)


class IdealGas(BaseEquation):
    """Ideal gas law calculations."""
    
    equation_id = "ideal_gas"
    name = "Ideal Gas Law (PV=nRT)"
    category = "Thermodynamics"
    description = "Calculate any one unknown from P, V, n, T using PV = nRT"
    reference = "Standard thermodynamics"
    learning_content = IDEAL_GAS_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("n", "Moles of gas - amount of substance", 
                              "dimensionless", "mol", ParameterType.INPUT, symbol="n",
                              tooltip="n = mass / molecular weight"),
            EquationParameter("T", "Absolute temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="T",
                              tooltip="Must use absolute scale internally (R or K)"),
            EquationParameter("P", "Absolute pressure", "pressure", "psia", 
                              ParameterType.INPUT, symbol="P",
                              tooltip="Absolute, not gauge pressure"),
            EquationParameter("V", "Volume occupied by gas", "volume", "ft**3", 
                              ParameterType.OUTPUT, symbol="V"),
        ]
    
    def _calculate(self, n: pint.Quantity, T: pint.Quantity,
                   P: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        n_val = n.to('mol').magnitude
        t = T.to('K').magnitude
        p = P.to('Pa').magnitude
        
        R = 8.314  # J/(mol·K)
        v = n_val * R * t / p  # m³
        v_ft3 = v * 35.3147
        
        return {"V": ureg.Quantity(v_ft3, "ft**3")}


# ===============================
# ANTOINE EQUATION LEARNING CONTENT  
# ===============================
ANTOINE_LEARNING = LearningContent(
    background_theory="""
The Antoine Equation is an empirical correlation for calculating the **saturation vapor pressure** 
of pure substances as a function of temperature. It's essential for VLE calculations.

**The Equation:** log₁₀(P) = A - B/(C + T)

Where:
- **P** = Vapor pressure (typically mmHg or bar)
- **T** = Temperature (typically °C)
- **A, B, C** = Compound-specific constants from tables

**Why It Matters:**
- Used to calculate bubble point and dew point temperatures
- Essential for distillation column design
- Foundation of all VLE calculations

The Antoine equation is more accurate than Clausius-Clapeyron for moderate temperature ranges 
(typically within 30-40°C of the normal boiling point).
""",
    key_concepts=[
        "Vapor-liquid equilibrium (VLE)",
        "Saturation (vapor) pressure increases exponentially with temperature",
        "Antoine constants are specific to each compound",
        "Units must match the tabulated constants"
    ],
    real_world_applications=[
        "Distillation column design and operation",
        "Flash drum calculations",
        "Evaporator sizing",
        "Environmental emission calculations",
        "Refrigeration system design"
    ],
    industry_examples=[
        "Determining overhead condenser duty in crude distillation",
        "Sizing vacuum systems for low-pressure distillation",
        "Calculating VOC emissions from storage tanks"
    ],
    common_mistakes=[
        "Using wrong temperature units (most constants are for °C)",
        "Mixing Antoine coefficients from different sources with different pressure units",
        "Applying outside the valid temperature range of the constants",
        "Confusing log₁₀ with natural log (ln)"
    ],
    pro_tips=[
        "NIST Chemistry WebBook is the most reliable source for Antoine constants",
        "Check the validity range - constants are only accurate within a specific temperature range",
        "Convert pressure units after calculation, not coefficients",
        "For extended ranges, consider Wagner or Riedel equations"
    ],
    variable_sources={
        "A": "Look up in NIST, Perry's, or CRC Handbook for your compound",
        "B": "Same source as A - always use A,B,C from the same source",
        "C": "Same source as A,B - these three are fitted together",
        "T": "Your operating temperature - check it's within valid range"
    },
    quiz_questions=[
        QuizQuestion(
            id="ant_q1",
            question="What type of pressure does the Antoine equation calculate?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Saturation (vapor) pressure", "Total system pressure", "Partial pressure", "Critical pressure"],
            correct_answer="Saturation (vapor) pressure",
            explanation="Antoine calculates the pressure at which a pure liquid boils at a given temperature.",
            difficulty=Difficulty.BEGINNER,
            hint="Think about what happens at the boiling point",
            points=10
        ),
        QuizQuestion(
            id="ant_q2",
            question="Why can't you mix Antoine constants from different sources?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Different sources use different units", "The constants are fitted together to minimize error", 
                     "Some sources have typos", "It's actually fine to mix them"],
            correct_answer="The constants are fitted together to minimize error",
            explanation="A, B, and C are regressed together. Using A from one source and B from another will give wrong results.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        ),
        QuizQuestion(
            id="ant_q3",
            question="Water has A=8.07, B=1730.6, C=233.4 (T in °C, P in mmHg). What's the vapor pressure at 100°C?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="760",
            explanation="log₁₀(P) = 8.07 - 1730.6/(233.4+100) = 2.88, so P = 10^2.88 ≈ 760 mmHg (1 atm at boiling!)",
            difficulty=Difficulty.INTERMEDIATE,
            hint="100°C is water's normal boiling point...",
            points=20
        )
    ],
    worked_example=WorkedExample(
        title="Vapor Pressure of Benzene",
        scenario="Calculate the vapor pressure of benzene at 50°C for storage tank vent design.",
        given_values={
            "A": "6.90565",
            "B": "1211.033",
            "C": "220.79",
            "T": "50°C",
            "Units": "P in mmHg, T in °C"
        },
        find=["Vapor pressure in mmHg and psia"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Apply Antoine Equation",
                description="Substitute values into Antoine equation",
                formula="log₁₀(P) = A - B/(C + T)",
                substitution="log₁₀(P) = 6.90565 - 1211.033/(220.79 + 50)",
                computation="log₁₀(P) = 6.90565 - 1211.033/270.79 = 6.90565 - 4.472 = 2.433",
                result="log₁₀(P) = 2.433"
            ),
            CalculationStep(
                step_number=2,
                title="Convert to Pressure",
                description="Take antilog to get pressure in mmHg",
                formula="P = 10^(log₁₀P)",
                substitution="P = 10^2.433",
                computation="P = 271 mmHg",
                result="P = 271 mmHg"
            ),
            CalculationStep(
                step_number=3,
                title="Convert Units",
                description="Convert to psia for US engineering use",
                formula="P(psia) = P(mmHg) × 14.696/760",
                substitution="P = 271 × 0.01934",
                computation="P = 5.24 psia",
                result="P = 5.24 psia",
                notes="This is above atmospheric, so benzene will readily evaporate at 50°C"
            )
        ],
        final_answer="P_vap = 271 mmHg = 5.24 psia",
        real_world_context="At 50°C, benzene has significant vapor pressure - tanks need proper venting!"
    ),
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=12,
    prerequisites=[],
    related_equations=["raoults_law", "clausius_clapeyron", "flash"],
    diagram_type="phase_diagram",
    diagram_labels={
        "curve": "Vapor Pressure Curve",
        "liquid_region": "Liquid Phase",
        "vapor_region": "Vapor Phase",
        "boiling_point": "Normal Boiling Point"
    }
)


# ===============================
# RAOULT'S LAW LEARNING CONTENT
# ===============================
RAOULTS_LAW_LEARNING = LearningContent(
    background_theory="""
Raoult's Law describes the vapor-liquid equilibrium (VLE) for **ideal solutions**. It states that 
the partial pressure of a component equals its vapor pressure times its mole fraction in the liquid.

**The Equation:** yᵢ·P = xᵢ·Pᵢˢᵃᵗ  →  yᵢ = xᵢ·Pᵢˢᵃᵗ/P

Where:
- **yᵢ** = Vapor mole fraction of component i
- **xᵢ** = Liquid mole fraction of component i  
- **P** = Total system pressure
- **Pᵢˢᵃᵗ** = Pure component vapor pressure at T (from Antoine)

**The K-value:**  K = y/x = Pˢᵃᵗ/P
- K > 1: Component concentrates in vapor (more volatile)
- K < 1: Component concentrates in liquid (less volatile)
- K = 1: Equal distribution (at component's boiling point at P)

Raoult's Law works for **ideal mixtures** - similar molecules like alkane mixtures.
""",
    key_concepts=[
        "Ideal solution behavior - similar molecules",
        "K-values (equilibrium ratios) for separation design",
        "Relative volatility (α = K₁/K₂) drives distillation",
        "Non-ideal systems require activity coefficients"
    ],
    real_world_applications=[
        "Binary distillation column design",
        "Flash drum vapor-liquid split",
        "Bubble point and dew point calculations",
        "Gasoline blending (hydrocarbon mixtures)"
    ],
    common_mistakes=[
        "Applying to highly non-ideal mixtures (alcohols + water)",
        "Forgetting that Pˢᵃᵗ depends on temperature",
        "Not verifying x and y sum to 1.0",
        "Using gauge instead of absolute pressure"
    ],
    pro_tips=[
        "Calculate relative volatility α = K₁/K₂ - this determines separation difficulty",
        "For non-ideal systems, multiply by activity coefficient: y·P = γ·x·Pˢᵃᵗ",
        "If K·x > 1, the calculation is inconsistent - check your inputs"
    ],
    quiz_questions=[
        QuizQuestion(
            id="raoult_q1",
            question="If K > 1 for a component, where does it concentrate?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["In the vapor phase", "In the liquid phase", "Equally distributed", "Cannot determine"],
            correct_answer="In the vapor phase",
            explanation="K = y/x, so K > 1 means y > x, the component is richer in vapor.",
            difficulty=Difficulty.BEGINNER,
            points=10
        ),
        QuizQuestion(
            id="raoult_q2",
            question="Component A has Pˢᵃᵗ = 2 atm and x = 0.3 at a total pressure of 1 atm. What is y?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="0.6",
            explanation="K = Pˢᵃᵗ/P = 2/1 = 2. y = K·x = 2 × 0.3 = 0.6",
            difficulty=Difficulty.INTERMEDIATE,
            hint="First calculate K, then multiply by x",
            points=15
        )
    ],
    worked_example=WorkedExample(
        title="Benzene-Toluene Flash Separation",
        scenario="A 50/50 benzene-toluene mixture is flashed at 1 atm and 95°C. Calculate vapor composition.",
        given_values={
            "x_benzene": "0.5",
            "Pˢᵃᵗ_benzene": "1180 mmHg (at 95°C)",
            "Pˢᵃᵗ_toluene": "480 mmHg (at 95°C)",
            "P_total": "760 mmHg (1 atm)"
        },
        find=["Vapor composition (y_benzene)"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Calculate K-values",
                description="K = Pˢᵃᵗ/P for each component",
                formula="K = Pˢᵃᵗ / P",
                substitution="K_benzene = 1180/760, K_toluene = 480/760",
                computation="K_benzene = 1.55, K_toluene = 0.63",
                result="K_B = 1.55, K_T = 0.63"
            ),
            CalculationStep(
                step_number=2,
                title="Apply Raoult's Law",
                description="Calculate vapor mole fractions",
                formula="y = K · x",
                substitution="y_benzene = 1.55 × 0.5 = 0.78",
                computation="y_toluene = 0.63 × 0.5 = 0.32",
                result="y_B = 0.78, y_T = 0.22 (normalized)",
                notes="Note: benzene (K>1) enriches in vapor"
            )
        ],
        final_answer="y_benzene = 0.78 (78% benzene in vapor)",
        real_world_context="This is why benzene comes off the top of a BTX splitter!"
    ),
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=10,
    prerequisites=["antoine"],
    related_equations=["antoine", "flash", "activity_coefficient"]
)


# ===============================
# CLAUSIUS-CLAPEYRON LEARNING CONTENT
# ===============================
CLAUSIUS_CLAPEYRON_LEARNING = LearningContent(
    background_theory="""
The Clausius-Clapeyron equation relates vapor pressure to temperature using the **heat of vaporization**.
It's derived from thermodynamic first principles and is useful when Antoine constants aren't available.

**The Equation:** ln(P₂/P₁) = (ΔHᵥₐₚ/R) × (1/T₁ - 1/T₂)

Where:
- **P₁, P₂** = Vapor pressures at T₁ and T₂
- **ΔHᵥₐₚ** = Heat (enthalpy) of vaporization
- **R** = Gas constant (8.314 J/mol·K)
- **T₁, T₂** = Absolute temperatures (Kelvin!)

**Key Insight:** A plot of ln(P) vs 1/T gives a straight line with slope = -ΔHᵥₐₚ/R

This is the theoretical basis for the Antoine equation (which adds the C parameter for better fit).
""",
    key_concepts=[
        "Heat of vaporization (latent heat) drives the P-T relationship",
        "Must use absolute temperature (K or R)",
        "Assumes ΔHᵥₐₚ is constant (valid for small T ranges)",
        "Natural log (ln), not log₁₀"
    ],
    real_world_applications=[
        "Estimating vapor pressure without Antoine data",
        "Calculating boiling point at different pressures",
        "Vacuum distillation design",
        "Understanding altitude effects on boiling"
    ],
    common_mistakes=[
        "Using °C or °F instead of absolute temperature",
        "Confusing units of ΔHᵥₐₚ (J/mol vs kJ/kg)",
        "Applying over too wide a temperature range",
        "Using log₁₀ instead of natural log"
    ],
    quiz_questions=[
        QuizQuestion(
            id="cc_q1",
            question="What happens to boiling point when you reduce pressure?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Boiling point decreases", "Boiling point increases", "No change", "Depends on the substance"],
            correct_answer="Boiling point decreases",
            explanation="Lower pressure = lower vapor pressure needed to boil = lower temperature.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    worked_example=WorkedExample(
        title="Water Boiling Point at Altitude",
        scenario="Calculate the boiling point of water at 0.8 atm (approximately Denver, CO altitude).",
        given_values={
            "P₁": "1 atm = 101.3 kPa",
            "T₁": "100°C = 373.15 K",
            "P₂": "0.8 atm = 81.1 kPa",
            "ΔHᵥₐₚ": "40,650 J/mol"
        },
        find=["Boiling point T₂"],
        steps=[
            CalculationStep(
                step_number=1,
                title="Rearrange for T₂",
                description="Solve Clausius-Clapeyron for the unknown temperature",
                formula="1/T₂ = 1/T₁ - (R/ΔHᵥₐₚ)·ln(P₂/P₁)",
                substitution="1/T₂ = 1/373.15 - (8.314/40650)·ln(0.8)",
                computation="1/T₂ = 0.002680 - (-4.56×10⁻⁵) = 0.002726",
                result="1/T₂ = 0.002726"
            ),
            CalculationStep(
                step_number=2,
                title="Calculate T₂",
                description="Take reciprocal and convert to Celsius",
                formula="T₂ = 1/0.002726 K",
                computation="T₂ = 366.9 K = 93.7°C",
                result="T₂ = 93.7°C",
                notes="Water boils about 6°C lower in Denver!"
            )
        ],
        final_answer="T₂ = 93.7°C (water boils at lower temperature at altitude)",
        real_world_context="This is why cooking times increase at high altitude - water boils at lower temperature!"
    ),
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["ideal_gas"],
    related_equations=["antoine", "heat_capacity_mixture"]
)


# ===============================
# FLASH CALCULATION LEARNING CONTENT
# ===============================
FLASH_LEARNING = LearningContent(
    background_theory="""
A **flash calculation** determines how a feed stream splits into vapor and liquid phases at 
equilibrium conditions. It's fundamental to understanding separators, flash drums, and distillation.

**The Rachford-Rice Equation:**  
Σ zᵢ(Kᵢ-1) / [1 + (V/F)(Kᵢ-1)] = 0

Where:
- **zᵢ** = Feed composition (mole fraction)
- **Kᵢ** = K-value for component i (= y/x from VLE)
- **V/F** = Vapor fraction (moles vapor per mole feed)

**Key Outputs:**
- V/F: Vapor fraction (0 = all liquid, 1 = all vapor)
- xᵢ: Liquid compositions
- yᵢ: Vapor compositions

Flash drums are used everywhere: crude oil processing, refrigeration, chemical plants.
""",
    key_concepts=[
        "Material balance + VLE = flash calculation",
        "V/F determines how much feed vaporizes",
        "K-values determine where components go",
        "If all Kᵢ < 1: subcooled liquid (no vapor)"
    ],
    real_world_applications=[
        "Flash separator design in oil/gas production",
        "Distillation feed tray conditions",
        "Refrigeration expansion valves",
        "LNG regasification"
    ],
    common_mistakes=[
        "Using K-values at wrong T/P conditions",
        "Not checking if flash is possible (bubble/dew point limits)",
        "Forgetting that K-values change with T and P",
        "Not normalizing compositions to sum = 1"
    ],
    quiz_questions=[
        QuizQuestion(
            id="flash_q1",
            question="If V/F = 0.5, what does this mean?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Half the feed vaporizes", "All becomes vapor", "All remains liquid", "System is at critical point"],
            correct_answer="Half the feed vaporizes",
            explanation="V/F = 0.5 means 50% of feed moles go to vapor phase.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=20,
    prerequisites=["raoults_law", "antoine"],
    related_equations=["raoults_law", "antoine"]
)


class AntoineVaporPressure(BaseEquation):
    """Antoine equation for vapor pressure."""
    
    equation_id = "antoine"
    name = "Antoine Vapor Pressure"
    category = "Thermodynamics"
    description = "Calculate saturation pressure using Antoine equation: log₁₀(P) = A - B/(C+T)"
    reference = "NIST Chemistry WebBook"
    learning_content = ANTOINE_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("A", "Antoine A coefficient - compound-specific constant", 
                              "dimensionless", "", ParameterType.INPUT, symbol="A",
                              tooltip="Look up in Antoine tables (NIST, Perry's)"),
            EquationParameter("B", "Antoine B coefficient", "dimensionless", "", 
                              ParameterType.INPUT, symbol="B"),
            EquationParameter("C", "Antoine C coefficient", "dimensionless", "", 
                              ParameterType.INPUT, symbol="C"),
            EquationParameter("T", "Temperature - must match units used for coefficients", 
                              "temperature", "degC", ParameterType.INPUT, symbol="T",
                              tooltip="Most Antoine tables use °C"),
            EquationParameter("P_vap", "Vapor (saturation) pressure", "pressure", "mmHg", 
                              ParameterType.OUTPUT, symbol="Pᵛᵃᵖ",
                              tooltip="Pressure at which liquid and vapor are in equilibrium"),
        ]
    
    def _calculate(self, A: pint.Quantity, B: pint.Quantity, C: pint.Quantity,
                   T: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        a = A.magnitude if hasattr(A, 'magnitude') else float(A)
        b = B.magnitude if hasattr(B, 'magnitude') else float(B)
        c = C.magnitude if hasattr(C, 'magnitude') else float(C)
        t = T.to('degC').magnitude
        
        log_p = a - b / (c + t)
        p = 10 ** log_p
        
        return {"P_vap": ureg.Quantity(p, "mmHg")}


class RaoultsLaw(BaseEquation):
    """Raoult's law for ideal VLE."""
    
    equation_id = "raoults_law"
    name = "Raoult's Law (Ideal VLE)"
    category = "Thermodynamics"
    description = "Calculate vapor composition for ideal solution: yᵢ·P = xᵢ·Pᵢˢᵃᵗ"
    reference = "Smith, Van Ness & Abbott"
    learning_content = RAOULTS_LAW_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x", "Liquid mole fraction of component", 
                              "dimensionless", "", ParameterType.INPUT, symbol="x", 
                              typical_range=(0, 1),
                              tooltip="Moles of component / total moles in liquid"),
            EquationParameter("P_sat", "Saturation pressure of pure component at T", 
                              "pressure", "psia", ParameterType.INPUT, symbol="Pˢᵃᵗ",
                              tooltip="From Antoine equation or tables"),
            EquationParameter("P_total", "Total system pressure", "pressure", "psia", 
                              ParameterType.INPUT, symbol="P"),
            EquationParameter("y", "Vapor mole fraction of component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y", 
                              typical_range=(0, 1)),
            EquationParameter("K", "K-value (equilibrium ratio) - y/x", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="K",
                              tooltip="K > 1: component concentrates in vapor"),
        ]
    
    def _calculate(self, x: pint.Quantity, P_sat: pint.Quantity,
                   P_total: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        x_val = x.magnitude if hasattr(x, 'magnitude') else float(x)
        psat = P_sat.to('psia').magnitude
        ptot = P_total.to('psia').magnitude
        
        k = psat / ptot
        y = x_val * k
        
        return {
            "y": ureg.Quantity(min(y, 1.0), ""),
            "K": ureg.Quantity(k, "")
        }


class ClausiusClapeyron(BaseEquation):
    """Clausius-Clapeyron equation for vapor pressure vs temperature."""
    
    equation_id = "clausius_clapeyron"
    name = "Clausius-Clapeyron Equation"
    category = "Thermodynamics"
    description = "Estimate vapor pressure at T2 given vapor pressure at T1 and heat of vaporization"
    reference = "Standard thermodynamics"
    learning_content = CLAUSIUS_CLAPEYRON_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Vapor pressure at temperature T1", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₁"),
            EquationParameter("T1", "Reference temperature (absolute)", 
                              "temperature", "degF", ParameterType.INPUT, symbol="T₁"),
            EquationParameter("T2", "Target temperature (absolute)", 
                              "temperature", "degF", ParameterType.INPUT, symbol="T₂"),
            EquationParameter("Hvap", "Heat of vaporization - latent heat", 
                              "energy", "BTU/lb", ParameterType.INPUT, symbol="ΔHᵥₐₚ",
                              tooltip="Water≈970 BTU/lb, varies with temperature"),
            EquationParameter("MW", "Molecular weight of substance", 
                              "dimensionless", "", ParameterType.INPUT, symbol="MW",
                              tooltip="g/mol, e.g., Water=18, Ethanol=46"),
            EquationParameter("P2", "Vapor pressure at temperature T2", 
                              "pressure", "psia", ParameterType.OUTPUT, symbol="P₂"),
        ]
    
    def _calculate(self, P1: pint.Quantity, T1: pint.Quantity, T2: pint.Quantity,
                   Hvap: pint.Quantity, MW: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('Pa').magnitude
        t1 = T1.to('K').magnitude
        t2 = T2.to('K').magnitude
        hvap = Hvap.to('J/kg').magnitude
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        
        R = 8.314  # J/(mol·K)
        hvap_molar = hvap * mw / 1000  # J/mol
        
        ln_p2_p1 = (hvap_molar / R) * (1/t1 - 1/t2)
        p2 = p1 * np.exp(ln_p2_p1)
        
        return {"P2": ureg.Quantity(p2 / 6894.76, "psia")}


class FlashCalculation(BaseEquation):
    """Simple flash drum vapor fraction calculation."""
    
    equation_id = "flash"
    name = "Flash Drum Vapor Fraction"
    category = "Thermodynamics"
    description = "Calculate vapor fraction in a flash drum using Rachford-Rice"
    reference = "Seader & Henley, Separation Process Principles"
    learning_content = FLASH_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("z", "Feed composition (mole fraction of light component)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="z",
                              typical_range=(0, 1)),
            EquationParameter("K1", "K-value of light component at drum conditions", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K₁",
                              tooltip="K = y/x = Psat/P for ideal mixtures"),
            EquationParameter("K2", "K-value of heavy component at drum conditions", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K₂"),
            EquationParameter("V_F", "Vapor fraction - moles vapor per mole feed", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="V/F",
                              typical_range=(0, 1)),
            EquationParameter("y1", "Vapor composition of light component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="y₁"),
            EquationParameter("x1", "Liquid composition of light component", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="x₁"),
        ]
    
    def _calculate(self, z: pint.Quantity, K1: pint.Quantity, K2: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        z1 = z.magnitude if hasattr(z, 'magnitude') else float(z)
        z2 = 1 - z1
        k1 = K1.magnitude if hasattr(K1, 'magnitude') else float(K1)
        k2 = K2.magnitude if hasattr(K2, 'magnitude') else float(K2)
        
        # Rachford-Rice: solve sum(zi(Ki-1)/(1+V/F(Ki-1))) = 0
        # For binary: V/F = (z1(K1-1) + z2(K2-1)) / ((K1-1)(K2-1))... simplified
        # Use bisection for robustness
        def rr(vf):
            return z1*(k1-1)/(1+vf*(k1-1)) + z2*(k2-1)/(1+vf*(k2-1))
        
        # Bisection
        vf_low, vf_high = 0.0, 1.0
        for _ in range(50):
            vf_mid = (vf_low + vf_high) / 2
            if rr(vf_mid) > 0:
                vf_low = vf_mid
            else:
                vf_high = vf_mid
        
        vf = (vf_low + vf_high) / 2
        vf = max(0, min(1, vf))
        
        x1 = z1 / (1 + vf*(k1-1))
        y1 = k1 * x1
        
        return {
            "V_F": ureg.Quantity(vf, ""),
            "y1": ureg.Quantity(y1, ""),
            "x1": ureg.Quantity(x1, "")
        }


class CompressorWork(BaseEquation):
    """Calculate ideal compressor work for gas compression."""
    
    equation_id = "compressor_work"
    name = "Compressor Work (Isentropic)"
    category = "Thermodynamics"
    description = "Calculate ideal isentropic work for gas compression"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("P1", "Suction pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Discharge pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₂"),
            EquationParameter("T1", "Suction temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="T₁"),
            EquationParameter("k", "Ratio of specific heats (Cp/Cv) - heat capacity ratio", 
                              "dimensionless", "", ParameterType.INPUT, symbol="k",
                              typical_range=(1.1, 1.7),
                              tooltip="Air=1.4, Methane=1.31, CO2=1.29"),
            EquationParameter("MW", "Molecular weight of gas", "dimensionless", "", 
                              ParameterType.INPUT, symbol="MW"),
            EquationParameter("m", "Mass flow rate", "mass_flow", "lb/hr", 
                              ParameterType.INPUT, symbol="ṁ"),
            EquationParameter("eta", "Isentropic efficiency", "dimensionless", "", 
                              ParameterType.INPUT, symbol="η", typical_range=(0.7, 0.9)),
            EquationParameter("W_ideal", "Ideal isentropic work", "power", "hp", 
                              ParameterType.OUTPUT, symbol="Wᵢₛ"),
            EquationParameter("W_actual", "Actual shaft work (accounting for efficiency)", 
                              "power", "hp", ParameterType.OUTPUT, symbol="Wₐₓ"),
            EquationParameter("T2_ideal", "Ideal discharge temperature", "temperature", "degF", 
                              ParameterType.OUTPUT, symbol="T₂ᵢₛ"),
        ]
    
    def _calculate(self, P1: pint.Quantity, P2: pint.Quantity, T1: pint.Quantity,
                   k: pint.Quantity, MW: pint.Quantity, m: pint.Quantity,
                   eta: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        p1 = P1.to('psia').magnitude
        p2 = P2.to('psia').magnitude
        t1 = T1.to('degR').magnitude
        k_val = k.magnitude if hasattr(k, 'magnitude') else float(k)
        mw = MW.magnitude if hasattr(MW, 'magnitude') else float(MW)
        m_val = m.to('lb/hr').magnitude
        eff = eta.magnitude if hasattr(eta, 'magnitude') else float(eta)
        
        R = 1545 / mw  # ft·lbf/(lb·°R)
        
        # Isentropic work: W = k/(k-1) * R * T1 * [(P2/P1)^((k-1)/k) - 1]
        ratio = (p2/p1)**((k_val-1)/k_val)
        w_spec = k_val/(k_val-1) * R * t1 * (ratio - 1)  # ft·lbf/lb
        
        w_ideal = m_val * w_spec / (3600 * 550)  # hp
        w_actual = w_ideal / eff
        
        # Discharge temperature
        t2_ideal = t1 * ratio - 459.67  # Convert back to °F
        
        return {
            "W_ideal": ureg.Quantity(w_ideal, "hp"),
            "W_actual": ureg.Quantity(w_actual, "hp"),
            "T2_ideal": ureg.Quantity(t2_ideal, "degF")
        }


THERMODYNAMICS_EQUATIONS = {
    'ideal_gas': IdealGas,
    'antoine': AntoineVaporPressure,
    'raoults_law': RaoultsLaw,
    'clausius_clapeyron': ClausiusClapeyron,
    'flash': FlashCalculation,
    'compressor_work': CompressorWork,
}


class VanDerWaals(BaseEquation):
    """Van der Waals equation of state for real gases."""
    
    equation_id = "van_der_waals"
    name = "Van der Waals Equation"
    category = "Thermodynamics"
    description = "(P + a/V²)(V - b) = RT - real gas behavior"
    reference = "Van der Waals, 1873"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("T", "Temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("P", "Pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P"),
            EquationParameter("Tc", "Critical temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="Tᶜ"),
            EquationParameter("Pc", "Critical pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="Pᶜ"),
            EquationParameter("Z", "Compressibility factor", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Z"),
            EquationParameter("V_molar", "Molar volume", "L/mol", "",
                              ParameterType.OUTPUT, symbol="Vₘ"),
        ]
    
    def _calculate(self, T: pint.Quantity, P: pint.Quantity, Tc: pint.Quantity, 
                   Pc: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        t = T.to('K').magnitude
        p = P.to('Pa').magnitude
        tc = Tc.to('K').magnitude
        pc = Pc.to('Pa').magnitude
        R = 8.314
        
        # Calculate a and b from critical properties
        a = 27 * R**2 * tc**2 / (64 * pc)
        b = R * tc / (8 * pc)
        
        # Solve cubic for V (iterative)
        v = R * t / p  # ideal gas initial guess
        for _ in range(20):
            v_new = (R * t / (p + a/v**2)) + b
            if abs(v_new - v) < 1e-10:
                break
            v = v_new
        
        z = p * v / (R * t)
        v_liters = v * 1000  # L/mol
        
        return {
            "Z": ureg.Quantity(z, ""),
            "V_molar": ureg.Quantity(v_liters, "L/mol")
        }


class ActivityCoefficient(BaseEquation):
    """Margules equation for activity coefficients."""
    
    equation_id = "activity_coefficient"
    name = "Margules Activity Coefficient"
    category = "Thermodynamics"
    description = "Calculate activity coefficients for non-ideal liquid mixtures"
    reference = "Smith, Van Ness & Abbott"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("x1", "Mole fraction of component 1", "dimensionless", "",
                              ParameterType.INPUT, symbol="x₁", typical_range=(0, 1)),
            EquationParameter("A12", "Margules parameter A12", "dimensionless", "",
                              ParameterType.INPUT, symbol="A₁₂", typical_range=(0, 3)),
            EquationParameter("A21", "Margules parameter A21", "dimensionless", "",
                              ParameterType.INPUT, symbol="A₂₁", typical_range=(0, 3)),
            EquationParameter("gamma1", "Activity coefficient of component 1", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="γ₁"),
            EquationParameter("gamma2", "Activity coefficient of component 2", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="γ₂"),
        ]
    
    def _calculate(self, x1: pint.Quantity, A12: pint.Quantity, A21: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        x = x1.magnitude if hasattr(x1, 'magnitude') else float(x1)
        a12 = A12.magnitude if hasattr(A12, 'magnitude') else float(A12)
        a21 = A21.magnitude if hasattr(A21, 'magnitude') else float(A21)
        
        x2 = 1 - x
        ln_gamma1 = x2**2 * (a12 + 2*(a21 - a12)*x)
        ln_gamma2 = x**2 * (a21 + 2*(a12 - a21)*x2)
        
        return {
            "gamma1": ureg.Quantity(np.exp(ln_gamma1), ""),
            "gamma2": ureg.Quantity(np.exp(ln_gamma2), "")
        }


class JouleThomson(BaseEquation):
    """Joule-Thomson coefficient for throttling calculations."""
    
    equation_id = "joule_thomson"
    name = "Joule-Thomson Effect"
    category = "Thermodynamics"
    description = "Temperature change during isenthalpic throttling"
    reference = "Perry's Chemical Engineers' Handbook"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("mu_JT", "Joule-Thomson coefficient", "degC/bar", "degF/psi",
                              ParameterType.INPUT, symbol="μⱼₜ", typical_range=(-1, 1)),
            EquationParameter("P1", "Inlet pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P₁"),
            EquationParameter("P2", "Outlet pressure", "pressure", "bar",
                              ParameterType.INPUT, symbol="P₂"),
            EquationParameter("T1", "Inlet temperature", "temperature", "degC",
                              ParameterType.INPUT, symbol="T₁"),
            EquationParameter("delta_T", "Temperature change", "delta_degC", "delta_degF",
                              ParameterType.OUTPUT, symbol="ΔT"),
            EquationParameter("T2", "Outlet temperature", "temperature", "degC",
                              ParameterType.OUTPUT, symbol="T₂"),
        ]
    
    def _calculate(self, mu_JT: pint.Quantity, P1: pint.Quantity, P2: pint.Quantity,
                   T1: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        mu = mu_JT.magnitude if hasattr(mu_JT, 'magnitude') else float(mu_JT)
        p1 = P1.to('bar').magnitude
        p2 = P2.to('bar').magnitude
        t1 = T1.to('degC').magnitude
        
        delta_t = mu * (p2 - p1)  # Usually negative (cooling)
        t2 = t1 + delta_t
        
        return {
            "delta_T": ureg.Quantity(delta_t, "delta_degC"),
            "T2": ureg.Quantity(t2, "degC")
        }


class HeatCapacityMixture(BaseEquation):
    """Heat capacity of ideal gas mixture."""
    
    equation_id = "heat_capacity_mixture"
    name = "Mixture Heat Capacity"
    category = "Thermodynamics"
    description = "Calculate Cp of ideal gas mixture from pure component values"
    reference = "Standard thermodynamics"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("y1", "Mole fraction of component 1", "dimensionless", "",
                              ParameterType.INPUT, symbol="y₁", typical_range=(0, 1)),
            EquationParameter("Cp1", "Heat capacity of component 1", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.INPUT, symbol="Cₚ₁"),
            EquationParameter("Cp2", "Heat capacity of component 2", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.INPUT, symbol="Cₚ₂"),
            EquationParameter("Cp_mix", "Mixture heat capacity", "J/(mol*K)", "BTU/(lbmol*R)",
                              ParameterType.OUTPUT, symbol="Cₚₘᵢₓ"),
        ]
    
    def _calculate(self, y1: pint.Quantity, Cp1: pint.Quantity, Cp2: pint.Quantity,
                   **kwargs) -> Dict[str, pint.Quantity]:
        y = y1.magnitude if hasattr(y1, 'magnitude') else float(y1)
        cp1 = Cp1.to('J/(mol*K)').magnitude
        cp2 = Cp2.to('J/(mol*K)').magnitude
        
        cp_mix = y * cp1 + (1-y) * cp2
        
        return {"Cp_mix": ureg.Quantity(cp_mix, "J/(mol*K)")}


# Update registry
THERMODYNAMICS_EQUATIONS = {
    'ideal_gas': IdealGas,
    'antoine': AntoineVaporPressure,
    'raoults_law': RaoultsLaw,
    'clausius_clapeyron': ClausiusClapeyron,
    'flash': FlashCalculation,
    'compressor_work': CompressorWork,
    'van_der_waals': VanDerWaals,
    'activity_coefficient': ActivityCoefficient,
    'joule_thomson': JouleThomson,
    'heat_capacity_mixture': HeatCapacityMixture,
}
