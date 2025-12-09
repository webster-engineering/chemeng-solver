"""
Process Control equations - PID tuning, control valve sizing, response analysis.
This is the primary focus module for the solver.
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
# ZIEGLER-NICHOLS LEARNING CONTENT
# ===============================
ZIEGLER_NICHOLS_LEARNING = LearningContent(
    background_theory="""
The **Ziegler-Nichols method** (1942) is the classic PID tuning technique. It uses the 
**ultimate gain (Ku)** and **ultimate period (Pu)** found by pushing the system to oscillation.

**Procedure:**
1. Set Ti = ∞, Td = 0 (P-only control)
2. Increase Kp until sustained oscillation occurs
3. Record: Ku = that Kp value, Pu = oscillation period

**Tuning Rules:**
| Controller | Kp | Ti | Td |
|---|---|---|---|
| P | 0.5·Ku | ∞ | 0 |
| PI | 0.45·Ku | Pu/1.2 | 0 |
| PID | 0.6·Ku | Pu/2 | Pu/8 |

**Characteristics:**
- Aggressive tuning - fast response but ~25% overshoot
- Good for set point changes, poor for disturbance rejection
- Quarter-decay ratio target
""",
    key_concepts=[
        "Ultimate gain Ku = controller gain at sustained oscillation",
        "Ultimate period Pu = time between peaks at Ku",
        "Produces aggressive tuning with ~25% overshoot",
        "Tyreus-Luyben is more conservative alternative"
    ],
    real_world_applications=[
        "Initial tuning starting point",
        "Flow and pressure loops (fast dynamics)",
        "Educational demonstration of feedback control",
        "Loops where speed matters more than overshoot"
    ],
    common_mistakes=[
        "Using this on integrating processes (level control)",
        "Expecting slug response - Z-N gives oscillatory response",
        "Not recognizing when process won't oscillate (overdamped)",
        "Applying to loops with significant noise"
    ],
    pro_tips=[
        "For less aggressive tuning, use Tyreus-Luyben (Kp=Ku/3.2)",
        "If you can't find Ku, use Cohen-Coon or Lambda tuning",
        "Many DCS systems have auto-tune that does this automatically"
    ],
    quiz_questions=[
        QuizQuestion(
            id="zn_q1",
            question="If Ku = 4 and Pu = 2 min, what is Kp for Z-N PID tuning?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="2.4",
            explanation="Kp = 0.6 × Ku = 0.6 × 4 = 2.4",
            difficulty=Difficulty.BEGINNER,
            hint="Kp = 0.6 × Ku for PID",
            points=10
        ),
        QuizQuestion(
            id="zn_q2",
            question="Ziegler-Nichols tuning typically produces what approximate overshoot?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["0%", "10%", "25%", "50%"],
            correct_answer="25%",
            explanation="Z-N targets quarter-decay ratio, which corresponds to ~25% overshoot.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=[],
    related_equations=["tyreus_luyben_pid", "lambda_tuning"],
    references=[
        "Ziegler, J.G. & Nichols, N.B. 'Optimum Settings for Automatic Controllers', Trans. ASME, 1942",
        "Åström, K.J. & Hägglund, T. 'PID Controllers: Theory, Design, and Tuning', ISA, 1995"
    ],
    derivation_summary="Empirical method developed at Taylor Instruments. Based on finding the point of marginal stability (sustained oscillation) and backing off from there using experimentally derived ratios.",
    limitations_assumptions=[
        "Results in ~25% overshoot - may be too aggressive for some processes",
        "Not suitable for integrating processes (level control)",
        "Process must be able to oscillate under P-only control",
        "Noise can make it difficult to identify true Ku and Pu"
    ],
    diagram_type="control_loop"
)


# ===============================
# LAMBDA TUNING LEARNING CONTENT
# ===============================
LAMBDA_TUNING_LEARNING = LearningContent(
    background_theory="""
**Lambda (λ) Tuning** is the most popular industrial PID tuning method. It gives you direct 
control over the closed-loop speed via the λ parameter.

**The Idea:** Choose your desired closed-loop time constant λ, then calculate controller gains.

**For First-Order Plus Dead Time (FOPDT) Process:**
- **Kc = τ / [K × (λ + θ)]**
- **Ti = τ** (integral time = process time constant)
- **Td = 0** (usually PI, or θ/2 for PID)

**Choosing λ:**
| λ Value | Response |
|---------|----------|
| λ = θ | Fast but aggressive |
| λ = τ | Moderate (good starting point) |
| λ = 3θ | Conservative, robust |

**Why It's Popular:**
- Single tuning parameter (λ)
- Directly controls speed vs robustness trade-off
- Works with process identification
""",
    key_concepts=[
        "λ = desired closed-loop time constant",
        "Larger λ = slower but more robust response",
        "Requires FOPDT model: K, τ, θ from step test",
        "Industry-standard method for most loops"
    ],
    real_world_applications=[
        "Temperature loop tuning (slow processes)",
        "Composition/analyzer loops",
        "Any loop where you have step test data",
        "Model predictive control initialization"
    ],
    common_mistakes=[
        "Using λ < θ (aggressive, can cause instability)",
        "Poor FOPDT model identification",
        "Forgetting to account for dead time θ"
    ],
    quiz_questions=[
        QuizQuestion(
            id="lam_q1",
            question="For a process with K=2, τ=10 min, θ=2 min, and λ=4 min, what is Kc?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="0.833",
            explanation="Kc = τ / [K(λ+θ)] = 10 / [2×(4+2)] = 10/12 = 0.833",
            difficulty=Difficulty.INTERMEDIATE,
            hint="Kc = τ / [K(λ+θ)]",
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=[],
    related_equations=["ziegler_nichols_pid", "simc_tuning"],
    references=[
        "Rivera, D.E., Morari, M. & Skogestad, S. 'Internal Model Control', Ind. Eng. Chem. Process Des. Dev., 1986",
        "Skogestad, S. 'Simple analytic rules for model reduction and PID controller tuning', J. Process Control, 2003"
    ],
    derivation_summary="Based on Internal Model Control theory. The tuning parameter λ directly sets the desired closed-loop time constant, providing intuitive control over the speed-robustness tradeoff.",
    limitations_assumptions=[
        "Requires FOPDT model identification from step test",
        "λ should generally not be smaller than dead time θ",
        "Performance depends on model accuracy",
        "Original IMC formulation assumes no model mismatch"
    ]
)


# ===============================
# CONTROL VALVE CV LEARNING CONTENT
# ===============================
CV_LIQUID_LEARNING = LearningContent(
    background_theory="""
**Cv (Valve Coefficient)** is THE valve sizing parameter. It defines how much 
flow a valve passes for a given pressure drop.

**Definition:** Cv = flow rate (gpm) that passes through the valve with 1 psi 
pressure drop and specific gravity = 1.0

**The Equation (Liquid):** Cv = Q × √(SG/ΔP)

Where:
- **Q** = Flow rate (gpm)
- **SG** = Specific gravity (water = 1.0)
- **ΔP** = Pressure drop across valve (psi)

**Valve Selection:**
1. Calculate required Cv at design conditions
2. Add margin (typically select valve where design Cv is ~70-80% of valve max Cv)
3. Verify valve can handle turndown
""",
    key_concepts=[
        "Cv = gpm at SG=1 and ΔP=1 psi",
        "Higher Cv = larger valve capacity",
        "Always add margin - don't size at 100% Cv",
        "Rangeability = Cv_max / Cv_min (typically 50:1)"
    ],
    real_world_applications=[
        "Control valve sizing for new installations",
        "Checking existing valve adequacy",
        "Troubleshooting valve capacity issues",
        "System hydraulic analysis"
    ],
    common_mistakes=[
        "Using gauge pressure instead of differential ΔP",
        "Forgetting specific gravity correction",
        "Not checking for choked flow conditions",
        "Sizing at 100% Cv (no margin for upset)"
    ],
    quiz_questions=[
        QuizQuestion(
            id="cv_q1",
            question="100 gpm water flow with 25 psi pressure drop. What Cv is needed?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="20",
            explanation="Cv = Q × √(SG/ΔP) = 100 × √(1/25) = 100 × 0.2 = 20",
            difficulty=Difficulty.BEGINNER,
            hint="Cv = Q × √(SG/ΔP)",
            points=10
        )
    ],
    difficulty=Difficulty.BEGINNER,
    estimated_time_minutes=12,
    prerequisites=[],
    related_equations=["cv_gas", "cv_travel"],
    references=[
        "ISA-75.01.01 - Flow Equations for Sizing Control Valves",
        "Masoneilan Control Valve Handbook",
        "Fisher Control Valve Handbook, Emerson"
    ],
    derivation_summary="Cv is defined empirically as the flow rate (gpm) of water at 60°F that will flow through a valve with a 1 psi pressure drop. The equation derives from Bernoulli with empirical corrections.",
    limitations_assumptions=[
        "Valid for non-flashing, non-cavitating liquid service",
        "Assumes fully turbulent flow through valve",
        "Does not account for choked flow conditions",
        "Pipe reducer corrections may be needed for high velocities"
    ]
)


# ===============================
# COHEN-COON LEARNING CONTENT
# ===============================
COHEN_COON_LEARNING = LearningContent(
    background_theory="""
**Cohen-Coon tuning** uses the process reaction curve (step test) to identify 
a First-Order Plus Dead Time (FOPDT) model, then applies tuning rules.

**The FOPDT Model:** G(s) = K × e^(-θs) / (τs + 1)

From step test, identify:
- **K** = Process gain (final value / step size)
- **τ** = Time constant (63.2% of final value)
- **θ** = Dead time (delay before response starts)

**Cohen-Coon PID Formulas:**
- Kp = (τ/Kθ) × (4/3 + θ/4τ)
- Ti = θ × (32 + 6θ/τ) / (13 + 8θ/τ)
- Td = 4θ / (11 + 2θ/τ)

Works well for processes with θ/τ < 1 (moderate dead time).
""",
    key_concepts=[
        "Uses process reaction curve (step test) data",
        "Requires FOPDT model: K, τ, θ",
        "Works well for θ/τ < 1",
        "More conservative than Z-N for many processes"
    ],
    common_mistakes=[
        "Poor model identification from noisy step test",
        "Applying when θ/τ > 1 (high dead time ratio)",
        "Confusing effective dead time with transportation lag"
    ],
    quiz_questions=[
        QuizQuestion(
            id="cc_q1",
            question="What does K represent in FOPDT model?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Process gain", "Controller gain", "Time constant", "Dead time"],
            correct_answer="Process gain",
            explanation="K = process gain = (final output change)/(input step size) at steady state.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=[],
    related_equations=["ziegler_nichols_pid", "lambda_tuning"],
    references=[
        "Cohen, G.H. & Coon, G.A. 'Theoretical Consideration of Retarded Control', Trans. ASME, 1953",
        "Ogunnaike, B.A. & Ray, W.H. 'Process Dynamics, Modeling, and Control', Oxford, 1994"
    ],
    derivation_summary="Uses process reaction curve to fit FOPDT model parameters K, τ, θ. Tuning formulas derived empirically to achieve quarter-decay ratio response with improved disturbance rejection.",
    limitations_assumptions=[
        "Best for processes with θ/τ < 1 (moderate dead time)",
        "Requires clean step test for model identification",
        "May give aggressive tuning for large dead time processes",
        "Assumes FOPDT adequately represents process dynamics"
    ]
)


# ===============================
# IMC PID LEARNING CONTENT
# ===============================
IMC_LEARNING = LearningContent(
    background_theory="""
**Internal Model Control (IMC)** PID tuning is similar to Lambda tuning.
It uses a process model and a single tuning parameter for desired closed-loop response.

**IMC for FOPDT Process:**
- Kc = τ / [K × (λ + θ)]
- Ti = τ
- Td = θ/2 (for PID) or 0 (for PI)

**Choosing λ (filter time constant):**
- λ = θ → Aggressive
- λ = τ → Moderate
- λ = max(0.1τ, 0.8θ) → Skogestad rule

IMC provides guaranteed stability if model is accurate.
""",
    key_concepts=[
        "Model-based tuning approach",
        "λ controls aggressiveness vs robustness",
        "Requires good FOPDT model",
        "Essentially same as Lambda tuning"
    ],
    common_mistakes=[
        "Using too small λ (aggressive, oscillatory)",
        "Model mismatch causing poor performance"
    ],
    quiz_questions=[
        QuizQuestion(
            id="imc_q1",
            question="In IMC tuning, what does larger λ give?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["Slower, more robust response", "Faster, more aggressive response", "Higher overshoot", "Lower integral time"],
            correct_answer="Slower, more robust response",
            explanation="Larger λ = larger closed-loop time constant = slower but more robust response.",
            difficulty=Difficulty.BEGINNER,
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=12,
    prerequisites=[],
    related_equations=["lambda_tuning", "ziegler_nichols_pid"],
    references=[
        "Morari, M. & Zafiriou, E. 'Robust Process Control', Prentice Hall, 1989",
        "Skogestad, S. 'Simple analytic rules for model reduction and PID controller tuning', J. Process Control, 2003"
    ],
    derivation_summary="IMC derives controller settings by inverting the process model and adding a filter. For FOPDT, the resulting PI/PID settings are equivalent to Lambda tuning.",
    limitations_assumptions=[
        "Assumes process model is reasonably accurate",
        "Performance degrades with model mismatch",
        "Filter time constant λ must be chosen appropriately",
        "Not suitable for unstable processes without modification"
    ]
)


# ===============================
# TYREUS-LUYBEN LEARNING CONTENT
# ===============================
TYREUS_LUYBEN_LEARNING = LearningContent(
    background_theory="""
**Tyreus-Luyben tuning** is a more conservative alternative to Ziegler-Nichols.
It uses the same Ku and Pu but gives less aggressive settings.

**Tyreus-Luyben Formulas:**
- Kp = Ku / 3.2
- Ti = 2.2 × Pu
- Td = Pu / 6.3

**Comparison to Z-N:**
| Parameter | Z-N | T-L |
|-----------|-----|-----|
| Kp | 0.6 Ku | 0.31 Ku |
| Ti | 0.5 Pu | 2.2 Pu |

T-L gives less overshoot and is preferred for chemical process control.
""",
    key_concepts=[
        "Uses same Ku, Pu as Ziegler-Nichols",
        "Kp ≈ half of Z-N (less aggressive)",
        "Higher integral time = slower reset",
        "Better for process control (less overshoot)"
    ],
    common_mistakes=[
        "Expecting same fast response as Z-N",
        "Confusing with Z-N formulas"
    ],
    quiz_questions=[
        QuizQuestion(
            id="tl_q1",
            question="If Ku = 5, what is Kp using Tyreus-Luyben?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="1.56",
            explanation="Kp = Ku/3.2 = 5/3.2 = 1.5625 ≈ 1.56",
            difficulty=Difficulty.BEGINNER,
            hint="Kp = Ku / 3.2",
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=10,
    prerequisites=[],
    related_equations=["ziegler_nichols_pid"],
    references=[
        "Tyreus, B.D. & Luyben, W.L. 'Tuning PI Controllers for Integrator/Dead Time Processes', Ind. Eng. Chem. Res., 1992",
        "Luyben, W.L. 'Essentials of Process Control', McGraw-Hill, 1997"
    ],
    derivation_summary="Uses same ultimate gain/period as Z-N but with more conservative ratios. Developed specifically for chemical process control where less overshoot is preferred.",
    limitations_assumptions=[
        "Slower response than Z-N - not suitable where speed is critical",
        "Still requires ability to find Ku and Pu",
        "Better for disturbance rejection than setpoint tracking"
    ]
)


# ===============================
# CV GAS LEARNING CONTENT
# ===============================
CV_GAS_LEARNING = LearningContent(
    background_theory="""
**Control valve sizing for gas/vapor** is more complex than liquids because
gas expands as pressure drops (compressible flow).

**Key Equation:**
Cv = Q / [963 × P1 × (1 - x/3) × √(x / (SG × T))]

Where:
- Q = Standard volumetric flow (scfm at 14.7 psia, 60°F)
- P1 = Upstream absolute pressure (psia)
- x = ΔP/P1 = Pressure drop ratio
- SG = Specific gravity (vs air = 1.0)
- T = Temperature (°R)

**Choked Flow:**
When x > 0.5, flow is choked (sonic velocity at vena contracta).
Flow cannot increase even with more ΔP.
""",
    key_concepts=[
        "Gas flow is compressible - different from liquid",
        "Choked flow when ΔP/P1 > 0.5",
        "Use absolute pressure and temperature",
        "SG relative to air = 1.0"
    ],
    common_mistakes=[
        "Using liquid Cv equation for gas",
        "Using gauge instead of absolute pressure",
        "Not checking for choked flow conditions"
    ],
    quiz_questions=[
        QuizQuestion(
            id="cvg_q1",
            question="At what pressure ratio is gas flow typically choked?",
            question_type=QuestionType.MULTIPLE_CHOICE,
            options=["ΔP/P1 > 0.5", "ΔP/P1 > 0.1", "ΔP/P1 > 0.9", "Never chokes"],
            correct_answer="ΔP/P1 > 0.5",
            explanation="When pressure ratio exceeds ~0.5, critical flow occurs and velocity reaches sonic.",
            difficulty=Difficulty.INTERMEDIATE,
            points=15
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=15,
    prerequisites=["cv_liquid"],
    related_equations=["cv_liquid"],
    references=[
        "ISA-75.01.01 - Flow Equations for Sizing Control Valves",
        "Fisher Controls Handbook, 'Sizing Control Valves for Gas and Vapor Flow'",
        "IEC 60534 - Industrial Process Control Valves"
    ],
    derivation_summary="Derived from compressible flow equations with the Cv coefficient. The (1-x/3) term approximates the expansion factor Y for subcritical flow.",
    limitations_assumptions=[
        "Flow is choked when x = ΔP/P1 > 0.5",
        "Use absolute pressure and temperature (psia, °R)",
        "Specific gravity is relative to air (MW/29)",
        "Simplified equation - full IEC 60534 more accurate for critical sizing"
    ]
)


# ===============================
# GAIN MARGIN LEARNING CONTENT
# ===============================
GAIN_MARGIN_LEARNING = LearningContent(
    background_theory="""
**Gain Margin (GM)** tells you how much the loop gain can increase before
the system becomes unstable.

**The Equation:** GM = Ku / Kc

Where:
- Ku = Ultimate gain (controller gain at sustained oscillation)
- Kc = Current controller gain

**Guidelines:**
- GM > 2 (linear) or > 6 dB is recommended
- GM = 1 means at stability boundary (oscillating)
- GM < 1 means unstable!

**In Decibels:** GM(dB) = 20 × log₁₀(GM)
""",
    key_concepts=[
        "GM = Ku/Kc = safety margin to instability",
        "GM > 2 (or > 6 dB) recommended",
        "GM = 1 means at stability boundary",
        "Related to robustness against model uncertainty"
    ],
    common_mistakes=[
        "Confusing gain margin with phase margin",
        "Operating with GM < 1.5 (too close to instability)",
        "Not accounting for process gain changes"
    ],
    quiz_questions=[
        QuizQuestion(
            id="gm_q1",
            question="If Ku = 10 and Kc = 2, what is the gain margin?",
            question_type=QuestionType.NUMERIC,
            options=[],
            correct_answer="5",
            explanation="GM = Ku/Kc = 10/2 = 5",
            difficulty=Difficulty.BEGINNER,
            hint="GM = Ku / Kc",
            points=10
        )
    ],
    difficulty=Difficulty.INTERMEDIATE,
    estimated_time_minutes=10,
    prerequisites=["ziegler_nichols_pid"],
    related_equations=["ziegler_nichols_pid", "tyreus_luyben_pid"],
    references=[
        "Åström, K.J. & Murray, R.M. 'Feedback Systems: An Introduction for Scientists and Engineers', Princeton, 2008",
        "Seborg, Edgar, Mellichamp & Doyle 'Process Dynamics and Control', Wiley, 2016"
    ],
    derivation_summary="Gain margin is the factor by which loop gain can increase before reaching the stability boundary (where Nyquist plot crosses -1). Derived from Nyquist stability criterion.",
    limitations_assumptions=[
        "Assumes linear system (gain margin may vary with operating point)",
        "Should check both gain AND phase margin for robustness",
        "Process gain changes with operating conditions may reduce effective GM",
        "GM > 2 (6 dB) is a guideline, not absolute rule"
    ]
)

class ZieglerNicholsPID(BaseEquation):
    """Ziegler-Nichols ultimate gain method for PID tuning."""
    
    equation_id = "ziegler_nichols_pid"
    name = "Ziegler-Nichols PID Tuning"
    category = "Process Control"
    description = "Calculate PID parameters using ultimate gain (Ku) and ultimate period (Pu)"
    reference = "Ziegler & Nichols, 1942"
    learning_content = ZIEGLER_NICHOLS_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Ku", "Ultimate gain - controller gain at sustained oscillation", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Kᵤ", 
                              typical_range=(0.1, 100),
                              tooltip="Found by increasing Kp until sustained oscillation occurs"),
            EquationParameter("Pu", "Ultimate period - oscillation period at ultimate gain", 
                              "time", "min", ParameterType.INPUT, symbol="Pᵤ", 
                              typical_range=(0.1, 1000),
                              tooltip="Measure peak-to-peak time at sustained oscillation"),
            EquationParameter("controller_type", "Controller type (1=P, 2=PI, 3=PID)", 
                              "dimensionless", "", ParameterType.INPUT, required=False,
                              tooltip="Enter 1 for P, 2 for PI, 3 for PID (default)"),
            EquationParameter("Kp", "Proportional gain - controller P term", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Kₚ"),
            EquationParameter("Ti", "Integral time - time for I action to equal P action", 
                              "time", "min", ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time - time for D action based on rate of change", 
                              "time", "min", ParameterType.OUTPUT, symbol="Tᵈ"),
        ]
    
    def _calculate(self, Ku: pint.Quantity, Pu: pint.Quantity,
                   controller_type: pint.Quantity = None, **kwargs) -> Dict[str, pint.Quantity]:
        ku_val = Ku.magnitude if hasattr(Ku, 'magnitude') else float(Ku)
        pu_val = Pu.to('min').magnitude
        ctrl = controller_type.magnitude if controller_type and hasattr(controller_type, 'magnitude') else 3
        
        if ctrl >= 3:  # PID
            kp = 0.6 * ku_val
            ti = 0.5 * pu_val
            td = 0.125 * pu_val
        elif ctrl >= 2:  # PI
            kp = 0.45 * ku_val
            ti = pu_val / 1.2
            td = 0.0
        else:  # P only
            kp = 0.5 * ku_val
            ti = float('inf')
            td = 0.0
        
        return {
            "Kp": ureg.Quantity(kp, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min")
        }


class CohenCoonPID(BaseEquation):
    """Cohen-Coon method for PID tuning using process reaction curve."""
    
    equation_id = "cohen_coon_pid"
    name = "Cohen-Coon PID Tuning"
    category = "Process Control"
    description = "PID tuning from first-order plus dead time (FOPDT) model parameters"
    reference = "Cohen & Coon, 1953"
    learning_content = COHEN_COON_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain - output change / input change at steady state", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K",
                              tooltip="Also called static gain or DC gain"),
            EquationParameter("tau", "Time constant - time to reach 63.2% of final value", 
                              "time", "min", ParameterType.INPUT, symbol="τ", 
                              typical_range=(0.1, 1000),
                              tooltip="From process reaction curve"),
            EquationParameter("theta", "Dead time - time delay before response begins", 
                              "time", "min", ParameterType.INPUT, symbol="θ", 
                              typical_range=(0, 100),
                              tooltip="Also called lag time or transport delay"),
            EquationParameter("Kp", "Proportional gain", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Kₚ"),
            EquationParameter("Ti", "Integral time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵈ"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau: pint.Quantity,
                   theta: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        k_val = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau_val = tau.to('min').magnitude
        theta_val = theta.to('min').magnitude
        
        if theta_val == 0:
            theta_val = 0.001
        
        r = theta_val / tau_val
        
        kp = (1/k_val) * (tau_val/theta_val) * (4/3 + r/4)
        ti = theta_val * (32 + 6*r) / (13 + 8*r)
        td = theta_val * 4 / (11 + 2*r)
        
        return {
            "Kp": ureg.Quantity(kp, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min")
        }


class IMCPID(BaseEquation):
    """Internal Model Control (IMC) based PID tuning."""
    
    equation_id = "imc_pid"
    name = "IMC-Based PID Tuning"
    category = "Process Control"
    description = "PID tuning with adjustable aggressiveness via λ (lambda) parameter"
    reference = "Rivera, Morari & Skogestad, 1986"
    learning_content = IMC_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain", "dimensionless", "",
                              ParameterType.INPUT, symbol="K"),
            EquationParameter("tau", "Process time constant", "time", "min",
                              ParameterType.INPUT, symbol="τ"),
            EquationParameter("theta", "Dead time (delay)", "time", "min",
                              ParameterType.INPUT, symbol="θ"),
            EquationParameter("lambda_", "IMC filter time constant - tuning parameter", 
                              "time", "min", ParameterType.INPUT, symbol="λ",
                              tooltip="Larger λ = more robust/slower, smaller = more aggressive"),
            EquationParameter("Kp", "Proportional gain", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Kₚ"),
            EquationParameter("Ti", "Integral time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵈ"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau: pint.Quantity, theta: pint.Quantity,
                   lambda_: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        k = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau_val = tau.to('min').magnitude
        theta_val = theta.to('min').magnitude
        lam = lambda_.to('min').magnitude
        
        kp = tau_val / (k * (lam + theta_val))
        ti = tau_val
        td = theta_val / 2
        
        return {
            "Kp": ureg.Quantity(kp, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min")
        }


class TyreusLuybenPID(BaseEquation):
    """Tyreus-Luyben PID tuning - more conservative than Z-N."""
    
    equation_id = "tyreus_luyben_pid"
    name = "Tyreus-Luyben PID Tuning"
    category = "Process Control"
    description = "Conservative PID tuning from ultimate gain/period - less aggressive than Z-N"
    reference = "Tyreus & Luyben, 1992"
    learning_content = TYREUS_LUYBEN_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Ku", "Ultimate gain", "dimensionless", "",
                              ParameterType.INPUT, symbol="Kᵤ",
                              tooltip="Controller gain at sustained oscillation"),
            EquationParameter("Pu", "Ultimate period", "time", "min",
                              ParameterType.INPUT, symbol="Pᵤ"),
            EquationParameter("Kp", "Proportional gain", "dimensionless", "",
                              ParameterType.OUTPUT, symbol="Kₚ"),
            EquationParameter("Ti", "Integral time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time", "time", "min",
                              ParameterType.OUTPUT, symbol="Tᵈ"),
        ]
    
    def _calculate(self, Ku: pint.Quantity, Pu: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        ku_val = Ku.magnitude if hasattr(Ku, 'magnitude') else float(Ku)
        pu_val = Pu.to('min').magnitude
        
        kp = ku_val / 3.2
        ti = 2.2 * pu_val
        td = pu_val / 6.3
        
        return {
            "Kp": ureg.Quantity(kp, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min")
        }


class ControlValveCvLiquid(BaseEquation):
    """Control valve Cv sizing for liquid service."""
    
    equation_id = "cv_liquid"
    name = "Control Valve Cv (Liquid)"
    category = "Process Control"
    description = "Calculate valve coefficient Cv for incompressible (liquid) flow"
    reference = "ISA-75.01.01"
    learning_content = CV_LIQUID_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Volumetric flow rate through valve", 
                              "volumetric_flow", "gpm", ParameterType.INPUT, symbol="Q", 
                              typical_range=(0.1, 100000)),
            EquationParameter("dP", "Pressure drop across valve", 
                              "pressure", "psi", ParameterType.INPUT, symbol="ΔP", 
                              typical_range=(1, 1000),
                              tooltip="Upstream pressure - downstream pressure"),
            EquationParameter("SG", "Specific gravity (relative to water at 60°F)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="SG", 
                              typical_range=(0.5, 2.0),
                              tooltip="Water=1.0, Gasoline≈0.72, Brine≈1.2"),
            EquationParameter("Cv", "Valve flow coefficient", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Cᵥ", 
                              typical_range=(0.1, 10000),
                              tooltip="Flow rate (gpm) with 1 psi drop and SG=1"),
        ]
    
    def _calculate(self, Q: pint.Quantity, dP: pint.Quantity,
                   SG: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gpm').magnitude
        dp = dP.to('psi').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        
        if dp <= 0:
            raise ValueError("Pressure drop must be positive")
        
        cv = q * np.sqrt(sg / dp)
        
        return {"Cv": ureg.Quantity(cv, "")}


class ControlValveCvGas(BaseEquation):
    """Control valve Cv sizing for gas/vapor service."""
    
    equation_id = "cv_gas"
    name = "Control Valve Cv (Gas)"
    category = "Process Control"
    description = "Calculate valve Cv for compressible (gas/vapor) flow"
    reference = "ISA-75.01.01"
    learning_content = CV_GAS_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Standard volumetric flow at 14.7 psia and 60°F", 
                              "volumetric_flow", "scfm", ParameterType.INPUT, symbol="Qₛ"),
            EquationParameter("P1", "Inlet pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₁",
                              tooltip="Upstream of valve; use gauge + 14.7 to get absolute"),
            EquationParameter("P2", "Outlet pressure (absolute)", 
                              "pressure", "psia", ParameterType.INPUT, symbol="P₂"),
            EquationParameter("T", "Flowing temperature", "temperature", "degF", 
                              ParameterType.INPUT, symbol="T"),
            EquationParameter("SG", "Specific gravity (relative to air = 1.0)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="SG",
                              tooltip="Natural gas≈0.6, CO2≈1.52, Propane≈1.52"),
            EquationParameter("Cv", "Valve flow coefficient", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Cᵥ"),
            EquationParameter("choked", "Is flow choked? (1=yes, 0=no)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Choked",
                              tooltip="Choked when ΔP/P1 > 0.5"),
        ]
    
    def _calculate(self, Q: pint.Quantity, P1: pint.Quantity, P2: pint.Quantity,
                   T: pint.Quantity, SG: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('scfm').magnitude
        p1 = P1.to('psia').magnitude
        p2 = P2.to('psia').magnitude
        t = T.to('degR').magnitude
        sg = SG.magnitude if hasattr(SG, 'magnitude') else float(SG)
        
        dp = p1 - p2
        if dp <= 0:
            raise ValueError("P1 must be greater than P2")
        
        x = dp / p1
        choked = 0
        if x > 0.5:  # Choked flow
            x = 0.5
            choked = 1
        
        cv = q / (963 * p1 * (1 - x/3) * np.sqrt(x / (sg * t)))
        
        return {
            "Cv": ureg.Quantity(cv, ""),
            "choked": ureg.Quantity(choked, "")
        }


class ControlLoopGainMargin(BaseEquation):
    """Calculate gain margin for control loop stability analysis."""
    
    equation_id = "gain_margin"
    name = "Gain Margin Calculation"
    category = "Process Control"
    description = "Calculate how much loop gain can increase before instability"
    reference = "Standard Control Theory"
    learning_content = GAIN_MARGIN_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Ku", "Ultimate gain - gain at which oscillation starts", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Kᵤ"),
            EquationParameter("Kc", "Current controller gain", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Kc",
                              tooltip="Currently tuned proportional gain"),
            EquationParameter("GM", "Gain margin - Ku/Kc", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="GM",
                              tooltip="Typically want GM > 2 for stability"),
            EquationParameter("GM_dB", "Gain margin in decibels", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="GM(dB)",
                              tooltip="20·log₁₀(GM), typically want > 6 dB"),
        ]
    
    def _calculate(self, Ku: pint.Quantity, Kc: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        ku = Ku.magnitude if hasattr(Ku, 'magnitude') else float(Ku)
        kc = Kc.magnitude if hasattr(Kc, 'magnitude') else float(Kc)
        
        gm = ku / kc
        gm_db = 20 * np.log10(gm)
        
        return {
            "GM": ureg.Quantity(gm, ""),
            "GM_dB": ureg.Quantity(gm_db, "")
        }


class FirstOrderResponse(BaseEquation):
    """First-order system step response analysis."""
    
    equation_id = "first_order_response"
    name = "First-Order Step Response"
    category = "Process Control"
    description = "Calculate response times for a first-order (single time constant) system"
    reference = "Seborg et al., Process Dynamics and Control"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain - steady-state output change per unit input", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K"),
            EquationParameter("tau", "Time constant - time to reach 63.2% of final change", 
                              "time", "min", ParameterType.INPUT, symbol="τ",
                              tooltip="Larger τ = slower response"),
            EquationParameter("theta", "Dead time - delay before any response occurs", 
                              "time", "min", ParameterType.INPUT, symbol="θ", required=False),
            EquationParameter("step_size", "Step input magnitude - size of step change", 
                              "dimensionless", "", ParameterType.INPUT, symbol="M"),
            EquationParameter("t63", "Time to 63.2% response (= τ + θ)", 
                              "time", "min", ParameterType.OUTPUT, symbol="t₆₃"),
            EquationParameter("t95", "Time to 95% response (= 3τ + θ)", 
                              "time", "min", ParameterType.OUTPUT, symbol="t₉₅"),
            EquationParameter("t99", "Time to 99% response (= 5τ + θ)", 
                              "time", "min", ParameterType.OUTPUT, symbol="t₉₉"),
            EquationParameter("final_value", "Final steady-state change (= K × M)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Δy∞"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau: pint.Quantity, step_size: pint.Quantity,
                   theta: pint.Quantity = None, **kwargs) -> Dict[str, pint.Quantity]:
        k = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau_val = tau.to('min').magnitude
        m = step_size.magnitude if hasattr(step_size, 'magnitude') else float(step_size)
        theta_val = theta.to('min').magnitude if theta else 0
        
        t63 = tau_val + theta_val
        t95 = 3 * tau_val + theta_val
        t99 = 5 * tau_val + theta_val
        final_val = k * m
        
        return {
            "t63": ureg.Quantity(t63, "min"),
            "t95": ureg.Quantity(t95, "min"),
            "t99": ureg.Quantity(t99, "min"),
            "final_value": ureg.Quantity(final_val, "")
        }


class SecondOrderResponse(BaseEquation):
    """Second-order system response analysis."""
    
    equation_id = "second_order_response"
    name = "Second-Order Step Response"
    category = "Process Control"
    description = "Analyze overshoot, settling time for second-order systems"
    reference = "Seborg et al., Process Dynamics and Control"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain", "dimensionless", "",
                              ParameterType.INPUT, symbol="K"),
            EquationParameter("tau", "Time constant", "time", "min",
                              ParameterType.INPUT, symbol="τ"),
            EquationParameter("zeta", "Damping ratio - determines oscillation behavior", 
                              "dimensionless", "", ParameterType.INPUT, symbol="ζ", 
                              typical_range=(0, 2),
                              tooltip="<1: underdamped (oscillates), =1: critical, >1: overdamped"),
            EquationParameter("response_type", "Response type: 0=under, 1=critical, 2=over", 
                              "dimensionless", "", ParameterType.OUTPUT),
            EquationParameter("overshoot", "Percent overshoot - first peak above final value", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="OS%",
                              tooltip="Only for underdamped (ζ<1)"),
            EquationParameter("decay_ratio", "Decay ratio - ratio of successive peaks", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="DR"),
            EquationParameter("rise_time", "Rise time - time from 10% to 90%", 
                              "time", "min", ParameterType.OUTPUT, symbol="tᵣ"),
            EquationParameter("settling_time", "2% settling time", 
                              "time", "min", ParameterType.OUTPUT, symbol="tₛ"),
            EquationParameter("period", "Oscillation period (underdamped only)", 
                              "time", "min", ParameterType.OUTPUT, symbol="T"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau: pint.Quantity,
                   zeta: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        k = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau_val = tau.to('min').magnitude
        z = zeta.magnitude if hasattr(zeta, 'magnitude') else float(zeta)
        
        if z < 1:
            resp_type = 0  # Underdamped
            wd = np.sqrt(1 - z**2) / tau_val
            overshoot = 100 * np.exp(-np.pi * z / np.sqrt(1 - z**2))
            decay_ratio = np.exp(-2 * np.pi * z / np.sqrt(1 - z**2))
            period = 2 * np.pi / wd
            rise_time = (np.pi - np.arccos(z)) / wd
        elif z == 1:
            resp_type = 1  # Critically damped
            overshoot = 0
            decay_ratio = 0
            period = float('inf')
            rise_time = 2.2 * tau_val
        else:
            resp_type = 2  # Overdamped
            overshoot = 0
            decay_ratio = 0
            period = float('inf')
            rise_time = 2.2 * tau_val * z
        
        settling_time = 4 * tau_val / z if z > 0 else float('inf')
        
        return {
            "response_type": ureg.Quantity(resp_type, ""),
            "overshoot": ureg.Quantity(overshoot, ""),
            "decay_ratio": ureg.Quantity(decay_ratio, ""),
            "rise_time": ureg.Quantity(rise_time, "min"),
            "settling_time": ureg.Quantity(settling_time, "min"),
            "period": ureg.Quantity(period, "min")
        }


class DecayRatio(BaseEquation):
    """Calculate damping ratio from decay ratio."""
    
    equation_id = "decay_ratio"
    name = "Decay Ratio Analysis"
    category = "Process Control"
    description = "Convert measured decay ratio to damping ratio"
    reference = "Standard control theory"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("DR", "Decay ratio - measured as (2nd peak)/(1st peak)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="DR", 
                              typical_range=(0, 1),
                              tooltip="Quarter-decay criterion: DR = 0.25 is common target"),
            EquationParameter("zeta", "Damping ratio", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="ζ"),
            EquationParameter("overshoot", "Percent overshoot for this decay ratio", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="OS%"),
        ]
    
    def _calculate(self, DR: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        dr = DR.magnitude if hasattr(DR, 'magnitude') else float(DR)
        
        if dr <= 0 or dr >= 1:
            raise ValueError("Decay ratio must be between 0 and 1")
        
        ln_dr = np.log(dr)
        zeta = -ln_dr / np.sqrt(4 * np.pi**2 + ln_dr**2)
        overshoot = np.sqrt(dr) * 100
        
        return {
            "zeta": ureg.Quantity(zeta, ""),
            "overshoot": ureg.Quantity(overshoot, "")
        }


class TransmitterScaling(BaseEquation):
    """Calculate transmitter output from process variable."""
    
    equation_id = "transmitter_scaling"
    name = "Transmitter Scaling (4-20mA)"
    category = "Process Control"
    description = "Convert between process variable and 4-20 mA signal"
    reference = "Standard instrumentation"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("PV", "Process variable - actual measured value", 
                              "dimensionless", "", ParameterType.INPUT, symbol="PV"),
            EquationParameter("LRV", "Lower range value - PV at 4 mA", 
                              "dimensionless", "", ParameterType.INPUT, symbol="LRV",
                              tooltip="0% of span"),
            EquationParameter("URV", "Upper range value - PV at 20 mA", 
                              "dimensionless", "", ParameterType.INPUT, symbol="URV",
                              tooltip="100% of span"),
            EquationParameter("mA", "Current output signal", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="mA",
                              tooltip="4-20 mA range"),
            EquationParameter("percent", "Percent of span", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="%"),
        ]
    
    def _calculate(self, PV: pint.Quantity, LRV: pint.Quantity, 
                   URV: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        pv = PV.magnitude if hasattr(PV, 'magnitude') else float(PV)
        lrv = LRV.magnitude if hasattr(LRV, 'magnitude') else float(LRV)
        urv = URV.magnitude if hasattr(URV, 'magnitude') else float(URV)
        
        span = urv - lrv
        percent = (pv - lrv) / span * 100
        ma = 4 + (pv - lrv) / span * 16
        
        return {
            "mA": ureg.Quantity(ma, ""),
            "percent": ureg.Quantity(percent, "")
        }


class LambdaTuning(BaseEquation):
    """Lambda tuning method for PID controllers - most popular in process industry."""
    
    equation_id = "lambda_tuning"
    name = "Lambda (λ) Tuning"
    category = "Process Control"
    description = "Industry-standard tuning with adjustable closed-loop time constant λ"
    reference = "Dahlin, 1968 / Lambda Tuning Method"
    learning_content = LAMBDA_TUNING_LEARNING
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain - steady-state ΔOutput/ΔInput", 
                              "dimensionless", "", ParameterType.INPUT, symbol="K",
                              tooltip="From step test: final change / step size"),
            EquationParameter("tau", "Process time constant - time to 63.2% response", 
                              "time", "min", ParameterType.INPUT, symbol="τ",
                              tooltip="From step test or process knowledge"),
            EquationParameter("theta", "Dead time (delay) - time before response starts", 
                              "time", "min", ParameterType.INPUT, symbol="θ",
                              tooltip="Also called transportation lag"),
            EquationParameter("lambda_", "Desired closed-loop time constant", 
                              "time", "min", ParameterType.INPUT, symbol="λ",
                              tooltip="Rule of thumb: λ ≥ 3θ for robustness, λ = τ for moderate speed"),
            EquationParameter("Kc", "Controller gain (Kp)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Kc"),
            EquationParameter("Ti", "Integral time", 
                              "time", "min", ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time (set to 0 for PI-only)", 
                              "time", "min", ParameterType.OUTPUT, symbol="Tᵈ"),
            EquationParameter("recommended_lambda", "Recommended λ (= max(τ, 3θ))", 
                              "time", "min", ParameterType.OUTPUT, symbol="λ_rec",
                              tooltip="Conservative starting point"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau: pint.Quantity, theta: pint.Quantity,
                   lambda_: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        k = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau_val = tau.to('min').magnitude
        theta_val = theta.to('min').magnitude
        lam = lambda_.to('min').magnitude
        
        # Lambda tuning for FOPDT: Kc = τ / (K(λ + θ))
        kc = tau_val / (k * (lam + theta_val))
        ti = tau_val  # Integral time = process time constant
        td = 0  # Usually PI-only, but can use θ/2 for PID
        
        # Recommended lambda for robustness
        rec_lambda = max(tau_val, 3 * theta_val)
        
        return {
            "Kc": ureg.Quantity(kc, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min"),
            "recommended_lambda": ureg.Quantity(rec_lambda, "min")
        }


class SIMCTuning(BaseEquation):
    """SIMC (Skogestad IMC) tuning rules - simple and effective."""
    
    equation_id = "simc_tuning"
    name = "SIMC Tuning (Skogestad)"
    category = "Process Control"
    description = "Simple IMC-based rules: fast yet robust PID tuning"
    reference = "Skogestad, 2003"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("K", "Process gain", "dimensionless", "", 
                              ParameterType.INPUT, symbol="K"),
            EquationParameter("tau1", "Dominant time constant", "time", "min", 
                              ParameterType.INPUT, symbol="τ₁"),
            EquationParameter("theta", "Dead time", "time", "min", 
                              ParameterType.INPUT, symbol="θ"),
            EquationParameter("tau_c", "Desired closed-loop time constant", 
                              "time", "min", ParameterType.INPUT, symbol="τc",
                              tooltip="SIMC recommends τc = θ for fast response"),
            EquationParameter("Kc", "Controller gain", "dimensionless", "", 
                              ParameterType.OUTPUT, symbol="Kc"),
            EquationParameter("Ti", "Integral time", "time", "min", 
                              ParameterType.OUTPUT, symbol="Tᵢ"),
            EquationParameter("Td", "Derivative time", "time", "min", 
                              ParameterType.OUTPUT, symbol="Tᵈ"),
        ]
    
    def _calculate(self, K: pint.Quantity, tau1: pint.Quantity, theta: pint.Quantity,
                   tau_c: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        k = K.magnitude if hasattr(K, 'magnitude') else float(K)
        tau1_val = tau1.to('min').magnitude
        theta_val = theta.to('min').magnitude
        tauc = tau_c.to('min').magnitude
        
        # SIMC PI rules
        kc = tau1_val / (k * (tauc + theta_val))
        ti = min(tau1_val, 4 * (tauc + theta_val))
        td = 0  # PI controller
        
        return {
            "Kc": ureg.Quantity(kc, ""),
            "Ti": ureg.Quantity(ti, "min"),
            "Td": ureg.Quantity(td, "min")
        }


class ControlLoopPerformance(BaseEquation):
    """Calculate control loop performance metrics."""
    
    equation_id = "loop_performance"
    name = "Control Loop Performance"
    category = "Process Control"
    description = "Calculate IAE, ISE, and other performance metrics from response data"
    reference = "Standard Control Theory"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("overshoot", "Percent overshoot from step test", 
                              "dimensionless", "", ParameterType.INPUT, symbol="OS%",
                              typical_range=(0, 100)),
            EquationParameter("settling_time", "Time to settle within 2% of final value", 
                              "time", "min", ParameterType.INPUT, symbol="tₛ"),
            EquationParameter("rise_time", "Time from 10% to 90% of final value", 
                              "time", "min", ParameterType.INPUT, symbol="tᵣ"),
            EquationParameter("period", "Oscillation period (if any)", 
                              "time", "min", ParameterType.INPUT, symbol="T", required=False),
            EquationParameter("damping_ratio", "Estimated damping ratio", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="ζ"),
            EquationParameter("natural_freq", "Estimated natural frequency", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="ωn"),
            EquationParameter("quality", "Loop quality assessment (1-10)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Quality"),
        ]
    
    def _calculate(self, overshoot: pint.Quantity, settling_time: pint.Quantity,
                   rise_time: pint.Quantity, period: pint.Quantity = None, 
                   **kwargs) -> Dict[str, pint.Quantity]:
        os = overshoot.magnitude if hasattr(overshoot, 'magnitude') else float(overshoot)
        ts = settling_time.to('min').magnitude
        tr = rise_time.to('min').magnitude
        
        # Estimate damping from overshoot
        if os > 0:
            ln_os = np.log(os / 100)
            zeta = -ln_os / np.sqrt(np.pi**2 + ln_os**2)
        else:
            zeta = 1.0  # Critically damped or overdamped
        
        # Natural frequency from settling time
        if zeta > 0 and ts > 0:
            omega_n = 4 / (zeta * ts)
        else:
            omega_n = 0
        
        # Quality score (1-10, 10 is best)
        # Penalize high overshoot and slow response
        quality = 10 - min(os/10, 5) - min(ts/tr/10, 5) if tr > 0 else 5
        quality = max(1, min(10, quality))
        
        return {
            "damping_ratio": ureg.Quantity(zeta, ""),
            "natural_freq": ureg.Quantity(omega_n, ""),
            "quality": ureg.Quantity(quality, "")
        }


class ControlValveTravel(BaseEquation):
    """Calculate valve travel/position from flow."""
    
    equation_id = "cv_travel"
    name = "Control Valve Travel"
    category = "Process Control"
    description = "Calculate valve percent open for equal percentage characteristic"
    reference = "ISA-75.01.01"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("Q", "Current flow rate", "volumetric_flow", "gpm", 
                              ParameterType.INPUT, symbol="Q"),
            EquationParameter("Q_max", "Maximum flow at full open", "volumetric_flow", "gpm", 
                              ParameterType.INPUT, symbol="Qmax"),
            EquationParameter("R", "Rangeability - ratio of max to min controllable flow", 
                              "dimensionless", "", ParameterType.INPUT, symbol="R",
                              typical_range=(25, 100), tooltip="Typically 50:1 for globe valves"),
            EquationParameter("characteristic", "Valve characteristic (1=linear, 2=eq%, 3=quick)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="Char",
                              tooltip="1=Linear, 2=Equal %, 3=Quick opening"),
            EquationParameter("travel", "Valve travel (percent open)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="%Open"),
            EquationParameter("Cv_rel", "Relative Cv (Cv/Cv_max)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Cv/Cvmax"),
        ]
    
    def _calculate(self, Q: pint.Quantity, Q_max: pint.Quantity, R: pint.Quantity,
                   characteristic: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        q = Q.to('gpm').magnitude
        qmax = Q_max.to('gpm').magnitude
        r = R.magnitude if hasattr(R, 'magnitude') else float(R)
        char = characteristic.magnitude if hasattr(characteristic, 'magnitude') else 2
        
        cv_rel = q / qmax
        
        if char >= 3:  # Quick opening
            travel = np.sqrt(cv_rel) * 100
        elif char >= 2:  # Equal percentage
            if cv_rel > 0:
                travel = (1 + np.log(cv_rel) / np.log(r)) * 100
            else:
                travel = 0
        else:  # Linear
            travel = cv_rel * 100
        
        travel = max(0, min(100, travel))
        
        return {
            "travel": ureg.Quantity(travel, ""),
            "Cv_rel": ureg.Quantity(cv_rel, "")
        }


class ProcessCapability(BaseEquation):
    """Calculate process capability indices Cp and Cpk."""
    
    equation_id = "process_capability"
    name = "Process Capability (Cp, Cpk)"
    category = "Process Control"
    description = "Statistical process capability indices for quality control"
    reference = "Six Sigma / SPC Standards"
    
    def get_parameters(self) -> List[EquationParameter]:
        return [
            EquationParameter("USL", "Upper specification limit", 
                              "dimensionless", "", ParameterType.INPUT, symbol="USL"),
            EquationParameter("LSL", "Lower specification limit", 
                              "dimensionless", "", ParameterType.INPUT, symbol="LSL"),
            EquationParameter("mean", "Process mean (average)", 
                              "dimensionless", "", ParameterType.INPUT, symbol="μ"),
            EquationParameter("std_dev", "Process standard deviation", 
                              "dimensionless", "", ParameterType.INPUT, symbol="σ"),
            EquationParameter("Cp", "Process capability - potential capability", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Cp",
                              tooltip="≥1.33 typically required"),
            EquationParameter("Cpk", "Process capability - actual capability (accounts for centering)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="Cpk",
                              tooltip="≥1.33 typically required"),
            EquationParameter("sigma_level", "Sigma level (Six Sigma metric)", 
                              "dimensionless", "", ParameterType.OUTPUT, symbol="σ-level"),
        ]
    
    def _calculate(self, USL: pint.Quantity, LSL: pint.Quantity, mean: pint.Quantity,
                   std_dev: pint.Quantity, **kwargs) -> Dict[str, pint.Quantity]:
        usl = USL.magnitude if hasattr(USL, 'magnitude') else float(USL)
        lsl = LSL.magnitude if hasattr(LSL, 'magnitude') else float(LSL)
        mu = mean.magnitude if hasattr(mean, 'magnitude') else float(mean)
        sigma = std_dev.magnitude if hasattr(std_dev, 'magnitude') else float(std_dev)
        
        cp = (usl - lsl) / (6 * sigma)
        cpu = (usl - mu) / (3 * sigma)
        cpl = (mu - lsl) / (3 * sigma)
        cpk = min(cpu, cpl)
        sigma_level = 3 * cpk
        
        return {
            "Cp": ureg.Quantity(cp, ""),
            "Cpk": ureg.Quantity(cpk, ""),
            "sigma_level": ureg.Quantity(sigma_level, "")
        }


# Registry of all process control equations
PROCESS_CONTROL_EQUATIONS = {
    'ziegler_nichols_pid': ZieglerNicholsPID,
    'cohen_coon_pid': CohenCoonPID,
    'imc_pid': IMCPID,
    'tyreus_luyben_pid': TyreusLuybenPID,
    'lambda_tuning': LambdaTuning,
    'simc_tuning': SIMCTuning,
    'cv_liquid': ControlValveCvLiquid,
    'cv_gas': ControlValveCvGas,
    'cv_travel': ControlValveTravel,
    'gain_margin': ControlLoopGainMargin,
    'loop_performance': ControlLoopPerformance,
    'process_capability': ProcessCapability,
    'first_order_response': FirstOrderResponse,
    'second_order_response': SecondOrderResponse,
    'decay_ratio': DecayRatio,
    'transmitter_scaling': TransmitterScaling,
}

