"""
Equation Verification Script
Tests all equations against reference/known values.
No external testing framework required.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import json
from typing import Dict, Any, List, Tuple
import traceback

# Import all equation modules
from equations.process_control import PROCESS_CONTROL_EQUATIONS
from equations.fluid_dynamics import FLUID_DYNAMICS_EQUATIONS
from equations.heat_transfer import HEAT_TRANSFER_EQUATIONS
from equations.thermodynamics import THERMODYNAMICS_EQUATIONS
from equations.mass_transfer import MASS_TRANSFER_EQUATIONS
from equations.distillation import DISTILLATION_EQUATIONS
from equations.reaction_kinetics import REACTION_KINETICS_EQUATIONS
from equations.piping import PIPING_EQUATIONS
from equations.pumps import PUMPS_EQUATIONS
from equations.vessels import VESSELS_EQUATIONS
from equations.safety import SAFETY_EQUATIONS
from equations.economics import ECONOMICS_EQUATIONS
from equations.instrumentation import INSTRUMENTATION_EQUATIONS
from equations.basic_math import BASIC_MATH_EQUATIONS

# Comprehensive test cases with reference values
TEST_CASES = {
    # ==================== Process Control ====================
    "ziegler_nichols_pid": {
        "module": "process_control",
        "description": "Ziegler-Nichols PID tuning (Ziegler & Nichols, 1942)",
        "cases": [
            {
                "inputs": {"Ku": 2.0, "Pu": (5.0, "min")},
                "expected": {"Kp": 1.2, "Ti": (2.5, "min"), "Td": (0.625, "min")},
                "tolerance": 0.01
            },
            {
                "inputs": {"Ku": 4.0, "Pu": (10.0, "s")},
                "expected": {"Kp": 2.4, "Ti": (5.0, "s"), "Td": (1.25, "s")},
                "tolerance": 0.01
            }
        ]
    },
    "cv_liquid": {
        "module": "process_control", 
        "description": "Control Valve Cv for liquid (ISA-75.01.01)",
        "cases": [
            {
                "inputs": {"Q": (100, "gpm"), "dP": (25, "psi"), "SG": 1.0},
                "expected": {"Cv": 20.0},
                "tolerance": 0.01
            },
            {
                "inputs": {"Q": (500, "gpm"), "dP": (100, "psi"), "SG": 0.85},
                "expected": {"Cv": 46.1},
                "tolerance": 0.02
            }
        ]
    },
    
    # ==================== Fluid Dynamics ====================
    "reynolds_number": {
        "module": "fluid_dynamics",
        "description": "Reynolds Number calculation",
        "cases": [
            {
                "inputs": {
                    "rho": (1000, "kg/m**3"), 
                    "V": (1, "m/s"), 
                    "D": (0.1, "m"), 
                    "mu": (0.001, "Pa*s")
                },
                "expected": {"Re": 100000},
                "tolerance": 0.001
            },
            {
                "inputs": {
                    "rho": (62.4, "lb/ft**3"),
                    "V": (10, "ft/s"),
                    "D": (4, "inch"),
                    "mu": (0.001, "Pa*s")
                },
                "expected": {"Re": 101400},  # Approximate from unit conversion
                "tolerance": 0.05
            }
        ]
    },
    "darcy_weisbach": {
        "module": "fluid_dynamics",
        "description": "Darcy-Weisbach pressure drop (Perry's 8th Ed)",
        "cases": [
            {
                "inputs": {
                    "f": 0.02,
                    "L": (100, "ft"),
                    "D": (4, "inch"),
                    "V": (10, "ft/s"),
                    "rho": (62.4, "lb/ft**3")
                },
                "expected": {"dP": (6.48, "psi")},
                "tolerance": 0.05
            }
        ]
    },
    
    # ==================== Heat Transfer ====================
    "lmtd": {
        "module": "heat_transfer",
        "description": "Log Mean Temperature Difference (Perry's)",
        "cases": [
            {
                "inputs": {
                    "T_hot_in": (300, "degF"),
                    "T_hot_out": (200, "degF"),
                    "T_cold_in": (100, "degF"),
                    "T_cold_out": (150, "degF")
                },
                # LMTD = (150-50)/ln(150/50) = 100/1.099 = 91.0 for parallel
                # For counterflow: dt1=300-150=150, dt2=200-100=100
                # LMTD = (150-100)/ln(150/100) = 50/0.405 = 123.3
                "expected": {"LMTD": (123.3, "delta_degF")},
                "tolerance": 0.02
            }
        ]
    },
    
    # ==================== Thermodynamics ====================
    "ideal_gas": {
        "module": "thermodynamics",
        "description": "Ideal Gas Law PV=nRT",
        "cases": [
            {
                "inputs": {
                    "P": (101.325, "kPa"),
                    "n": (1, "mol"),
                    "T": (273.15, "K")
                },
                "expected": {"V": (22.41, "L")},  # Molar volume at STP
                "tolerance": 0.01
            }
        ]
    },
    
    # ==================== Safety ====================
    "relief_valve_api": {
        "module": "safety",
        "description": "API 520 Relief Valve Sizing",
        "cases": [
            {
                "inputs": {
                    "W": (10000, "lb/hr"),
                    "P1": (150, "psia"),
                    "k": 1.4,
                    "T": (520, "degR"),
                    "Z": 1.0,
                    "M": 29
                },
                # API 520 formula check
                "expected": {"A": (0.5, "inch**2")},  # Approximate
                "tolerance": 0.20  # Higher tolerance - verify manually
            }
        ]
    },
    
    # ==================== Basic Math ====================
    "quadratic_solver": {
        "module": "basic_math",
        "description": "Quadratic equation solver",
        "cases": [
            {
                "inputs": {"a": 1, "b": -5, "c": 6},
                "expected": {"x1": 3.0, "x2": 2.0},
                "tolerance": 0.0001
            },
            {
                "inputs": {"a": 1, "b": 0, "c": -4},
                "expected": {"x1": 2.0, "x2": -2.0},
                "tolerance": 0.0001
            }
        ]
    }
}


def get_equation_registry():
    """Get all equations organized by module."""
    return {
        "process_control": PROCESS_CONTROL_EQUATIONS,
        "fluid_dynamics": FLUID_DYNAMICS_EQUATIONS,
        "heat_transfer": HEAT_TRANSFER_EQUATIONS,
        "thermodynamics": THERMODYNAMICS_EQUATIONS,
        "mass_transfer": MASS_TRANSFER_EQUATIONS,
        "distillation": DISTILLATION_EQUATIONS,
        "reaction_kinetics": REACTION_KINETICS_EQUATIONS,
        "piping": PIPING_EQUATIONS,
        "pumps": PUMPS_EQUATIONS,
        "vessels": VESSELS_EQUATIONS,
        "safety": SAFETY_EQUATIONS,
        "economics": ECONOMICS_EQUATIONS,
        "instrumentation": INSTRUMENTATION_EQUATIONS,
        "basic_math": BASIC_MATH_EQUATIONS,
    }


def run_test_case(equation_obj, inputs: Dict, expected: Dict, tolerance: float) -> Tuple[bool, str, Dict]:
    """Run a single test case and return (passed, message, actual_outputs)."""
    try:
        result = equation_obj.calculate(inputs)
        
        if not result.success:
            return False, f"Calculation failed: {result.error}", {}
        
        actual_outputs = {}
        errors = []
        
        for output_name, expected_val in expected.items():
            # Handle expected value with unit
            if isinstance(expected_val, tuple):
                exp_value, exp_unit = expected_val
                actual = result.get_output_value(output_name, exp_unit)
            else:
                exp_value = expected_val
                actual = result.get_output_value(output_name)
            
            actual_outputs[output_name] = actual
            
            # Check tolerance
            if exp_value != 0:
                rel_error = abs(actual - exp_value) / abs(exp_value)
            else:
                rel_error = abs(actual - exp_value)
            
            if rel_error > tolerance:
                errors.append(f"{output_name}: expected {exp_value}, got {actual:.6g} (error: {rel_error*100:.2f}%)")
        
        if errors:
            return False, "; ".join(errors), actual_outputs
        else:
            return True, "PASS", actual_outputs
            
    except Exception as e:
        return False, f"Exception: {str(e)}", {}


def verify_all_equations():
    """Run all verification tests."""
    registry = get_equation_registry()
    
    results = {
        "total": 0,
        "passed": 0,
        "failed": 0,
        "errors": 0,
        "details": []
    }
    
    print("=" * 70)
    print("EQUATION SOLVER VERIFICATION REPORT")
    print("=" * 70)
    print()
    
    # Run test cases
    for eq_id, test_config in TEST_CASES.items():
        module_name = test_config["module"]
        description = test_config["description"]
        
        print(f"\n[{eq_id}] {description}")
        print("-" * 50)
        
        # Find equation in registry
        if module_name not in registry:
            print(f"  ERROR: Module '{module_name}' not found")
            results["errors"] += 1
            continue
            
        equation_obj = None
        for eq in registry[module_name]:
            if eq.equation_id == eq_id:
                equation_obj = eq
                break
        
        if not equation_obj:
            print(f"  ERROR: Equation '{eq_id}' not found in {module_name}")
            results["errors"] += 1
            continue
        
        # Run all test cases for this equation
        for i, case in enumerate(test_config["cases"], 1):
            results["total"] += 1
            inputs = case["inputs"]
            expected = case["expected"]
            tolerance = case.get("tolerance", 0.01)
            
            passed, message, actuals = run_test_case(equation_obj, inputs, expected, tolerance)
            
            if passed:
                results["passed"] += 1
                status = "✓ PASS"
            else:
                results["failed"] += 1
                status = "✗ FAIL"
            
            print(f"  Case {i}: {status}")
            if not passed:
                print(f"    {message}")
            
            results["details"].append({
                "equation": eq_id,
                "case": i,
                "passed": passed,
                "message": message,
                "inputs": str(inputs),
                "expected": str(expected),
                "actual": str(actuals)
            })
    
    # List all available equations
    print("\n" + "=" * 70)
    print("EQUATION INVENTORY")
    print("=" * 70)
    
    total_equations = 0
    tested_equations = set(TEST_CASES.keys())
    
    for module_name, equations in registry.items():
        print(f"\n{module_name.upper().replace('_', ' ')} ({len(equations)} equations)")
        for eq in equations:
            tested = "✓" if eq.equation_id in tested_equations else "○"
            print(f"  {tested} {eq.equation_id}: {eq.name}")
            total_equations += 1
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal equations available: {total_equations}")
    print(f"Equations with test cases: {len(TEST_CASES)}")
    print(f"Test coverage: {len(TEST_CASES)/total_equations*100:.1f}%")
    print()
    print(f"Test cases run: {results['total']}")
    print(f"  ✓ Passed: {results['passed']}")
    print(f"  ✗ Failed: {results['failed']}")
    print(f"  ! Errors: {results['errors']}")
    
    if results['failed'] == 0 and results['errors'] == 0:
        print("\n✓ ALL TESTS PASSED")
    else:
        print(f"\n✗ {results['failed'] + results['errors']} ISSUES FOUND")
    
    return results


if __name__ == "__main__":
    results = verify_all_equations()
    
    # Exit with appropriate code
    if results['failed'] > 0 or results['errors'] > 0:
        sys.exit(1)
    sys.exit(0)
