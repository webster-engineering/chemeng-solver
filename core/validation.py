"""
Validation system for dimensional analysis and benchmark comparison.
"""

import json
from pathlib import Path
from enum import Enum
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Tuple
import pint

from .units import ureg, UnitConverter, get_converter


class ValidationStatus(Enum):
    """Validation result status."""
    PASS = "pass"
    WARNING = "warning"
    FAIL = "fail"
    SKIPPED = "skipped"


@dataclass
class ValidationMessage:
    """A single validation message."""
    status: ValidationStatus
    message: str
    details: str = ""


@dataclass
class ValidationResult:
    """Complete validation result for a calculation."""
    overall_status: ValidationStatus
    messages: List[ValidationMessage] = field(default_factory=list)
    dimension_check: Optional[ValidationMessage] = None
    range_check: Optional[ValidationMessage] = None
    benchmark_check: Optional[ValidationMessage] = None
    
    def add_message(self, status: ValidationStatus, message: str, details: str = "") -> None:
        """Add a validation message."""
        msg = ValidationMessage(status, message, details)
        self.messages.append(msg)
        
        # Update overall status (worst case)
        if status == ValidationStatus.FAIL:
            self.overall_status = ValidationStatus.FAIL
        elif status == ValidationStatus.WARNING and self.overall_status != ValidationStatus.FAIL:
            self.overall_status = ValidationStatus.WARNING
    
    def is_valid(self) -> bool:
        """Check if validation passed (no failures)."""
        return self.overall_status != ValidationStatus.FAIL


class Validator:
    """
    Multi-layer validation system for chemical engineering calculations.
    """
    
    def __init__(self, benchmarks_file: Optional[Path] = None):
        """
        Initialize validator.
        
        Args:
            benchmarks_file: Path to JSON file with benchmark values
        """
        self.converter = get_converter()
        self.benchmarks: Dict[str, Any] = {}
        
        if benchmarks_file is None:
            benchmarks_file = Path(__file__).parent.parent / "benchmarks" / "reference_values.json"
        
        if benchmarks_file.exists():
            self.load_benchmarks(benchmarks_file)
    
    def load_benchmarks(self, filepath: Path) -> None:
        """Load benchmark values from JSON file."""
        try:
            with open(filepath, 'r') as f:
                self.benchmarks = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load benchmarks: {e}")
    
    def validate_calculation(
        self,
        equation_id: str,
        inputs: Dict[str, pint.Quantity],
        outputs: Dict[str, pint.Quantity],
        expected_dimensions: Optional[Dict[str, str]] = None,
        valid_ranges: Optional[Dict[str, Tuple[float, float, str]]] = None
    ) -> ValidationResult:
        """
        Perform complete validation on a calculation.
        
        Args:
            equation_id: Identifier for the equation (for benchmark lookup)
            inputs: Dictionary of input quantities
            outputs: Dictionary of output quantities
            expected_dimensions: Expected dimension for each output {name: dimension_string}
            valid_ranges: Valid ranges for outputs {name: (min, max, unit)}
        
        Returns:
            ValidationResult with all checks
        """
        result = ValidationResult(overall_status=ValidationStatus.PASS)
        
        # Dimensional analysis
        if expected_dimensions:
            result.dimension_check = self._check_dimensions(outputs, expected_dimensions)
            if result.dimension_check.status != ValidationStatus.PASS:
                result.add_message(
                    result.dimension_check.status,
                    result.dimension_check.message,
                    result.dimension_check.details
                )
        
        # Range checks
        if valid_ranges:
            result.range_check = self._check_ranges(outputs, valid_ranges)
            if result.range_check.status != ValidationStatus.PASS:
                result.add_message(
                    result.range_check.status,
                    result.range_check.message,
                    result.range_check.details
                )
        
        # Benchmark comparison
        if equation_id in self.benchmarks:
            result.benchmark_check = self._check_benchmark(equation_id, inputs, outputs)
            if result.benchmark_check.status != ValidationStatus.PASS:
                result.add_message(
                    result.benchmark_check.status,
                    result.benchmark_check.message,
                    result.benchmark_check.details
                )
        
        # Physical sanity checks
        sanity_messages = self._sanity_checks(outputs)
        for msg in sanity_messages:
            result.add_message(msg.status, msg.message, msg.details)
        
        return result
    
    def _check_dimensions(
        self,
        outputs: Dict[str, pint.Quantity],
        expected: Dict[str, str]
    ) -> ValidationMessage:
        """Check dimensional consistency of outputs."""
        errors = []
        
        for name, quantity in outputs.items():
            if name in expected:
                expected_dim = expected[name]
                if not self.converter.check_dimensionality(quantity, expected_dim):
                    actual_dim = self.converter.get_dimension_string(quantity)
                    errors.append(f"{name}: expected {expected_dim}, got {actual_dim}")
        
        if errors:
            return ValidationMessage(
                status=ValidationStatus.FAIL,
                message="Dimensional analysis failed",
                details="; ".join(errors)
            )
        
        return ValidationMessage(
            status=ValidationStatus.PASS,
            message="Dimensional analysis passed"
        )
    
    def _check_ranges(
        self,
        outputs: Dict[str, pint.Quantity],
        valid_ranges: Dict[str, Tuple[float, float, str]]
    ) -> ValidationMessage:
        """Check if outputs fall within valid ranges."""
        warnings = []
        errors = []
        
        for name, quantity in outputs.items():
            if name in valid_ranges:
                min_val, max_val, unit = valid_ranges[name]
                try:
                    value = quantity.to(unit).magnitude
                    
                    if value < min_val or value > max_val:
                        msg = f"{name}: {value:.4g} {unit} outside range [{min_val}, {max_val}]"
                        # Distinguish between slightly out of range (warning) and way out (error)
                        range_size = max_val - min_val if max_val > min_val else 1
                        if value < min_val - 0.5 * range_size or value > max_val + 0.5 * range_size:
                            errors.append(msg)
                        else:
                            warnings.append(msg)
                except pint.DimensionalityError:
                    errors.append(f"{name}: cannot convert to {unit} for range check")
        
        if errors:
            return ValidationMessage(
                status=ValidationStatus.FAIL,
                message="Output values outside valid range",
                details="; ".join(errors + warnings)
            )
        elif warnings:
            return ValidationMessage(
                status=ValidationStatus.WARNING,
                message="Output values near range limits",
                details="; ".join(warnings)
            )
        
        return ValidationMessage(
            status=ValidationStatus.PASS,
            message="All outputs within valid ranges"
        )
    
    def _check_benchmark(
        self,
        equation_id: str,
        inputs: Dict[str, pint.Quantity],
        outputs: Dict[str, pint.Quantity]
    ) -> ValidationMessage:
        """Compare results against benchmark values."""
        benchmark = self.benchmarks.get(equation_id, {})
        cases = benchmark.get('cases', [])
        
        if not cases:
            return ValidationMessage(
                status=ValidationStatus.SKIPPED,
                message="No benchmark cases available"
            )
        
        # Find matching benchmark case
        for case in cases:
            case_inputs = case.get('inputs', {})
            if self._inputs_match(inputs, case_inputs):
                return self._compare_outputs(outputs, case.get('outputs', {}), case.get('tolerance', 0.05))
        
        return ValidationMessage(
            status=ValidationStatus.SKIPPED,
            message="No matching benchmark case found for given inputs"
        )
    
    def _inputs_match(
        self,
        actual: Dict[str, pint.Quantity],
        expected: Dict[str, dict]
    ) -> bool:
        """Check if actual inputs match expected benchmark inputs."""
        for name, expected_val in expected.items():
            if name not in actual:
                return False
            
            try:
                expected_q = ureg.Quantity(expected_val['value'], expected_val['unit'])
                actual_q = actual[name]
                
                # Convert to same units and compare
                if not actual_q.dimensionality == expected_q.dimensionality:
                    return False
                
                actual_converted = actual_q.to(expected_q.units).magnitude
                tolerance = expected_val.get('tolerance', 0.01)
                
                if abs(actual_converted - expected_q.magnitude) > tolerance * abs(expected_q.magnitude):
                    return False
            except Exception:
                return False
        
        return True
    
    def _compare_outputs(
        self,
        actual: Dict[str, pint.Quantity],
        expected: Dict[str, dict],
        tolerance: float
    ) -> ValidationMessage:
        """Compare actual outputs against expected benchmark values."""
        comparisons = []
        all_pass = True
        
        for name, expected_val in expected.items():
            if name not in actual:
                comparisons.append(f"{name}: missing")
                all_pass = False
                continue
            
            try:
                expected_q = ureg.Quantity(expected_val['value'], expected_val['unit'])
                actual_q = actual[name].to(expected_q.units)
                
                error = abs(actual_q.magnitude - expected_q.magnitude) / abs(expected_q.magnitude) if expected_q.magnitude != 0 else abs(actual_q.magnitude)
                
                if error <= tolerance:
                    comparisons.append(f"{name}: {actual_q.magnitude:.4g} ≈ {expected_q.magnitude:.4g} {expected_q.units} (within {tolerance*100:.1f}%)")
                else:
                    comparisons.append(f"{name}: {actual_q.magnitude:.4g} ≠ {expected_q.magnitude:.4g} {expected_q.units} (error: {error*100:.1f}%)")
                    all_pass = False
            except Exception as e:
                comparisons.append(f"{name}: comparison error - {e}")
                all_pass = False
        
        if all_pass:
            return ValidationMessage(
                status=ValidationStatus.PASS,
                message="Results match benchmark values",
                details="; ".join(comparisons)
            )
        else:
            return ValidationMessage(
                status=ValidationStatus.WARNING,
                message="Results differ from benchmark values",
                details="; ".join(comparisons)
            )
    
    def _sanity_checks(self, outputs: Dict[str, pint.Quantity]) -> List[ValidationMessage]:
        """Perform physical sanity checks on outputs."""
        messages = []
        
        for name, quantity in outputs.items():
            # Check for negative values where they don't make sense
            name_lower = name.lower()
            
            if quantity.magnitude < 0:
                # These should never be negative
                if any(word in name_lower for word in ['pressure', 'flow', 'cv', 'area', 'volume', 'diameter', 'length', 'mass']):
                    messages.append(ValidationMessage(
                        status=ValidationStatus.FAIL,
                        message=f"Negative value for {name}",
                        details=f"{name} = {quantity.magnitude:.4g} {quantity.units} (should be positive)"
                    ))
            
            # Check for unreasonably large values
            if quantity.magnitude > 1e15:
                messages.append(ValidationMessage(
                    status=ValidationStatus.WARNING,
                    message=f"Unusually large value for {name}",
                    details=f"{name} = {quantity.magnitude:.4g} {quantity.units}"
                ))
            
            # Check for NaN or Inf
            import math
            if math.isnan(quantity.magnitude) or math.isinf(quantity.magnitude):
                messages.append(ValidationMessage(
                    status=ValidationStatus.FAIL,
                    message=f"Invalid value for {name}",
                    details=f"{name} is NaN or Infinity"
                ))
        
        return messages
    
    def validate_input(
        self,
        name: str,
        value: pint.Quantity,
        expected_dimension: str,
        valid_range: Optional[Tuple[float, float, str]] = None
    ) -> ValidationResult:
        """
        Validate a single input value.
        
        Args:
            name: Parameter name
            value: Input quantity
            expected_dimension: Expected dimension string
            valid_range: Optional (min, max, unit) tuple
        
        Returns:
            ValidationResult
        """
        result = ValidationResult(overall_status=ValidationStatus.PASS)
        
        # Check dimension
        if not self.converter.check_dimensionality(value, expected_dimension):
            actual_dim = self.converter.get_dimension_string(value)
            result.add_message(
                ValidationStatus.FAIL,
                f"Invalid dimension for {name}",
                f"Expected {expected_dimension}, got {actual_dim}"
            )
        
        # Check range
        if valid_range:
            min_val, max_val, unit = valid_range
            try:
                val = value.to(unit).magnitude
                if val < min_val or val > max_val:
                    result.add_message(
                        ValidationStatus.WARNING,
                        f"{name} outside typical range",
                        f"Value {val:.4g} {unit} outside [{min_val}, {max_val}]"
                    )
            except pint.DimensionalityError:
                result.add_message(
                    ValidationStatus.FAIL,
                    f"Cannot validate range for {name}",
                    f"Unit conversion failed"
                )
        
        return result
