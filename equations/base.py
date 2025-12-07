"""
Base equation class that all equations inherit from.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Tuple
from enum import Enum
import pint

from core.units import ureg, get_converter
from core.validation import Validator, ValidationResult
from core.uncertainty import MonteCarloSimulator, UncertainParameter, UncertaintyResult


class ParameterType(Enum):
    INPUT = "input"
    OUTPUT = "output"


@dataclass
class EquationParameter:
    """Definition of an equation parameter."""
    name: str
    description: str
    unit_category: str
    default_unit: str
    param_type: ParameterType
    symbol: str = ""
    typical_range: Optional[Tuple[float, float]] = None
    required: bool = True
    tooltip: str = ""
    expected_dimension: str = ""
    
    def __post_init__(self):
        if not self.symbol:
            self.symbol = self.name
        if not self.tooltip:
            self.tooltip = self.description


@dataclass
class EquationResult:
    """Result from an equation calculation."""
    success: bool
    outputs: Dict[str, pint.Quantity] = field(default_factory=dict)
    validation: Optional[ValidationResult] = None
    uncertainty: Optional[Dict[str, UncertaintyResult]] = None
    error_message: str = ""
    calculation_steps: List[str] = field(default_factory=list)
    
    def get_output_value(self, name: str, unit: Optional[str] = None) -> Optional[float]:
        if name not in self.outputs:
            return None
        quantity = self.outputs[name]
        if unit:
            try:
                return quantity.to(unit).magnitude
            except pint.DimensionalityError:
                return None
        return quantity.magnitude
    
    def format_outputs(self) -> Dict[str, str]:
        converter = get_converter()
        return {name: converter.format_quantity(q) for name, q in self.outputs.items()}


class BaseEquation(ABC):
    """Abstract base class for all equations."""
    
    equation_id: str = ""
    name: str = ""
    category: str = ""
    description: str = ""
    reference: str = ""
    
    def __init__(self):
        self.converter = get_converter()
        self.validator = Validator()
        self._parameters: Optional[List[EquationParameter]] = None
    
    @property
    def parameters(self) -> List[EquationParameter]:
        if self._parameters is None:
            self._parameters = self.get_parameters()
        return self._parameters
    
    @property
    def inputs(self) -> List[EquationParameter]:
        return [p for p in self.parameters if p.param_type == ParameterType.INPUT]
    
    @property
    def outputs(self) -> List[EquationParameter]:
        return [p for p in self.parameters if p.param_type == ParameterType.OUTPUT]
    
    @abstractmethod
    def get_parameters(self) -> List[EquationParameter]:
        pass
    
    @abstractmethod
    def _calculate(self, **kwargs: pint.Quantity) -> Dict[str, pint.Quantity]:
        pass
    
    def calculate(self, inputs: Dict[str, Any], validate: bool = True,
                  with_uncertainty: bool = False,
                  uncertain_params: Optional[List[UncertainParameter]] = None) -> EquationResult:
        try:
            input_quantities = self._process_inputs(inputs)
            output_quantities = self._calculate(**input_quantities)
            
            validation_result = None
            if validate:
                validation_result = self._validate(input_quantities, output_quantities)
            
            uncertainty_results = None
            if with_uncertainty:
                uncertainty_results = self._run_uncertainty(input_quantities, uncertain_params)
            
            return EquationResult(success=True, outputs=output_quantities,
                                  validation=validation_result, uncertainty=uncertainty_results)
        except Exception as e:
            return EquationResult(success=False, error_message=str(e))
    
    def _process_inputs(self, inputs: Dict[str, Any]) -> Dict[str, pint.Quantity]:
        quantities = {}
        for param in self.inputs:
            if param.name not in inputs:
                if param.required:
                    raise ValueError(f"Missing required input: {param.name}")
                continue
            value = inputs[param.name]
            if isinstance(value, pint.Quantity):
                quantities[param.name] = value
            elif isinstance(value, tuple) and len(value) == 2:
                quantities[param.name] = ureg.Quantity(value[0], value[1])
            else:
                quantities[param.name] = ureg.Quantity(value, param.default_unit)
        return quantities
    
    def _validate(self, inputs: Dict[str, pint.Quantity],
                  outputs: Dict[str, pint.Quantity]) -> ValidationResult:
        expected_dims = {}
        valid_ranges = {}
        for param in self.outputs:
            if param.expected_dimension:
                expected_dims[param.name] = param.expected_dimension
            if param.typical_range:
                valid_ranges[param.name] = (*param.typical_range, param.default_unit)
        return self.validator.validate_calculation(self.equation_id, inputs, outputs,
                                                    expected_dims, valid_ranges)
    
    def _run_uncertainty(self, input_quantities: Dict[str, pint.Quantity],
                         custom_params: Optional[List[UncertainParameter]]) -> Dict[str, UncertaintyResult]:
        simulator = MonteCarloSimulator()
        if custom_params:
            uncertain_params = custom_params
            fixed_params = {n: q for n, q in input_quantities.items()
                           if n not in [p.name for p in custom_params]}
        else:
            uncertain_params = [UncertainParameter(name=n, nominal_value=q.magnitude,
                               unit=str(q.units), relative_std=0.05)
                               for n, q in input_quantities.items()]
            fixed_params = {}
        return simulator.run_simulation(self._calculate, uncertain_params, fixed_params, keep_samples=True)
    
    def get_info(self) -> Dict[str, Any]:
        return {
            'id': self.equation_id, 'name': self.name, 'category': self.category,
            'description': self.description, 'reference': self.reference,
            'inputs': [{'name': p.name, 'description': p.description, 'unit': p.default_unit}
                       for p in self.inputs],
            'outputs': [{'name': p.name, 'description': p.description, 'unit': p.default_unit}
                        for p in self.outputs]
        }
