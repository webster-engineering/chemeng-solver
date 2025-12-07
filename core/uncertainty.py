"""
Monte Carlo simulation engine for uncertainty analysis.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Callable, Any, Optional, Tuple
from enum import Enum
import pint

from .units import ureg, get_converter
from config import get_settings


class DistributionType(Enum):
    """Supported probability distributions for uncertainty."""
    NORMAL = "normal"
    UNIFORM = "uniform"
    TRIANGULAR = "triangular"
    LOGNORMAL = "lognormal"


@dataclass
class UncertainParameter:
    """Definition of an uncertain input parameter."""
    name: str
    nominal_value: float
    unit: str
    distribution: DistributionType = DistributionType.NORMAL
    
    # For normal distribution
    std_dev: Optional[float] = None  # absolute standard deviation
    relative_std: Optional[float] = None  # relative (as fraction of nominal)
    
    # For uniform distribution
    low: Optional[float] = None
    high: Optional[float] = None
    
    # For triangular distribution
    mode: Optional[float] = None
    
    def __post_init__(self):
        """Set default uncertainty if not specified."""
        if self.distribution == DistributionType.NORMAL:
            if self.std_dev is None and self.relative_std is None:
                self.relative_std = 0.05  # Default 5% relative uncertainty
    
    def generate_samples(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """
        Generate n random samples from the distribution.
        
        Args:
            n: Number of samples
            rng: NumPy random generator
        
        Returns:
            Array of samples
        """
        if self.distribution == DistributionType.NORMAL:
            std = self.std_dev if self.std_dev is not None else self.nominal_value * self.relative_std
            return rng.normal(self.nominal_value, std, n)
        
        elif self.distribution == DistributionType.UNIFORM:
            low = self.low if self.low is not None else self.nominal_value * 0.9
            high = self.high if self.high is not None else self.nominal_value * 1.1
            return rng.uniform(low, high, n)
        
        elif self.distribution == DistributionType.TRIANGULAR:
            low = self.low if self.low is not None else self.nominal_value * 0.9
            high = self.high if self.high is not None else self.nominal_value * 1.1
            mode = self.mode if self.mode is not None else self.nominal_value
            return rng.triangular(low, mode, high, n)
        
        elif self.distribution == DistributionType.LOGNORMAL:
            # Convert to log-normal parameters
            std = self.std_dev if self.std_dev is not None else self.nominal_value * (self.relative_std or 0.05)
            mean = self.nominal_value
            sigma = np.sqrt(np.log(1 + (std/mean)**2))
            mu = np.log(mean) - sigma**2/2
            return rng.lognormal(mu, sigma, n)
        
        else:
            raise ValueError(f"Unknown distribution: {self.distribution}")


@dataclass
class UncertaintyResult:
    """Results from Monte Carlo uncertainty analysis."""
    parameter_name: str
    unit: str
    
    # Basic statistics
    mean: float
    std_dev: float
    median: float
    
    # Confidence interval
    confidence_level: float
    ci_lower: float
    ci_upper: float
    
    # Distribution info
    min_value: float
    max_value: float
    percentile_5: float
    percentile_95: float
    
    # Raw samples (optional, for plotting)
    samples: Optional[np.ndarray] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary (excluding raw samples)."""
        return {
            'parameter_name': self.parameter_name,
            'unit': self.unit,
            'mean': self.mean,
            'std_dev': self.std_dev,
            'median': self.median,
            'confidence_level': self.confidence_level,
            'ci_lower': self.ci_lower,
            'ci_upper': self.ci_upper,
            'min_value': self.min_value,
            'max_value': self.max_value,
            'percentile_5': self.percentile_5,
            'percentile_95': self.percentile_95
        }
    
    def format_result(self, precision: int = 4) -> str:
        """Format result as string with uncertainty."""
        return f"{self.mean:.{precision}g} ± {self.std_dev:.{precision}g} {self.unit} (95% CI: [{self.ci_lower:.{precision}g}, {self.ci_upper:.{precision}g}])"


@dataclass
class SensitivityResult:
    """Sensitivity analysis results showing which inputs most affect outputs."""
    output_name: str
    input_sensitivities: Dict[str, float] = field(default_factory=dict)
    
    def get_ranked_inputs(self) -> List[Tuple[str, float]]:
        """Get inputs ranked by sensitivity (highest first)."""
        return sorted(self.input_sensitivities.items(), key=lambda x: abs(x[1]), reverse=True)


class MonteCarloSimulator:
    """
    Monte Carlo simulation engine for uncertainty propagation.
    """
    
    def __init__(self, n_samples: Optional[int] = None, seed: Optional[int] = None):
        """
        Initialize simulator.
        
        Args:
            n_samples: Number of Monte Carlo samples (default from settings)
            seed: Random seed for reproducibility
        """
        settings = get_settings()
        self.n_samples = n_samples or settings.monte_carlo_samples
        self.confidence_level = settings.confidence_level
        self.rng = np.random.default_rng(seed)
        self.converter = get_converter()
    
    def run_simulation(
        self,
        calculation_func: Callable[..., Dict[str, pint.Quantity]],
        uncertain_params: List[UncertainParameter],
        fixed_params: Optional[Dict[str, pint.Quantity]] = None,
        keep_samples: bool = False
    ) -> Dict[str, UncertaintyResult]:
        """
        Run Monte Carlo simulation for uncertainty propagation.
        
        Args:
            calculation_func: Function that takes parameter dict and returns output dict
            uncertain_params: List of uncertain input parameters
            fixed_params: Dictionary of fixed (non-uncertain) parameters
            keep_samples: Whether to keep raw samples in results
        
        Returns:
            Dictionary of UncertaintyResult for each output
        """
        fixed_params = fixed_params or {}
        
        # Generate samples for all uncertain parameters
        param_samples = {}
        for param in uncertain_params:
            param_samples[param.name] = param.generate_samples(self.n_samples, self.rng)
        
        # Run calculations for each sample
        output_samples: Dict[str, List[float]] = {}
        output_units: Dict[str, str] = {}
        
        for i in range(self.n_samples):
            # Build input dict for this sample
            inputs = dict(fixed_params)
            for param in uncertain_params:
                value = param_samples[param.name][i]
                inputs[param.name] = ureg.Quantity(value, param.unit)
            
            try:
                # Run calculation
                outputs = calculation_func(**inputs)
                
                # Collect output values
                for name, quantity in outputs.items():
                    if name not in output_samples:
                        output_samples[name] = []
                        output_units[name] = str(quantity.units)
                    output_samples[name].append(quantity.magnitude)
            except Exception:
                # Skip failed calculations (will reduce effective sample size)
                continue
        
        # Calculate statistics for each output
        results = {}
        alpha = 1 - self.confidence_level
        
        for name, samples in output_samples.items():
            samples_arr = np.array(samples)
            
            results[name] = UncertaintyResult(
                parameter_name=name,
                unit=output_units[name],
                mean=float(np.mean(samples_arr)),
                std_dev=float(np.std(samples_arr, ddof=1)),
                median=float(np.median(samples_arr)),
                confidence_level=self.confidence_level,
                ci_lower=float(np.percentile(samples_arr, alpha/2 * 100)),
                ci_upper=float(np.percentile(samples_arr, (1 - alpha/2) * 100)),
                min_value=float(np.min(samples_arr)),
                max_value=float(np.max(samples_arr)),
                percentile_5=float(np.percentile(samples_arr, 5)),
                percentile_95=float(np.percentile(samples_arr, 95)),
                samples=samples_arr if keep_samples else None
            )
        
        return results
    
    def run_sensitivity_analysis(
        self,
        calculation_func: Callable[..., Dict[str, pint.Quantity]],
        uncertain_params: List[UncertainParameter],
        fixed_params: Optional[Dict[str, pint.Quantity]] = None
    ) -> Dict[str, SensitivityResult]:
        """
        Perform sensitivity analysis to identify dominant uncertainty contributors.
        
        Uses correlation-based sensitivity where sensitivity is the correlation
        coefficient between each input and output.
        
        Args:
            calculation_func: Function that takes parameter dict and returns output dict
            uncertain_params: List of uncertain input parameters
            fixed_params: Dictionary of fixed parameters
        
        Returns:
            Dictionary of SensitivityResult for each output
        """
        fixed_params = fixed_params or {}
        
        # Generate samples
        param_samples = {}
        for param in uncertain_params:
            param_samples[param.name] = param.generate_samples(self.n_samples, self.rng)
        
        # Run calculations and collect results
        input_arrays = {name: [] for name in param_samples}
        output_arrays: Dict[str, List[float]] = {}
        
        for i in range(self.n_samples):
            inputs = dict(fixed_params)
            for param in uncertain_params:
                value = param_samples[param.name][i]
                inputs[param.name] = ureg.Quantity(value, param.unit)
                input_arrays[param.name].append(value)
            
            try:
                outputs = calculation_func(**inputs)
                for name, quantity in outputs.items():
                    if name not in output_arrays:
                        output_arrays[name] = []
                    output_arrays[name].append(quantity.magnitude)
            except Exception:
                # Remove last input values if calculation failed
                for arr in input_arrays.values():
                    if arr:
                        arr.pop()
        
        # Calculate correlation coefficients
        results = {}
        
        for output_name, output_vals in output_arrays.items():
            output_arr = np.array(output_vals)
            sensitivities = {}
            
            for input_name, input_vals in input_arrays.items():
                input_arr = np.array(input_vals[:len(output_vals)])
                
                # Pearson correlation coefficient
                if len(input_arr) > 1:
                    corr = np.corrcoef(input_arr, output_arr)[0, 1]
                    sensitivities[input_name] = float(corr) if not np.isnan(corr) else 0.0
                else:
                    sensitivities[input_name] = 0.0
            
            results[output_name] = SensitivityResult(
                output_name=output_name,
                input_sensitivities=sensitivities
            )
        
        return results
    
    def quick_uncertainty(
        self,
        value: float,
        unit: str,
        relative_uncertainty: float = 0.05,
        n_samples: int = 1000
    ) -> UncertaintyResult:
        """
        Quick uncertainty estimate for a single value with assumed normal distribution.
        
        Args:
            value: Nominal value
            unit: Unit string
            relative_uncertainty: Relative uncertainty (default 5%)
            n_samples: Number of samples for this quick estimate
        
        Returns:
            UncertaintyResult
        """
        std_dev = value * relative_uncertainty
        samples = self.rng.normal(value, std_dev, n_samples)
        
        alpha = 1 - self.confidence_level
        
        return UncertaintyResult(
            parameter_name="value",
            unit=unit,
            mean=float(np.mean(samples)),
            std_dev=float(np.std(samples, ddof=1)),
            median=float(np.median(samples)),
            confidence_level=self.confidence_level,
            ci_lower=float(np.percentile(samples, alpha/2 * 100)),
            ci_upper=float(np.percentile(samples, (1 - alpha/2) * 100)),
            min_value=float(np.min(samples)),
            max_value=float(np.max(samples)),
            percentile_5=float(np.percentile(samples, 5)),
            percentile_95=float(np.percentile(samples, 95)),
            samples=samples
        )
