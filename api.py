"""
FastAPI REST API for Chemical Engineering Equation Solver.
Exposes all equations as RESTful endpoints.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Dict, Any, List, Optional
import uvicorn

# Import equation registries
from equations.process_control import PROCESS_CONTROL_EQUATIONS
from equations.fluid_dynamics import FLUID_DYNAMICS_EQUATIONS
from equations.heat_transfer import HEAT_TRANSFER_EQUATIONS
from equations.mass_transfer import MASS_TRANSFER_EQUATIONS
from equations.reaction_kinetics import REACTION_KINETICS_EQUATIONS
from equations.thermodynamics import THERMODYNAMICS_EQUATIONS
from equations.basic_math import BASIC_MATH_EQUATIONS
from equations.piping import PIPING_EQUATIONS
from equations.instrumentation import INSTRUMENTATION_EQUATIONS
from equations.distillation import DISTILLATION_EQUATIONS
from equations.pumps import PUMPS_EQUATIONS
from equations.vessels import VESSELS_EQUATIONS
from equations.safety import SAFETY_EQUATIONS
from equations.economics import ECONOMICS_EQUATIONS

# All equation categories
EQUATION_REGISTRY = {
    "process_control": PROCESS_CONTROL_EQUATIONS,
    "fluid_dynamics": FLUID_DYNAMICS_EQUATIONS,
    "piping": PIPING_EQUATIONS,
    "heat_transfer": HEAT_TRANSFER_EQUATIONS,
    "thermodynamics": THERMODYNAMICS_EQUATIONS,
    "mass_transfer": MASS_TRANSFER_EQUATIONS,
    "distillation": DISTILLATION_EQUATIONS,
    "reaction_kinetics": REACTION_KINETICS_EQUATIONS,
    "pumps": PUMPS_EQUATIONS,
    "vessels": VESSELS_EQUATIONS,
    "instrumentation": INSTRUMENTATION_EQUATIONS,
    "safety": SAFETY_EQUATIONS,
    "economics": ECONOMICS_EQUATIONS,
    "basic_math": BASIC_MATH_EQUATIONS,
}

# Create FastAPI app
app = FastAPI(
    title="Chemical Engineering Solver API",
    description="REST API for chemical engineering calculations",
    version="1.0.0"
)

# Enable CORS for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files directory for frontend
static_dir = Path(__file__).parent / "web" / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")


# Pydantic models for request/response
class CalculationRequest(BaseModel):
    inputs: Dict[str, Any]
    
    class Config:
        json_schema_extra = {
            "example": {
                "inputs": {
                    "Ku": 2.0,
                    "Pu": {"value": 5.0, "unit": "min"}
                }
            }
        }


class ParameterInfo(BaseModel):
    name: str
    description: str
    symbol: str
    unit: str
    param_type: str
    required: bool
    tooltip: Optional[str] = None
    typical_range: Optional[List[float]] = None


class EquationInfo(BaseModel):
    id: str
    name: str
    category: str
    description: str
    reference: str
    parameters: List[ParameterInfo]


class CalculationResult(BaseModel):
    success: bool
    outputs: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    validation_status: Optional[str] = None
    validation_messages: Optional[List[str]] = None


class CategoryInfo(BaseModel):
    name: str
    equation_count: int
    equations: List[Dict[str, str]]


# API Endpoints

@app.get("/")
async def root():
    """Serve the main web interface."""
    index_path = Path(__file__).parent / "web" / "index.html"
    if index_path.exists():
        return FileResponse(index_path)
    return {"message": "Chemical Engineering Solver API", "docs": "/docs"}


@app.get("/api/categories", response_model=List[CategoryInfo])
async def get_categories():
    """Get all available equation categories."""
    categories = []
    for cat_id, equations in EQUATION_REGISTRY.items():
        cat_name = cat_id.replace("_", " ").title()
        eq_list = []
        for eq_id, eq_class in equations.items():
            eq = eq_class()
            eq_list.append({"id": eq_id, "name": eq.name})
        categories.append(CategoryInfo(
            name=cat_name,
            equation_count=len(equations),
            equations=eq_list
        ))
    return categories


@app.get("/api/equations/{category}", response_model=List[Dict[str, str]])
async def get_equations_in_category(category: str):
    """Get all equations in a specific category."""
    if category not in EQUATION_REGISTRY:
        raise HTTPException(status_code=404, detail=f"Category '{category}' not found")
    
    equations = []
    for eq_id, eq_class in EQUATION_REGISTRY[category].items():
        eq = eq_class()
        equations.append({
            "id": eq_id,
            "name": eq.name,
            "description": eq.description
        })
    return equations


@app.get("/api/equations/{category}/{equation_id}", response_model=EquationInfo)
async def get_equation_details(category: str, equation_id: str):
    """Get detailed information about a specific equation."""
    if category not in EQUATION_REGISTRY:
        raise HTTPException(status_code=404, detail=f"Category '{category}' not found")
    
    if equation_id not in EQUATION_REGISTRY[category]:
        raise HTTPException(status_code=404, detail=f"Equation '{equation_id}' not found")
    
    eq_class = EQUATION_REGISTRY[category][equation_id]
    eq = eq_class()
    
    parameters = []
    for param in eq.parameters:
        parameters.append(ParameterInfo(
            name=param.name,
            description=param.description,
            symbol=param.symbol or param.name,
            unit=param.default_unit or "",
            param_type=param.param_type.value,
            required=param.required,
            tooltip=param.tooltip,
            typical_range=list(param.typical_range) if param.typical_range else None
        ))
    
    return EquationInfo(
        id=eq.equation_id,
        name=eq.name,
        category=eq.category,
        description=eq.description,
        reference=eq.reference,
        parameters=parameters
    )


@app.post("/api/calculate/{category}/{equation_id}", response_model=CalculationResult)
async def calculate(category: str, equation_id: str, request: CalculationRequest):
    """Execute a calculation with the given inputs."""
    if category not in EQUATION_REGISTRY:
        raise HTTPException(status_code=404, detail=f"Category '{category}' not found")
    
    if equation_id not in EQUATION_REGISTRY[category]:
        raise HTTPException(status_code=404, detail=f"Equation '{equation_id}' not found")
    
    eq_class = EQUATION_REGISTRY[category][equation_id]
    eq = eq_class()
    
    # Parse inputs - handle both simple values and {value, unit} objects
    parsed_inputs = {}
    for key, value in request.inputs.items():
        if isinstance(value, dict) and "value" in value:
            parsed_inputs[key] = (value["value"], value.get("unit", ""))
        else:
            parsed_inputs[key] = value
    
    try:
        result = eq.calculate(parsed_inputs, validate=True)
        
        # Format outputs
        outputs = {}
        for name, qty in result.outputs.items():
            if hasattr(qty, 'magnitude'):
                outputs[name] = {
                    "value": float(qty.magnitude),
                    "unit": str(qty.units) if not qty.dimensionless else ""
                }
            else:
                outputs[name] = {"value": float(qty), "unit": ""}
        
        # Format validation
        validation_status = None
        validation_messages = []
        if result.validation:
            validation_status = result.validation.overall_status.value
            validation_messages = [msg.message for msg in result.validation.messages]
        
        return CalculationResult(
            success=result.success,
            outputs=outputs,
            validation_status=validation_status,
            validation_messages=validation_messages
        )
    
    except Exception as e:
        return CalculationResult(
            success=False,
            error=str(e)
        )


@app.get("/api/search")
async def search_equations(q: str):
    """Search for equations by name or description."""
    results = []
    query = q.lower()
    
    for cat_id, equations in EQUATION_REGISTRY.items():
        for eq_id, eq_class in equations.items():
            eq = eq_class()
            if query in eq.name.lower() or query in eq.description.lower():
                results.append({
                    "category": cat_id,
                    "id": eq_id,
                    "name": eq.name,
                    "description": eq.description
                })
    
    return results


# Run server
if __name__ == "__main__":
    print("Starting Chemical Engineering Solver API...")
    print("API Documentation: http://localhost:8000/docs")
    print("Web Interface: http://localhost:8000")
    uvicorn.run(app, host="0.0.0.0", port=8000)
