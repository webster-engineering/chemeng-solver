"""
Learning Mode data structures for educational content.
Provides dataclasses for quiz questions, worked examples, and learning content.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from enum import Enum


class QuestionType(Enum):
    """Types of quiz questions."""
    MULTIPLE_CHOICE = "multiple_choice"
    NUMERIC = "numeric"
    TRUE_FALSE = "true_false"
    FILL_BLANK = "fill_blank"


class Difficulty(Enum):
    """Difficulty levels for equations and quizzes."""
    BEGINNER = "beginner"
    INTERMEDIATE = "intermediate"
    ADVANCED = "advanced"


@dataclass
class QuizQuestion:
    """A single quiz question for learning mode."""
    id: str
    question: str
    question_type: QuestionType
    correct_answer: str
    explanation: str
    difficulty: Difficulty = Difficulty.BEGINNER
    options: List[str] = field(default_factory=list)  # For multiple choice
    hint: str = ""
    points: int = 10
    
    def validate_answer(self, user_answer: str) -> bool:
        """Check if the user's answer is correct."""
        if self.question_type == QuestionType.NUMERIC:
            try:
                user_val = float(user_answer)
                correct_val = float(self.correct_answer)
                # Allow 1% tolerance for numeric answers
                return abs(user_val - correct_val) / max(abs(correct_val), 1e-10) < 0.01
            except ValueError:
                return False
        elif self.question_type == QuestionType.TRUE_FALSE:
            return user_answer.lower().strip() == self.correct_answer.lower().strip()
        else:
            return user_answer.strip().lower() == self.correct_answer.strip().lower()


@dataclass
class CalculationStep:
    """A single step in a worked example or explanation."""
    step_number: int
    title: str
    description: str
    formula: str = ""  # LaTeX or plain text formula
    substitution: str = ""  # Values substituted into formula
    computation: str = ""  # The actual calculation
    result: str = ""
    notes: str = ""  # Additional tips or warnings


@dataclass
class WorkedExample:
    """A complete worked example with step-by-step solution."""
    title: str
    scenario: str
    given_values: Dict[str, str]  # param_name -> "value unit"
    find: List[str]  # What we're solving for
    steps: List[CalculationStep] = field(default_factory=list)
    final_answer: str = ""
    real_world_context: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "title": self.title,
            "scenario": self.scenario,
            "given_values": self.given_values,
            "find": self.find,
            "steps": [
                {
                    "step_number": s.step_number,
                    "title": s.title,
                    "description": s.description,
                    "formula": s.formula,
                    "substitution": s.substitution,
                    "computation": s.computation,
                    "result": s.result,
                    "notes": s.notes
                }
                for s in self.steps
            ],
            "final_answer": self.final_answer,
            "real_world_context": self.real_world_context
        }


@dataclass
class LearningContent:
    """Complete learning content for an equation."""
    # Theory and background
    background_theory: str = ""
    key_concepts: List[str] = field(default_factory=list)
    
    # Real-world applications
    real_world_applications: List[str] = field(default_factory=list)
    industry_examples: List[str] = field(default_factory=list)
    
    # Common mistakes and tips
    common_mistakes: List[str] = field(default_factory=list)
    pro_tips: List[str] = field(default_factory=list)
    
    # New enhanced fields for Introduction section
    references: List[str] = field(default_factory=list)  # Academic/industry sources
    derivation_summary: str = ""  # Brief explanation of equation origin
    limitations_assumptions: List[str] = field(default_factory=list)  # When it applies/doesn't
    
    # Variable information
    variable_sources: Dict[str, str] = field(default_factory=dict)  # param -> where to find
    
    # Quiz and assessment
    quiz_questions: List[QuizQuestion] = field(default_factory=list)
    
    # Worked example
    worked_example: Optional[WorkedExample] = None
    
    # Metadata
    difficulty: Difficulty = Difficulty.BEGINNER
    estimated_time_minutes: int = 10
    prerequisites: List[str] = field(default_factory=list)  # equation_ids
    related_equations: List[str] = field(default_factory=list)
    
    # Diagram information
    diagram_type: Optional[str] = None  # 'gas_tank', 'heat_exchanger', etc.
    diagram_labels: Dict[str, str] = field(default_factory=dict)  # label_id -> description
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "background_theory": self.background_theory,
            "key_concepts": self.key_concepts,
            "real_world_applications": self.real_world_applications,
            "industry_examples": self.industry_examples,
            "common_mistakes": self.common_mistakes,
            "pro_tips": self.pro_tips,
            "references": self.references,
            "derivation_summary": self.derivation_summary,
            "limitations_assumptions": self.limitations_assumptions,
            "variable_sources": self.variable_sources,
            "quiz_questions": [
                {
                    "id": q.id,
                    "question": q.question,
                    "question_type": q.question_type.value,
                    "correct_answer": q.correct_answer,
                    "explanation": q.explanation,
                    "difficulty": q.difficulty.value,
                    "options": q.options,
                    "hint": q.hint,
                    "points": q.points
                }
                for q in self.quiz_questions
            ],
            "worked_example": self.worked_example.to_dict() if self.worked_example else None,
            "difficulty": self.difficulty.value,
            "estimated_time_minutes": self.estimated_time_minutes,
            "prerequisites": self.prerequisites,
            "related_equations": self.related_equations,
            "diagram_type": self.diagram_type,
            "diagram_labels": self.diagram_labels
        }


@dataclass
class UserProgress:
    """Track user's learning progress for an equation."""
    equation_id: str
    completed_steps: List[int] = field(default_factory=list)  # Step numbers completed
    quiz_scores: Dict[str, bool] = field(default_factory=dict)  # question_id -> correct
    attempts: int = 0
    best_score: int = 0
    completed: bool = False
    last_accessed: str = ""  # ISO timestamp
    
    @property
    def quiz_percentage(self) -> float:
        """Calculate quiz score percentage."""
        if not self.quiz_scores:
            return 0.0
        correct = sum(1 for v in self.quiz_scores.values() if v)
        return (correct / len(self.quiz_scores)) * 100
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "equation_id": self.equation_id,
            "completed_steps": self.completed_steps,
            "quiz_scores": self.quiz_scores,
            "attempts": self.attempts,
            "best_score": self.best_score,
            "completed": self.completed,
            "last_accessed": self.last_accessed,
            "quiz_percentage": self.quiz_percentage
        }


@dataclass
class LearningStats:
    """Overall learning statistics for gamification."""
    total_equations_completed: int = 0
    total_quizzes_passed: int = 0
    current_streak: int = 0
    longest_streak: int = 0
    total_xp: int = 0
    badges: List[str] = field(default_factory=list)
    last_activity_date: str = ""
    
    # XP rewards
    XP_PER_EQUATION = 100
    XP_PER_QUIZ = 50
    XP_STREAK_BONUS = 25
    
    def add_xp(self, amount: int) -> None:
        """Add XP with streak bonus."""
        bonus = self.XP_STREAK_BONUS if self.current_streak > 0 else 0
        self.total_xp += amount + bonus
    
    def check_badges(self) -> List[str]:
        """Check and award new badges."""
        new_badges = []
        
        badge_criteria = {
            "first_calculation": self.total_equations_completed >= 1,
            "equation_explorer": self.total_equations_completed >= 10,
            "quiz_master": self.total_quizzes_passed >= 10,
            "week_warrior": self.longest_streak >= 7,
            "century_club": self.total_xp >= 1000,
        }
        
        for badge_id, earned in badge_criteria.items():
            if earned and badge_id not in self.badges:
                self.badges.append(badge_id)
                new_badges.append(badge_id)
        
        return new_badges
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "total_equations_completed": self.total_equations_completed,
            "total_quizzes_passed": self.total_quizzes_passed,
            "current_streak": self.current_streak,
            "longest_streak": self.longest_streak,
            "total_xp": self.total_xp,
            "badges": self.badges,
            "last_activity_date": self.last_activity_date
        }
