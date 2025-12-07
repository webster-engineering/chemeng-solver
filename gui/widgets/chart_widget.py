"""
Chart widget for plotting response curves and distributions.
"""

from PySide6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np


class ChartWidget(QWidget):
    """Widget for displaying matplotlib charts."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.figure = Figure(figsize=(6, 4), dpi=100, facecolor='#1e1e2e')
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        
        self.ax = self.figure.add_subplot(111)
        self._style_axis()
    
    def _style_axis(self):
        """Apply dark theme styling to axis."""
        self.ax.set_facecolor('#1e1e2e')
        self.ax.tick_params(colors='#cdd6f4', which='both')
        self.ax.xaxis.label.set_color('#cdd6f4')
        self.ax.yaxis.label.set_color('#cdd6f4')
        self.ax.title.set_color('#89b4fa')
        
        for spine in self.ax.spines.values():
            spine.set_color('#45475a')
        
        self.ax.grid(True, color='#313244', alpha=0.5)
    
    def plot_step_response(self, K, tau, zeta=1.0, theta=0, t_max=None):
        """Plot step response for first or second order system."""
        self.ax.clear()
        self._style_axis()
        
        if t_max is None:
            t_max = 5 * tau * max(1, zeta)
        
        t = np.linspace(0, t_max + theta, 500)
        y = np.zeros_like(t)
        
        mask = t >= theta
        t_shifted = t[mask] - theta
        
        if zeta >= 1:  # First order or overdamped
            y[mask] = K * (1 - np.exp(-t_shifted / tau))
        else:  # Underdamped
            wd = np.sqrt(1 - zeta**2) / tau
            y[mask] = K * (1 - np.exp(-zeta * t_shifted / tau) * 
                          (np.cos(wd * t_shifted) + 
                           zeta / np.sqrt(1 - zeta**2) * np.sin(wd * t_shifted)))
        
        self.ax.plot(t, y, color='#89b4fa', linewidth=2)
        self.ax.axhline(y=K, color='#a6adc8', linestyle='--', alpha=0.5, label='Final Value')
        
        if theta > 0:
            self.ax.axvline(x=theta, color='#f9e2af', linestyle=':', alpha=0.5, label='Dead Time')
        
        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Response')
        self.ax.set_title('Step Response')
        self.ax.legend(facecolor='#313244', edgecolor='#45475a', labelcolor='#cdd6f4')
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def plot_histogram(self, data, xlabel='Value', title='Distribution'):
        """Plot histogram for uncertainty distribution."""
        self.ax.clear()
        self._style_axis()
        
        n, bins, patches = self.ax.hist(data, bins=50, color='#89b4fa', 
                                         alpha=0.7, edgecolor='#1e1e2e')
        
        mean = np.mean(data)
        self.ax.axvline(x=mean, color='#f38ba8', linestyle='-', linewidth=2, label=f'Mean: {mean:.4g}')
        
        percentile_5 = np.percentile(data, 5)
        percentile_95 = np.percentile(data, 95)
        self.ax.axvline(x=percentile_5, color='#f9e2af', linestyle='--', label='5th %ile')
        self.ax.axvline(x=percentile_95, color='#f9e2af', linestyle='--', label='95th %ile')
        
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel('Frequency')
        self.ax.set_title(title)
        self.ax.legend(facecolor='#313244', edgecolor='#45475a', labelcolor='#cdd6f4')
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def plot_sensitivity(self, sensitivities: dict, title='Sensitivity Analysis'):
        """Plot tornado chart for sensitivity analysis."""
        self.ax.clear()
        self._style_axis()
        
        sorted_items = sorted(sensitivities.items(), key=lambda x: abs(x[1]))
        names = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]
        
        colors = ['#a6e3a1' if v >= 0 else '#f38ba8' for v in values]
        
        y_pos = np.arange(len(names))
        self.ax.barh(y_pos, values, color=colors, alpha=0.8)
        self.ax.set_yticks(y_pos)
        self.ax.set_yticklabels(names)
        self.ax.axvline(x=0, color='#cdd6f4', linewidth=0.5)
        
        self.ax.set_xlabel('Correlation Coefficient')
        self.ax.set_title(title)
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def clear(self):
        """Clear the chart."""
        self.ax.clear()
        self._style_axis()
        self.canvas.draw()
