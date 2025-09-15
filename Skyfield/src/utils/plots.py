import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
from matplotlib.ticker import FuncFormatter

def set_plot_style():
    """
    Set the global plot style for matplotlib.
    This function configures the default aesthetics for plots.
    """
    plt.rcParams['axes.grid'] = False  # Disable grid lines
    plt.rcParams['lines.linewidth'] = 2.5  # Set line width
    plt.rcParams['grid.alpha'] = 0.5  # Set grid transparency
    plt.rcParams['savefig.dpi'] = 300  # High resolution for saved figures
    # Use sans-serif fonts which are more reliably available
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Serif"]
    plt.rcParams["legend.loc"] = 'lower center'  # Default legend location
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["font.size"] = 11
    # Set default colormap to plasma
    plt.rcParams["image.cmap"] = "plasma"

def ra_formatter(x, pos):
    hours = int(x)
    minutes = int((x - hours) * 60)
    seconds = int(((x - hours) * 60 - minutes) * 60)
    return f'{hours:02d}:{minutes:02d}:{seconds:02d}'