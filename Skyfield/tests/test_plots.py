import ast
from pathlib import Path
import runpy
import subprocess
import sys

import numpy as np
import pytest

from ra_formatting import ra_formatter


@pytest.fixture(scope='module')
def plots_module():
    path = Path(__file__).parents[1] / 'src' / 'utils' / 'plots.py'
    print(f'Reading from {path}...')
    return runpy.run_path(str(path))


@pytest.mark.parametrize(
    ('hours', 'expected'),
    [
        (0.0, '00:00:00'),
        (1.25, '01:15:00'),
        (1 + 1 / 60, '01:01:00'),
        (2 + 1 / 60, '02:01:00'),
        (1 + 59.8 / 60, '01:59:48'),
        (23 + 59 / 60 + 59.8 / 3600, '24:00:00'),
        (24.0, '24:00:00'),
        (-0.5, '-00:30:00'),
    ],
)
def test_ra_formatter_handles_sexagesimal_boundaries(plots_module, hours, expected):
    assert plots_module['ra_formatter'](hours, 0) == expected


@pytest.mark.parametrize('direction', [-np.inf, np.inf])
def test_ra_formatter_is_stable_around_an_exact_minute(plots_module, direction):
    hours = np.nextafter(1 + 1 / 60, direction)
    assert plots_module['ra_formatter'](hours, 0) == '01:01:00'


def test_existing_plots_utility_reexports_the_lightweight_formatter(plots_module):
    assert plots_module['ra_formatter'] is ra_formatter


@pytest.mark.parametrize(
    'path',
    sorted(Path(__file__).parents[1].glob('skyfield_*.py')) + [
        Path(__file__).parents[1] / 'targets_master.py',
        Path(__file__).parents[1] / 'targets_random.py',
    ],
    ids=lambda path: path.name,
)
def test_historical_script_callback_imports_preserve_shared_identity(path):
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    assert not any(isinstance(node, ast.FunctionDef) and node.name == 'ra_formatter' for node in tree.body)
    imports = [node for node in tree.body if isinstance(node, ast.ImportFrom) and node.module == 'ra_formatting']
    assert len(imports) == 1
    namespace = {}

    exec(compile(ast.Module(body=imports, type_ignores=[]), str(path), 'exec'), namespace)

    assert namespace['ra_formatter'] is ra_formatter
    assert namespace['ra_formatter'](1 + 1 / 60, 0) == '01:01:00'


def test_formatter_renders_correct_matplotlib_tick_labels():
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.ticker import FixedLocator, FuncFormatter

    figure = Figure(figsize=(6, 3))
    canvas = FigureCanvasAgg(figure)
    axes = figure.subplots()
    axes.plot([0.0, 1.0, 2.1], [0.0, 1.0, 0.0])
    axes.xaxis.set_major_locator(FixedLocator([0.0, 1 + 1 / 60, 2 + 1 / 60]))
    axes.xaxis.set_major_formatter(FuncFormatter(ra_formatter))

    canvas.draw()

    assert [label.get_text() for label in axes.get_xticklabels()] == ['00:00:00', '01:01:00', '02:01:00']
    assert np.any(np.asarray(canvas.buffer_rgba())[:, :, :3] != 255)


def test_lightweight_import_does_not_initialize_satellite_analysis():
    result = subprocess.run(
        [sys.executable, '-c',
         'import sys; from ra_formatting import ra_formatter; '
         'assert not {"src", "skyfield", "lumos"}.intersection(sys.modules); '
         'print(ra_formatter(1 + 1 / 60, 0))'],
        cwd=Path(__file__).parents[1], check=True, text=True, capture_output=True,
    )
    assert result.stdout == '01:01:00\n'
