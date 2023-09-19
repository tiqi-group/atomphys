import pkg_resources
import pint

_ureg = pint.UnitRegistry(system="SI", auto_reduce_dimensions=True)
_ureg.enable_contexts('spectroscopy', 'Gaussian')
_ureg.define('_2pi = 2 * pi')
_ureg.default_format = "~0.3gP"

pint.set_application_registry(_ureg)
_ureg = pint.get_application_registry()

from .atom import Atom  # noqa: E402
# from .laser import Laser  # noqa: E402
from .state import State  # noqa: E402
from .transition import Transition  # noqa: E402
from .plot import plot_atom  # noqa: E402

__version__ = pkg_resources.get_distribution("atomphys").version

__all__ = [
    "Atom",
    "State",
    "Transition",
    "plot_atom"
]
