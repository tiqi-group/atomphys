from .atom import Atom
from .state import State
from .transition import Transition
from .electric_field import ElectricField, GaussianBeam
from .data_utils.load_data import from_json, from_nist
from importlib.metadata import version
import pint

_ureg = pint.UnitRegistry(system="SI", auto_reduce_dimensions=True)
_ureg.enable_contexts("spectroscopy", "Gaussian")
_ureg.define("_2pi = 2 * pi")
_ureg.default_format = "~0.3gP"

pint.set_application_registry(_ureg)
_ureg = pint.get_application_registry()
__version__ = version("atomphys")
__all__ = [
    "Atom",
    "State",
    "Transition",
    "ElectricField",
    "GaussianBeam",
    "from_json",
    "from_nist",
]
