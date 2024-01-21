# TODO: remove duplicated function between the modules

from ._maximize_buckling import maximize_buckling_load
from ._minimize_panel_weight import minimize_panel_weight


__all__ = ['maximize_buckling_load', 'minimize_panel_weight']