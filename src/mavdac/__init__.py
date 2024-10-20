from .mavdac import *  # type: ignore # noqa: F403
from .util import coeffs_from_cogs, run_mavdac

__doc__ = mavdac.__doc__  # type: ignore # noqa: F405
__all__ = [
    "coeffs_from_cogs",
    "run_mavdac",
]
if hasattr(mavdac, "__all__"):  # type: ignore # noqa: F405
    __all__ += mavdac.__all__  # type: ignore # noqa: F405
