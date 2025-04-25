# import libs
from typing import List, Dict, Optional
# local
from .docs import VLE


def vle(components: List[str], model_source: Optional[Dict] = None, **kwargs) -> VLE:
    '''
    VLE model for vapor-liquid equilibrium (VLE) calculations for multi-component systems.

    Parameters
    ----------
    components : list
        List of component names.
    model_source : dict
        Model source parameters.
    kwargs : dict
        Additional parameters for the model.

    Returns
    -------
    VLE
        VLE model object.
    '''
    # init
    return VLE(components, model_source, **kwargs)
