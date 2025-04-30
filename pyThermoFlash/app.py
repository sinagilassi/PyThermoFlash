# import libs
from typing import List, Dict, Optional
# local
from .docs import VLE
from .utils import model_source_checker


def vle(components: List[str],
        model_source: Optional[Dict] = None,
        **kwargs) -> VLE:
    '''
    VLE model for vapor-liquid equilibrium (VLE) calculations for
    multi-component systems.

    Parameters
    ----------
    components : list
        List of component names.
    model_source : dict, optional
        Model source parameters defined as a dictionary.
    kwargs : dict
        Additional parameters for the model.

    Returns
    -------
    VLE
        VLE model object.

    Notes
    -----
    The model source can be a dictionary containing the following keys:

    - datasource: str
        The data source for the model.
    - equationsource: str
        The equation source for the model.

    These two sources are generated with `PythermoDBLink`, please refer to the
    documentation for more details.

    ```python
    # model source example
    model_source = {
        "datasource": datasource,
        "equationsource": equationsource
    }
    ```
    '''
    try:
        # NOTE: check if components are valid
        if not isinstance(components, list) or not all(isinstance(c, str) for c in components):
            raise ValueError("Components must be a list of strings.")

        # NOTE: check if model_source is valid
        if model_source is not None:
            if not model_source_checker(model_source):
                raise ValueError("Invalid model source.")

        # NOTE
        # init
        return VLE(components, model_source, **kwargs)
    except Exception as e:
        raise RuntimeError(f"Failed to create VLE model: {e}") from e
