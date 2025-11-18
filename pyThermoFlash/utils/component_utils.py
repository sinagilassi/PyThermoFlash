# import libs
import logging
from typing import (
    List,
    Literal,
    Dict
)
from pythermodb_settings.models import (
    Component,
)
# locals

# NOTE: logger
logger = logging.getLogger(__name__)


def set_feed_specification(
    components: List[Component],
    component_key: str = 'Name-State',
) -> Dict[str, float]:
    """
    Set feed specification for a list of components with their mole fractions.

    Parameters
    ----------
    components : List[Component]
        List of Component objects, each with name, formula, state, and mole_fraction attributes.
    component_key : str, optional
        Key to use for feed specification. Options are 'Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State' or 'Formula-Name-State'.

    Returns
    -------
    Dict[str, float]
        Dictionary with component identifiers as keys and their mole fractions as values.
    """
    try:
        # NOTE: Initialize feed specification dictionary
        feed_spec = {}

        # NOTE: Iterate over components to set feed specification
        for i, component in enumerate(components):
            # set
            name_ = component.name
            formula_ = component.formula
            state_ = component.state

            # Check if mole_fraction is provided, otherwise skip
            if component.mole_fraction is None:
                logging.warning(
                    f"Component {name_} does not have a mole fraction defined. Skipping.")
                continue

            # std component key
            component_key = component_key.strip().lower()

            # NOTE: Set feed specification
            if component_key == 'Name-State'.lower():
                feed_spec[f"{name_}-{state_}"] = component.mole_fraction
            elif component_key == 'Formula-State'.lower():
                feed_spec[f"{formula_}-{state_}"] = component.mole_fraction
            elif component_key == 'Name'.lower():
                feed_spec[name_] = component.mole_fraction
            elif component_key == 'Formula'.lower():
                feed_spec[formula_] = component.mole_fraction
            elif component_key == 'Name-Formula-State'.lower():
                feed_spec[f"{name_}-{formula_}-{state_}"] = component.mole_fraction
            elif component_key == 'Formula-Name-State'.lower():
                feed_spec[f"{formula_}-{name_}-{state_}"] = component.mole_fraction
            else:
                # raise ValueError("Invalid component_key. Use 'name' or 'formula'.")
                logging.error(
                    f"Invalid component_key: {component_key}. Use 'name' or 'formula'.")
                raise ValueError(
                    f"Invalid component_key: {component_key}. Use 'name' or 'formula'.")

        return feed_spec
    except Exception as e:
        logging.error(f"Failed to set feed specification: {e}")
        raise Exception(f"Failed to set feed specification: {e}") from e
