# LLE
# -----
# packages/modules
import time  # Import the time module
from typing import List, Dict, Optional, Literal, Any
import numpy as np
import pycuc
# local
from .equilibria import Equilibria
from .source import Source


class LLE(Equilibria):
    """
    LLE class for handling liquid-liquid equilibrium calculations.
    """

    def __init__(self,
                 components: List[str],
                 model_source: Optional[Dict] = None, **kwargs):
        """
        Initialize the LLE class.
        """
        self.__model_source = model_source
        self.__datasource = {
        } if model_source is None else model_source['datasource']
        self.__equationsource = {
        } if model_source is None else model_source['equationsource']

        # NOTE: init class
        Equilibria.__init__(self, components, model_source)

    @property
    def model_source(self) -> Dict:
        '''
        Get the model source property.

        Returns
        -------
        dict
            The model source dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__model_source is None:
            return {}
        return self.__model_source

    @property
    def datasource(self) -> Dict:
        '''
        Get the datasource property.

        Returns
        -------
        dict
            The datasource dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__datasource is None:
            return {}
        return self.__datasource

    @property
    def equationsource(self) -> Dict:
        '''
        Get the equationsource property.

        Returns
        -------
        dict
            The equationsource dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__equationsource is None:
            return {}
        return self.__equationsource

    def flash_isothermal(self):
        pass
