# SOURCE
# import libs
from typing import List, Dict, Optional
# local
from ..config import DATASOURCE, EQUATIONSOURCE


class Source:
    '''
    Source class for the pyThermoFlash package.
    '''
    # NOTE: variables

    def __init__(self, model_source: Optional[Dict] = None, **kwargs):
        '''Initialize the Source class.'''
        # set
        self.model_source = model_source

        # NOTE: source
        # set
        self._datasource, self._equationsource = self.set_source(
            model_source=model_source)

    def __repr__(self):
        des = "Source class for the pyThermoFlash package."
        return des

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
        if self._datasource is None:
            return {}
        return self._datasource

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
        if self._equationsource is None:
            return {}
        return self._equationsource

    def set_source(self, model_source: Dict):
        '''
        Set the model source.

        Parameters
        ----------
        model_source : dict
            The model source dictionary.
        '''
        # NOTE: source
        # datasource
        _datasource = {
        } if model_source is None else model_source[DATASOURCE]

        # equationsource
        _equationsource = {
        } if model_source is None else model_source[EQUATIONSOURCE]

        # res
        return _datasource, _equationsource

    def eq_extractor(self, component_name: str, prop_name: str):
        '''
        Extracts the equilibrium property from the model source.

        Parameters
        ----------
        component_name : str
            The name of the component.
        prop_name : str
            The name of the property to extract.

        Returns
        -------
        dict
            The extracted property.
        '''
        if self.equationsource is None:
            raise ValueError("Equation source is not defined.")

        # NOTE: check component
        if component_name not in self.equationsource.keys():
            raise ValueError(
                f"Component '{component_name}' not found in model source.")

        # NOTE: check property
        if prop_name not in self.equationsource[component_name].keys():
            raise ValueError(
                f"Property '{prop_name}' not found in model source registered for {component_name}.")

        return self.equationsource[component_name][prop_name]

    def data_extractor(self, component_name: str, prop_name: str):
        '''
        Extracts the data property from the model source.

        Parameters
        ----------
        component_name : str
            The name of the component.
        prop_name : str
            The name of the property to extract.

        Returns
        -------
        dict
            The extracted property.
        '''
        if self.datasource is None:
            return None

        # NOTE: check component
        if component_name not in self.datasource.keys():
            raise ValueError(
                f"Component '{component_name}' not found in model datasource.")

        # NOTE: check property
        if prop_name not in self.datasource[component_name].keys():
            raise ValueError(
                f"Property '{prop_name}' not found in model datasource registered for {component_name}.")

        return self.datasource[component_name][prop_name]

    def check_args(self, component_name: str, args):
        '''
        Checks equation args

        Parameters
        ----------
        component_name : str
            The name of the component.
        args : tuple
            equation args
        '''
        try:
            # required args
            required_args = []

            # datasource list
            datasource_component_list = list(
                self.datasource[component_name].keys())

            # NOTE: default args
            datasource_component_list.append("P")
            datasource_component_list.append("T")

            # check args within datasource
            for arg_key, arg_value in args.items():
                # symbol
                if arg_value['symbol'] in datasource_component_list:
                    # update
                    required_args.append(arg_value)
                else:
                    raise Exception('Args not in datasource!')

            # res
            return required_args

        except Exception as e:
            raise Exception('Finding args failed!, ', e)

    def build_args(self, component_name: str,
                   args,
                   ignore_symbols: List[str] = ["T", "P"]):
        '''
        Builds args

        Parameters
        ----------
        component_name : str
            The name of the component.
        args : tuple
            equation args
        ignore_symbols : list
            list of symbols to ignore, default is ["T", "P"]
        '''
        try:
            # res
            res = {}
            for arg in args:
                # symbol
                symbol = arg['symbol']

                # NOTE: check if symbol is in ignore symbols
                if ignore_symbols is not None:
                    # check in ignore symbols
                    if symbol not in ignore_symbols:
                        # check in component database
                        for key, value in self.datasource.items():
                            if symbol == key:
                                res[symbol] = value
                else:
                    # check in component database
                    for key, value in self.datasource.items():
                        if symbol == key:
                            res[symbol] = value
            return res
        except Exception as e:
            raise Exception('Building args failed!, ', e)
