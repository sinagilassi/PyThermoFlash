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
        self.__model_source = model_source

        # NOTE: source
        # set
        self.__datasource = {
        } if model_source is None else model_source[DATASOURCE]
        self.__equationsource = {
        } if model_source is None else model_source[EQUATIONSOURCE]

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
        if self.__model_source is None:
            raise ValueError("Model source is not defined.")

        if prop_name not in self.__model_source[EQUATIONSOURCE]:
            raise ValueError(
                f"Property '{prop_name}' not found in model source.")

        return self.__model_source[EQUATIONSOURCE][component_name][prop_name]

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

        if prop_name not in self.datasource.keys():
            raise ValueError(
                f"Property '{prop_name}' not found in model datasource.")

        return self.datasource[component_name][prop_name]

    def check_args(self, args):
        '''
        Checks equation args

        Parameters
        ----------
        args : tuple
            equation args
        '''
        try:
            # required args
            required_args = []

            # datasource list
            datasource_component_list = list(self.datasource.keys())
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

    def build_args(self, args, ignore_symbols: List[str] = ["T", "P"]):
        '''
        Builds args

        Parameters
        ----------
        args : tuple
            equation args
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
