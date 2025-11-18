# ACTIVITY
# ---------
# import libs
import numpy as np
import pyThermoModels as ptm
from typing import List, Dict, Optional, List, Any
from pyThermoDB import TableMatrixData
# local


class Activity:
    def __init__(
        self,
        datasource: Optional[Dict[str, Any]] = None,
        equationsource: Optional[Dict[str, Any]] = None
    ):
        '''
        Activity model class for pyThermoFlash.

        Parameters
        ----------
        datasource : dict, optional
            Data source for the model, containing data for components and equations.
        equationsource : dict, optional
            Equation source for the model, containing equations for components.
        '''
        # set datasource
        self.datasource = datasource
        # set equationsource
        self.equationsource = equationsource

    def __repr__(self):
        return "Activity model class for pyThermoFlash"

    def NRTL(
        self,
        components: List[str],
        z_i_comp: Dict[str, float],
        temperature: float,
        **kwargs
    ):
        '''
        NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        components : list
            List of component names.
        z_i_comp : dict
            Dictionary of component names and their respective mole fractions.
        temperature : float
            Temperature in Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: check src
            # extract activity model inputs
            activity_inputs = kwargs.get('activity_inputs', None)
            nrtl_datasource = kwargs.get('NRTL', None)

            # check
            if activity_inputs is None:
                # ! check nrtl inputs in datasource
                activity_inputs = nrtl_datasource

            # check if activity_inputs is None
            if activity_inputs is None:
                raise ValueError(
                    "No valid source provided for activity model (NRTL) inputs.")

            # SECTION: check if activity_inputs is a dictionary
            if activity_inputs is not None:
                # check if activity_inputs is a dictionary
                if not isinstance(activity_inputs, dict):
                    raise ValueError(
                        "activity_inputs must be a dictionary.")
                # check if activity_inputs is empty
                if len(activity_inputs) == 0:
                    raise ValueError(
                        "activity_inputs cannot be empty.")

                # NOTE: update activity_inputs with nrtl_datasource
                if nrtl_datasource is not None:
                    activity_inputs.update(nrtl_datasource)

            # ! set initial values
            # a_ij = None
            # b_ij = None
            # c_ij = None
            # d_ij = None
            # dg_ij = None
            # alpha_ij = None
            # tau_ij = None

            # SECTION: extract activity model inputs
            # NOTE: method 1
            # ! Δg_ij, interaction energy parameter
            dg_ij_src = activity_inputs.get('dg_ij', None)
            if dg_ij_src is None:
                dg_ij_src = activity_inputs.get('dg', None)

            # NOTE: method 2
            # ! constants a, b, c, and d
            a_ij_src = activity_inputs.get('a_ij', None)
            if a_ij_src is None:
                a_ij_src = activity_inputs.get('a', None)
            b_ij_src = activity_inputs.get('b_ij', None)
            if b_ij_src is None:
                b_ij_src = activity_inputs.get('b', None)
            c_ij_src = activity_inputs.get('c_ij', None)
            if c_ij_src is None:
                c_ij_src = activity_inputs.get('c', None)
            d_ij_src = activity_inputs.get('d_ij', None)
            if d_ij_src is None:
                d_ij_src = activity_inputs.get('d', None)

            # NOTE: α_ij, non-randomness parameter
            alpha_ij_src = activity_inputs.get('alpha_ij', None)
            if alpha_ij_src is None:
                alpha_ij_src = activity_inputs.get('alpha', None)

            # SECTION: init NRTL model
            # activity model
            activity = ptm.activity(
                components=components, model_name='NRTL')
            # set
            activity_nrtl = activity.nrtl
            # check
            if activity_nrtl is None:
                raise ValueError(
                    "Failed to initialize NRTL activity model. Please check the components and model name.")

            # NOTE: check method
            tau_ij_cal_method = 0
            if dg_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if a_ij_src is None or b_ij_src is None or c_ij_src is None or d_ij_src is None:
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dg_ij_src is not None:
                # use dg_ij
                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', components)
                elif isinstance(dg_ij_src, list):
                    dg_ij = np.array(dg_ij_src)
                elif isinstance(dg_ij_src, np.ndarray):
                    dg_ij = dg_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: extract data
            # α_ij, non-randomness parameter
            if isinstance(alpha_ij_src, TableMatrixData):
                alpha_ij = alpha_ij_src.mat('alpha', components)
            elif isinstance(alpha_ij_src, list):
                alpha_ij = np.array(alpha_ij_src)
            elif isinstance(alpha_ij_src, np.ndarray):
                alpha_ij = alpha_ij_src
            else:
                raise ValueError(
                    "Invalid source for non-randomness parameter (α_ij). Must be TableMatrixData, list of lists, or numpy array.")

            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            if tau_ij_cal_method == 1:
                # check
                tau_ij, _ = activity_nrtl.cal_tau_ij_M1(
                    temperature=temperature,
                    dg_ij=dg_ij
                )
            elif tau_ij_cal_method == 2:
                tau_ij, _ = activity_nrtl.cal_tau_ij_M2(
                    temperature=temperature,
                    a_ij=a_ij,
                    b_ij=b_ij,
                    c_ij=c_ij,
                    d_ij=d_ij
                )
            else:
                raise ValueError(
                    "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # NOTE: nrtl inputs
            inputs_ = {
                'mole_fraction': z_i_comp,
                "tau_ij": tau_ij,
                "alpha_ij": alpha_ij
            }

            # NOTE: calculate activity
            res_, _ = activity_nrtl.cal(model_input=inputs_)

            # res format
            # res = {
            #     'property_name': 'activity coefficients',
            #     'components': components,
            #     'mole_fraction': xi,
            #     'value': AcCo_i,
            #     'unit': 1,
            #     'symbol': "AcCo_i",
            #     'message': message,
            # }

            # res
            return res_
        except Exception as e:
            raise Exception(f"Failed to calculate NRTL activity: {e}") from e

    def UNIQUAC(
        self,
            components: List[str],
            z_i_comp: Dict[str, float],
            temperature: float,
            **kwargs
    ):
        '''
        UNIQUAC activity model for calculating activity coefficients.

        Parameters
        ----------
        components : list
            List of component names.
        z_i_comp : dict
            Dictionary of component names and their respective mole fractions.
        temperature : float
            Temperature in Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: check src
            # extract activity model inputs
            activity_inputs = kwargs.get('activity_inputs', None)
            uniquac_datasource = kwargs.get('UNIQUAC', None)

            # check
            if activity_inputs is None:
                # ! check nrtl inputs in datasource
                activity_inputs = uniquac_datasource

            # check if activity_inputs is None
            if activity_inputs is None:
                raise ValueError(
                    "No valid source provided for activity model (UNIQUAC) inputs.")

            # NOTE: check if activity_inputs is a dictionary
            if activity_inputs is not None:
                # check if activity_inputs is a dictionary
                if not isinstance(activity_inputs, dict):
                    raise ValueError(
                        "activity_inputs must be a dictionary.")
                # check if activity_inputs is empty
                if len(activity_inputs) == 0:
                    raise ValueError(
                        "activity_inputs cannot be empty.")

            # NOTE: update activity_inputs with uniquac_datasource
            if uniquac_datasource is not None:
                activity_inputs.update(uniquac_datasource)

            # NOTE: method 1
            # Δg_ij, interaction energy parameter
            dU_ij_src = activity_inputs.get('dU_ij', None)
            if dU_ij_src is None:
                dU_ij_src = activity_inputs.get('dU', None)

            # NOTE: method 2
            # constants a, b, c, and d
            a_ij_src = activity_inputs.get('a_ij', None)
            if a_ij_src is None:
                a_ij_src = activity_inputs.get('a', None)
            b_ij_src = activity_inputs.get('b_ij', None)
            if b_ij_src is None:
                b_ij_src = activity_inputs.get('b', None)
            c_ij_src = activity_inputs.get('c_ij', None)
            if c_ij_src is None:
                c_ij_src = activity_inputs.get('c', None)
            d_ij_src = activity_inputs.get('d_ij', None)
            if d_ij_src is None:
                d_ij_src = activity_inputs.get('d', None)

            # NOTE: r_i, relative van der Waals volume of component i
            r_i_src = activity_inputs.get('r_i', None)
            if r_i_src is None:
                r_i_src = activity_inputs.get('r', None)

            # final check if r_i is provided
            if r_i_src is None:
                raise ValueError("No valid source provided for r_i.")

            # check if r_i is a list or numpy array
            if isinstance(r_i_src, list):
                r_i = np.array(r_i_src)
            elif isinstance(r_i_src, np.ndarray):
                r_i = r_i_src
            else:
                raise ValueError(
                    "Invalid source for r_i. Must be a list or numpy array.")

            # NOTE: q_i, relative van der Waals area of component i
            q_i_src = activity_inputs.get('q_i', None)
            if q_i_src is None:
                q_i_src = activity_inputs.get('q', None)

            # final check if q_i is provided
            if q_i_src is None:
                raise ValueError("No valid source provided for q_i.")

            # check if q_i is a list or numpy array
            if isinstance(q_i_src, list):
                q_i = np.array(q_i_src)
            elif isinstance(q_i_src, np.ndarray):
                q_i = q_i_src
            else:
                raise ValueError(
                    "Invalid source for q_i. Must be a list or numpy array.")

            # SECTION: init NRTL model
            # activity model
            activity = ptm.activity(
                components=components, model_name='NRTL')
            # set
            activity_uniquac = activity.uniquac

            # check
            if activity_uniquac is None:
                raise ValueError(
                    "Failed to initialize UNIQUAC activity model. Please check the components and model name.")

            # NOTE: check method
            tau_ij_cal_method = 0
            if dU_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if (a_ij_src is None or
                    b_ij_src is None or
                    c_ij_src is None or
                        d_ij_src is None):
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dU_ij_src is not None:
                # use dU_ij
                if isinstance(dU_ij_src, TableMatrixData):
                    dU_ij = dU_ij_src.mat('dU', components)
                elif isinstance(dU_ij_src, list):
                    dU_ij = np.array(dU_ij_src)
                elif isinstance(dU_ij_src, np.ndarray):
                    dU_ij = dU_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: calculate tau_ij
            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            if tau_ij_cal_method == 1:
                # check
                if isinstance(dU_ij, np.ndarray):
                    tau_ij, _ = activity_uniquac.cal_tau_ij_M1(
                        temperature=temperature,
                        dU_ij=dU_ij)
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be numpy array.")
            elif tau_ij_cal_method == 2:
                # check
                if (isinstance(a_ij, np.ndarray) and
                    isinstance(b_ij, np.ndarray) and
                    isinstance(c_ij, np.ndarray) and
                        isinstance(d_ij, np.ndarray)):
                    # calculate tau_ij
                    tau_ij, _ = activity_uniquac.cal_tau_ij_M2(
                        temperature=temperature,
                        a_ij=a_ij,
                        b_ij=b_ij,
                        c_ij=c_ij,
                        d_ij=d_ij)
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij, b_ij, c_ij, d_ij). Must be numpy array.")
            else:
                raise ValueError(
                    "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # NOTE: nrtl inputs
            inputs_ = {
                'mole_fraction': z_i_comp,
                "tau_ij": tau_ij,
                "r_i": r_i,
                "q_i": q_i
            }

            # NOTE: calculate activity
            res_, _ = activity_uniquac.cal(model_input=inputs_)

            # res format
            # res = {
            #     'property_name': 'activity coefficients',
            #     'components': components,
            #     'mole_fraction': xi,
            #     'value': AcCo_i,
            #     'unit': 1,
            #     'symbol': "AcCo_i",
            #     'message': message,
            # }

            # res
            return res_
        except Exception as e:
            raise Exception(
                f"Failed to calculate UNIQUAC activity: {e}") from e
