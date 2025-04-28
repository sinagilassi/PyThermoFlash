# ACTIVITY
# ---------
# import libs
import numpy as np
import pyThermoModels as ptm
from typing import List, Dict, Optional
from pyThermoDB import TableMatrixData

# local


class Activity:
    def __init__(self):
        pass

    def __repr__(self):
        return "Activity model class for pyThermoFlash"

    def NRTL(self,
             components: List[str],
             z_i_comp: Dict[str, float],
             temperature: float, **kwargs):
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
            # NOTE: method 1
            # Δg_ij, interaction energy parameter
            dg_ij_src = kwargs.get(
                'interaction-energy-parameter', None) or kwargs.get('dg_ij', None)

            # NOTE: method 2
            # constants a, b, and c
            a_ij_src = kwargs.get('a_ij', None)
            b_ij_src = kwargs.get('b_ij', None)
            c_ij_src = kwargs.get('c_ij', None)

            # NOTE: α_ij, non-randomness parameter
            alpha_ij_src = kwargs.get(
                'non-randomness-parameter', None) or kwargs.get('alpha_ij', None)
            if alpha_ij_src is None:
                raise ValueError(
                    "No valid source provided for non-randomness parameter (α_ij).")

            # SECTION: init NRTL model
            # activity model
            activity = ptm.activity(
                components=components, model_name='NRTL')
            # set
            activity_nrtl = activity.nrtl

            # NOTE: check method
            if dg_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if a_ij_src is None or b_ij_src is None or c_ij_src is None:
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")
                # calculate dg_ij
                dg_ij, _ = activity_nrtl.cal_dg_ij_M1(
                    temperature=temperature,
                    a_ij=a_ij_src, b_ij=b_ij_src, c_ij=c_ij_src)
            elif dg_ij_src is not None:
                # use dg_ij
                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', components)
                elif isinstance(dg_ij_src, List[List[float]]):
                    dg_ij = np.array(dg_ij_src)
                elif isinstance(dg_ij_src, np.ndarray):
                    dg_ij = dg_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: extract data
            # α_ij, non-randomness parameter
            if isinstance(alpha_ij_src, TableMatrixData):
                alpha_ij = alpha_ij_src.mat('alpha', components)
            elif isinstance(alpha_ij_src, List[List[float]]):
                alpha_ij = np.array(alpha_ij_src)
            elif isinstance(alpha_ij_src, np.ndarray):
                alpha_ij = alpha_ij_src
            else:
                raise ValueError(
                    "Invalid source for non-randomness parameter (α_ij). Must be TableMatrixData, list of lists, or numpy array.")

            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            tau_ij, _ = activity_nrtl.cal_tau_ij_M1(
                temperature=temperature, dg_ij=dg_ij)

            # NOTE: nrtl inputs
            inputs_ = {
                'mole_fraction': z_i_comp,
                "tau_ij": tau_ij,
                "alpha_ij": alpha_ij
            }

            # NOTE: calculate activity
            res_, _ = activity_nrtl.cal(model_input=inputs_)

            # res
            return res_
        except Exception as e:
            raise Exception(f"Failed to calculate NRTL activity: {e}") from e
