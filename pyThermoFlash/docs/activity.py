# ACTIVITY
# ---------
# import libs
import pyThermoModels as ptm
from typing import List, Dict, Optional
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
            dg_ij_src = kwargs.get('interaction-energy-parameter', None)

            # NOTE: method 2
            # constants a, b, and c
            a_ij_src = kwargs.get('A-parameter', None)
            b_ij_src = kwargs.get('B-parameter', None)
            c_ij_src = kwargs.get('C-parameter', None)

            # NOTE: α_ij, non-randomness parameter
            alpha_ij_src = kwargs.get('non-randomness-parameter', None)

            # SECTION: init NRTL model
            # activity model
            activity = ptm.activity(
                components=components, model_name='NRTL')
            # set
            activity_nrtl = activity.nrtl

            # NOTE: check method
            if dg_ij_src is not None:
                # dg_ij
                for component_i in components:
                    for component_j in components:
                        # check
                        if component_i == component_j:
                            continue
                        else:
                            dg_ij_str = f"dg | {component_i} | {component_j}"
                            dg_ij = dg_ij_src.ijs(dg_ij_str)
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: extract data
            # α_ij
            alpha_ij_str = f"alpha | {components[0]} | {components[1]}"
            alpha_ij = alpha_ij_src.ijs(alpha_ij_str)

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
