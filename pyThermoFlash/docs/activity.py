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
            # NOTE: check src
            # Δg_ij, interaction energy parameter
            dg_ij_src = kwargs.get('interaction-energy-parameter', None)
            # α_ij, non-randomness parameter
            alpha_ij_src = kwargs.get('non-randomness-parameter', None)

            # SECTION: extract data
            # dg_ij
            dg_ij_str = f"dg | {components[0]} | {components[1]}"
            dg_ij = dg_ij_src.ijs(dg_ij_str)
            # α_ij
            alpha_ij_str = f"alpha | {components[0]} | {components[1]}"
            alpha_ij = alpha_ij_src.ijs(alpha_ij_str)

            # SECTION: init NRTL model
            # activity model
            activity = ptm.activity(
                components=components, model_name='NRTL')
            # set
            activity_nrtl = activity.nrtl

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
