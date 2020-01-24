# -*- coding: utf-8 -*-

from .Lgm_Wrap import Lgm_Set_Lgm_B_cdip, Lgm_Set_Lgm_B_edip, Lgm_Set_Lgm_B_OP77, Lgm_Set_Lgm_B_T89
from .Lgm_Wrap import Lgm_Set_Lgm_B_T89c, Lgm_Set_Lgm_B_Dungey, Lgm_Set_Lgm_B_T96, Lgm_Set_Lgm_B_T01S


Bfield_dict = {'Lgm_B_cdip': Lgm_Set_Lgm_B_cdip,
                'Lgm_B_edip': Lgm_Set_Lgm_B_edip,
                'Lgm_B_OP77': Lgm_Set_Lgm_B_OP77,
                'Lgm_B_T89': Lgm_Set_Lgm_B_T89,
                'Lgm_B_T89c': Lgm_Set_Lgm_B_T89c,
                'Lgm_B_T96': Lgm_Set_Lgm_B_T96,
                'Lgm_B_T01S': Lgm_Set_Lgm_B_T01S,
                'Lgm_B_Dungey': Lgm_Set_Lgm_B_Dungey}
