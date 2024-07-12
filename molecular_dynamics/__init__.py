
import matplotlib.pyplot as plt

from .molecular_dynamics_pbc import MolecularDynamicsPBC
from .molecular_dynamics_bhw import MolecularDynamicsBHW


plt.rcParams['animation.embed_limit'] = 2**128 #allow for large memory use of animations

__all__ = ["MolecularDynamicsPBC", 
           "MolecularDynamicsBHW"
           ]