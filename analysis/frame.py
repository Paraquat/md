import scipy as sp
from pdb import set_trace


class frame:
  """Stores a frame from a MD trajectory"""

  def __init__(self, nat):
    self.step = 0
    self.species = sp.zeros(nat)
    self.r = sp.zeros((nat, 3), dtype='int')
    self.v = sp.zeros((nat, 3))
    self.f = sp.zeros((nat, 3))
    self.lat = sp.zeros((3, 3))
    self.stress = sp.zeros((3, 3))
    self.ke = 0.
    self.pe = 0.
    self.E = 0.
    self.T = 0.
    self.P = 0.
