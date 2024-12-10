from typing import Optional

from ..serial import Serial
from .general import ObjectId
from .seqid import SeqId

# --- SeqLoc -------------------------------------------------------------------

class SeqLoc(Serial):
    pass

class NullLoc(SeqLoc):
    pass

class EmptySeqLoc(SeqLoc):
    pass

class WholeSeqLoc(SeqLoc):
    def __init__(self, sequence_id: SeqId) -> None: ...
    @property
    def sequence_id(self) -> SeqId: ...

class SeqIntervalLoc(SeqLoc):
    @property
    def sequence_id(self) -> SeqId: ...
    @property
    def start(self) -> int: ...
    @property
    def stop(self) -> int: ...

class PackedSeqLoc(SeqLoc):
    pass

class PointLoc(SeqLoc):
    pass

class PackedPointsLoc(SeqLoc):
    pass

class MixLoc(SeqLoc):
    pass

class EquivalentLoc(SeqLoc):
    pass

class BondLoc(SeqLoc):
    pass

class FeatureLoc(SeqLoc):
    pass
