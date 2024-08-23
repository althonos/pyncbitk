from typing import Optional

from ..serial import Serial
from .general import ObjectId

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

# --- SeqId --------------------------------------------------------------------

class SeqId(Serial):
    def __eq__(self, other: object) -> bool: ...
    @staticmethod
    def parse(text: str) -> SeqId: ...

class LocalId(SeqId):
    def __init__(self, object_id: ObjectId) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def object_id(self) -> ObjectId: ...

class RefSeqId(SeqId):
    pass

class GenBankId(SeqId):
    def __repr__(self) -> str: ...
    @property
    def id(self) -> TextSeqId: ...

class ProteinDataBankId(SeqId):
    pass

class GeneralId(SeqId):
    pass

# --- TextSeqId ----------------------------------------------------------------

class TextSeqId(Serial):
    def __init__(
        self,
        accession: str,
        *,
        name: Optional[str] = None,
        version: int = 0,
        release: Optional[str] = None,
        allow_dot_version: bool = True,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def accession(self) -> Optional[str]: ...
    @property
    def name(self) -> Optional[str]: ...
    @property
    def version(self) -> Optional[int]: ...
    @property
    def release(self) -> Optional[str]: ...
