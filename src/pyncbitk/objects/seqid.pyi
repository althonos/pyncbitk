from typing import Optional

from ..serial import Serial
from .general import ObjectId

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

class OtherId(SeqId):
    def __repr__(self) -> str: ...
    @property
    def id(self) -> TextSeqId: ...

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
