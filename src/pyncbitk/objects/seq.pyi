from typing import Iterable, List

from ..serial import Serial
from .seqid import SeqId
from .seqinst import SeqInst
from .seqdesc import SeqDesc, SeqDescSet

class BioSeq(Serial):
    def __init__(
        self, 
        instance: SeqInst, 
        id: SeqId, 
        *ids: SeqId,
        descriptions: Iterable[SeqDesc] = (),
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def id(self) -> SeqId: ...
    @property
    def ids(self) -> List[SeqId]: ...
    @property
    def instance(self) -> SeqInst: ...
    @property
    def descriptions(self) -> SeqDescSet: ...