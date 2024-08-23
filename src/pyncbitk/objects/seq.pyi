from ..serial import Serial
from .seqloc import SeqId
from .seqinst import SeqInst

class BioSeq(Serial):
    def __init__(self, instance: SeqInst, id: SeqId, *ids: SeqId) -> None: ...
