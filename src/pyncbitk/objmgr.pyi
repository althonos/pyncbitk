from typing import Optional, ContextManager

from .objects.seq import BioSeq
from .objects.seqid import SeqId

class BioSeqHandle:   
    def __call__(self) -> BioSeq: ...
    @property
    def scope(self) -> Scope: ...
    @property
    def id(self) -> Optional[SeqId]: ...


class ObjectManager:
    def __init__(self) -> None: ...
    def scope(self) -> Scope: ...


class Scope(ContextManager["Scope"]):
    def __init__(self, manager: ObjectManager) -> None: ...
    def __enter__(self) -> Scope: ...
    def __exit__(self, exc_ty, exc_val, traceback) -> bool: ...
    def __contains__(self, key: object) -> bool: ...
    def __getitem__(self, key: SeqId) -> BioSeqHandle: ...
    def close(self) -> None: ...
    def add_bioseq(self, seq: BioSeq) -> BioSeqHandle: ...
