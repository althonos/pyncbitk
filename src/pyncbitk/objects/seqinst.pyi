from typing import Optional

from ..serial import Serial
from .seqdata import SeqData

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

Molecule = Literal["dna", "rna", "protein", "nucleotide", "other"]
Topology = Literal["linear", "circular", "tandem", "other"]
Strandedness = Literal["single", "double", "mixed", "other"]

class SeqInst(Serial):
    def __init__(
        self,
        *,
        topology: Topology = "linear",
        strandedness: Optional[Strandedness] = None,
        molecule: Optional[Molecule] = None,
        length: Optional[int] = None,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def length(self) -> Optional[int]: ...
    @property
    def molecule(self) -> Optional[Molecule]: ...
    @property
    def topology(self) -> Optional[Topology]: ...
    @property
    def strandedness(self) -> Optional[Strandedness]: ...
    @property
    def data(self) -> Optional[SeqData]: ...

class VirtualInst(SeqInst):
    pass

class ContinuousInst(SeqInst):
    def __init__(
        self,
        data: SeqData,
        *,
        topology: Topology = "linear",
        strandedness: Optional[Strandedness] = None,
        molecule: Optional[Molecule] = None,
        length: Optional[int] = None,
    ) -> None: ...
    @property
    def data(self) -> SeqData: ...

class SegmentedInst(SeqInst):
    pass

class ConstructedInst(SeqInst):
    pass

class RefInst(SeqInst):
    pass

class ConsensusInst(SeqInst):
    pass

class MapInst(SeqInst):
    pass

class DeltaInst(SeqInst):
    def to_continuous(self) -> ContinuousInst: ...
