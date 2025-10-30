import typing
from typing import Iterable, Sized, Optional, Union

from .objects.general import ObjectId
from .objects.seqloc import SeqLoc
from .objects.seqid import SeqId
from .objects.seqalign import SeqAlign, SeqAlignSet
from .objects.seq import BioSeq
from .objects.seqset import BioSeqSet
from .objmgr import Scope
from .objtools import DatabaseReader

class SearchQuery:
    def __init__(self, seqloc: SeqLoc, scope: Scope) -> None: ...
    @property
    def seqloc(self) -> SeqLoc: ...
    @property
    def length(self) -> int: ...
    @property
    def scope(self) -> Scope: ...

class SearchQueryVector(Sized):
    def __init__(self, queries: Iterable[SearchQuery] = ()): ...
    def __len__(self) -> int: ...

class SearchResultsSet(typing.Sequence[SearchResult]):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> SearchResult: ...

class SearchResults:
    @property
    def query_id(self) -> SeqId: ...
    @property
    def alignments(self) -> SeqAlignSet: ...

class _BlastOptions(typing.TypedDict):
    evalue: Optional[float] = None,
    gapped: Optional[bool] = None,
    window_size: Optional[int] = None,
    max_target_sequences: Optional[int] = None,
    xdrop_gap: Optional[float] = None,
    culling_limit: Optional[int] = None,
    percent_identity: Optional[float] = None,

class Blast:
    @staticmethod
    def tasks() -> List[str]: ...
    def __init__(
        self,
        **kwargs: _BlastOptions,
    ) -> None: ...
    def __repr__(self) -> str: ...
    def program(self) -> str: ...
    @property
    def window_size(self) -> int: ...
    @window_size.setter
    def window_size(self, window_size: int) -> None: ...
    @property
    def off_diagonal_range(self) -> int: ...
    @property
    def xdrop_gap(self) -> float: ...
    @xdrop_grap.setter
    def xdrop_gap(self, xrdop_gap: float) -> None: ...
    @property
    def evalue(self) -> float: ...
    @evalue.setter
    def evalue(self, evalue: float) -> None: ...
    @property
    def percent_identity(self) -> float: ...
    @percent_identity.setter
    def percent_identity(self, percent_identity: float) -> None: ...
    @property
    def coverage_hsp(self) -> float: ...
    @property
    def gapped(self) -> bool: ...
    @gapped.setter
    def gapped(self, gapped: bool) -> None: ...
    @property
    def culling_limit(self) -> int: ...
    @culling_limit.setter
    def culling_limited(self, culling_limit: int) -> None: ...
    @property
    def database_size(self) -> int: ...
    @property
    def search_space(self) -> int: ...
    @property
    def max_target_sequences(self) -> int: ...
    @max_target_sequences.setter
    def max_target_sequences(self, max_target_sequences: int) -> None: ...
    def run(
        self,
        queries: Union[
            BioSeq, BioSeqSet, SearchQuery, SearchQueryVector, Iterable[SearchQuery]
        ],
        subjects: Union[
            BioSeq,
            BioSeqSet,
            SearchQuery,
            SearchQueryVector,
            Iterable[SearchQuery],
            DatabaseReader,
        ],
        pairwise: bool = False,
    ) -> SearchResultsSet: ...

class NucleotideBlast(Blast):
    pass

class ProteinBlast(Blast):
    pass

class MappingBlast(Blast):
    pass

class BlastP(ProteinBlast):
    def __init__(
        self,
        *,
        word_threshold: Optional[float] = None,
        word_size: Optional[int] = None,
        **kwargs: _BlastOptions,
    ) -> None: ...
    @property
    def word_threshold(self) -> float: ...
    @word_threshold.setter
    def word_threshold(self, word_threshold: float) -> None: ...
    @property
    def word_size(self) -> int: ...
    @word_size.setter
    def word_size(self, word_size: int) -> None: ...

class BlastN(NucleotideBlast):
    def __init__(
        self,
        *,
        dust_filtering: Optional[bool] = None,
        penalty: Optional[int] = None,
        reward: Optional[int] = None,
        **kwargs: _BlastOptions,
    ) -> None: ...
    @property
    def dust_filtering(self) -> bool: ...
    @dust_filtering.setter
    def dust_filtering(self, dust_filtering: bool) -> None: ...
    @property
    def penalty(self) -> int: ...
    @penalty.setter
    def penalty(self, penalty: int) -> None: ...
    @property
    def reward(self) -> int: ...
    @reward.setter
    def reward(self, reward: int) -> None: ...

class BlastX(NucleotideBlast):
    def __init__(
        self,
        *,
        query_genetic_code: int = 1,
        max_intron_length: int = 0,
        **kwargs: _BlastOptions,
    ) -> None: ...
    @property
    def max_intron_length(self) -> int: ...
    @max_intron_length.setter
    def max_intron_length(self, max_intron_length: int) -> None: ...
    @property
    def query_genetic_code(self) -> int: ...
    @query_genetic_code.setter
    def query_genetic_code(self, query_genetic_code: int) -> None: ...

class TBlastN(ProteinBlast):
    def __init__(
        self,
        *,
        database_genetic_code: int = 1,
        max_intron_length: int = 0,
        **kwargs: _BlastOptions,
    ) -> None: ...
    @property
    def max_intron_length(self) -> int: ...
    @max_intron_length.setter
    def max_intron_length(self, max_intron_length: int) -> None: ...
    @property
    def database_genetic_code(self) -> int: ...
    @database_genetic_code.setter
    def database_genetic_code(self, database_genetic_code: int) -> None: ...
