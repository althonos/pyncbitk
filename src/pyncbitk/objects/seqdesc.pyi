import typing
from typing import SupportsIndex, Tuple, Type, Union, Sized, Iterable

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

from ..serial import Serial

# --- SeqDesc ------------------------------------------------------------------

class SeqDesc(Serial):
    pass

class NameDesc(SeqDesc):
    def __init__(self, data: str) -> None: ...
    def __str__(self) -> str: ...

class TitleDesc(SeqDesc):
    def __init__(self, data: str) -> None: ...
    def __str__(self) -> str: ...

class RegionDesc(SeqDesc):
    def __init__(self, data: str) -> None: ...
    def __str__(self) -> str: ...

# --- SeqDescSet ---------------------------------------------------------------

class SeqDescSet(Serial, Sized, Iterable[SeqDesc]):
    def __init__(self, items: Iterable[SeqDesc] = ()): ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[SeqDesc]: ...