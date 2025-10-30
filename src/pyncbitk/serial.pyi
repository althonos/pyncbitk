try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

Format = Literal["asntext", "asnbinary", "xml", "json"]

class Serial:
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def dumps(
        self, format: Format = "asntext", indent: bool = True, eol: bool = True
    ) -> bytes: ...
