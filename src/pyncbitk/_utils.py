import typing
import functools

_T = typing.TypeVar("_T")


class peekable(typing.Iterator[_T], typing.Generic[_T]):
    """Turn an iterable into a peekable iterable.

    Example:
        >>> it = peekable(range(3))
        >>> it.peek()
        0
        >>> next(it)
        0
        >>> next(it)
        1
        >>> it.peek()
        2

    """

    _sentinel = object()

    def __init__(self, iterable: typing.Iterable[_T]):
        self.it = iter(iterable)
        self.peeked = self._sentinel

    def __iter__(self) -> "peekable[_T]":
        return self

    def __next__(self) -> _T:
        if self.peeked is not self._sentinel:
            self.peeked, item = self._sentinel, self.peeked
        else:
            item = next(self.it)
        return item  # type: ignore

    def peek(self) -> _T:
        if self.peeked is self._sentinel:
            self.peeked = next(self.it)
        return self.peeked  # type: ignore


def is_iterable(obj: object) -> bool:
    try:
        iter(obj)  # type: ignore
    except TypeError:
        return False
    else:
        return True
