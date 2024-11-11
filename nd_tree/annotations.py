from typing import Any, Sequence, Tuple, Union
import numpy as np
import numpy.typing as npt

Shape = Tuple[int, ...]

Array = npt.NDArray[Any]
ArrayLike = npt.ArrayLike
BoolArray = npt.NDArray[np.bool_]
IntArray = npt.NDArray[np.integer[Any]]
RealArray = npt.NDArray[np.floating[Any]]
ComplexArray = npt.NDArray[Union[np.floating[Any], np.complexfloating[Any, Any]]]

IntSequence = Union[int, np.integer[Any], Sequence[int], IntArray]
CPPIntSequence = Union[Sequence[int], IntArray]
RealSequence = Union[float, np.floating[Any], Sequence[float], RealArray]
CPPRealSequence = Union[Sequence[float], RealArray]