from typing import Dict, List, Tuple, Union, overload
from ..annotations import CPPRealSequence, CPPIntSequence, IntArray, RealArray

# 2D case : QuadTree

class QuadTreeFloat:
    box     : List[float]
    size    : int

    def __init__(self, database: RealArray, vmin: CPPRealSequence, vmax: CPPRealSequence):
        ...

    def find_k_nearest(self, query: RealArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: RealArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

class QuadTreeDouble(QuadTreeFloat):
    ...

class QuadTreeInt:
    box     : List[int]
    size    : int

    def __init__(self, database: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence):
        ...

    def find_k_nearest(self, query: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

@overload
def build_quad_tree(database: RealArray, vmin: CPPRealSequence, vmax: CPPRealSequence
                    ) -> Union[QuadTreeDouble, QuadTreeFloat]:
    ...

@overload
def build_quad_tree(database: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence) -> QuadTreeInt:
    ...

def build_quad_tree(database: Union[RealArray, IntArray], vmin: Union[CPPRealSequence, CPPIntSequence],
                    vmax: Union[CPPRealSequence, CPPIntSequence]
                   ) -> Union[QuadTreeDouble, QuadTreeFloat, QuadTreeInt]:
    ...

class QuadStackFloat:
    trees : Dict[int, QuadTreeFloat]

    def __init__(self, database: RealArray, indices: IntArray, vmin: CPPRealSequence,
                 vmax: CPPRealSequence):
        ...

    def find_k_nearest(self, query: RealArray, indices: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: RealArray, indices: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

class QuadStackDouble(QuadStackFloat):
    trees : Dict[int, QuadTreeDouble]

class QuadStackInt:
    trees : Dict[int, QuadTreeInt]

    def __init__(self, database: IntArray, indices: IntArray, vmin: CPPIntSequence,
                 vmax: CPPIntSequence):
        ...

    def find_k_nearest(self, query: IntArray, indices: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: IntArray, indices: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

@overload
def build_quad_stack(database: RealArray, indices: IntArray, vmin: CPPRealSequence, vmax: CPPRealSequence
                     ) -> Union[QuadStackDouble, QuadStackFloat]:
    ...

@overload
def build_quad_stack(database: IntArray, indices: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence
                     ) -> QuadStackInt:
    ...

def build_quad_stack(database: Union[RealArray, IntArray], indices: IntArray,
                     vmin: Union[CPPRealSequence, CPPIntSequence],
                     vmax: Union[CPPRealSequence, CPPIntSequence]
                    ) -> Union[QuadStackDouble, QuadStackFloat, QuadStackInt]:
    ...

# 3D case : Octree

class OctreeFloat:
    box     : List[float]
    size    : int

    def __init__(self, database: RealArray, vmin: CPPRealSequence, vmax: CPPRealSequence):
        ...

    def find_k_nearest(self, query: RealArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: RealArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

class OctreeDouble(OctreeFloat):
    ...

class OctreeInt:
    box     : List[int]
    size    : int

    def __init__(self, database: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence):
        ...

    def find_k_nearest(self, query: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

@overload
def build_octree(database: RealArray, vmin: CPPRealSequence, vmax: CPPRealSequence
                 ) -> Union[OctreeDouble, OctreeFloat]:
    ...

@overload
def build_octree(database: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence) -> OctreeInt:
    ...

def build_octree(database: Union[RealArray, IntArray], vmin: Union[CPPRealSequence, CPPIntSequence],
                 vmax: Union[CPPRealSequence, CPPIntSequence]
                 ) -> Union[OctreeDouble, OctreeFloat, OctreeInt]:
    ...

class OctStackFloat:
    trees : Dict[int, OctreeFloat]

    def __init__(self, database: RealArray, indices: IntArray, vmin: CPPRealSequence,
                 vmax: CPPRealSequence):
        ...

    def find_k_nearest(self, query: RealArray, indices: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: RealArray, indices: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

class OctStackDouble(OctStackFloat):
    trees : Dict[int, OctreeDouble]

class OctStackInt:
    trees : Dict[int, OctreeInt]

    def __init__(self, database: IntArray, indices: IntArray, vmin: CPPIntSequence,
                 vmax: CPPIntSequence):
        ...

    def find_k_nearest(self, query: IntArray, indices: IntArray, k: int=1, num_threads: int=1
                       ) -> Tuple[RealArray, IntArray]:
        ...

    def find_range(self, query: IntArray, indices: IntArray, range: float, num_threads: int=1
                   ) -> Tuple[List[List[float]], List[List[int]]]:
        ...

@overload
def build_oct_stack(database: RealArray, indices: IntArray, vmin: CPPRealSequence, vmax: CPPRealSequence
                   ) -> Union[OctStackDouble, OctStackFloat]:
    ...

@overload
def build_oct_stack(database: IntArray, indices: IntArray, vmin: CPPIntSequence, vmax: CPPIntSequence
                    ) -> OctStackInt:
    ...

def build_oct_stack(database: Union[RealArray, IntArray], indices: IntArray,
                   vmin: Union[CPPRealSequence, CPPIntSequence],
                   vmax: Union[CPPRealSequence, CPPIntSequence]
                   ) -> Union[OctStackDouble, OctStackFloat, OctStackInt]:
    ...
