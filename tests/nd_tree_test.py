import numpy as np
import pytest
from nd_tree import build_quad_tree, build_octree
from nd_tree.annotations import RealArray, Shape

class TestNDTree():
    @pytest.fixture(params=[10.0,])
    def length(self, request: pytest.FixtureRequest) -> float:
        return request.param

    @pytest.fixture(params=[3,])
    def num_neighbours(self, request: pytest.FixtureRequest) -> int:
        return request.param

    def project_to_box(self, point: RealArray, vmin: RealArray, vmax: RealArray
                       ) -> RealArray:
        return np.clip(point, vmin, vmax)

    def distance(self, box1: RealArray, box2: RealArray) -> RealArray:
        pt1 = self.project_to_box(box1[..., 0, :], box2[..., 0, :], box2[..., 1, :])
        vec = pt1 - self.project_to_box(pt1, box1[..., 0, :], box1[..., 1, :])
        return np.sqrt(np.sum(vec**2, axis=-1))

    def points_on_sphere(self, rng: np.random.Generator, length: float, shape: Shape
                         ) -> RealArray:
        points = rng.normal(size=shape)
        radii = np.sum(points**2, axis=-1)
        return length * points / np.sqrt(radii)[..., None]

    def points_in_cube(self, rng: np.random.Generator, vmin: RealArray, vmax: RealArray,
                       shape: Shape) -> RealArray:
        return rng.uniform(vmin, vmax, shape)

    @pytest.fixture
    def lines(self, rng: np.random.Generator, vmin: RealArray, vmax: RealArray,
              length: float, num_lines: int) -> RealArray:
        first = self.points_in_cube(rng, vmin, vmax, (num_lines, vmin.size))
        second = first + self.points_on_sphere(rng, length, (num_lines, vmin.size))
        return np.sort(np.stack((first, second), axis=-2), axis=-2)

    @pytest.fixture
    def queries(self, rng: np.random.Generator, vmin: RealArray, vmax: RealArray,
                num_queries: int) -> RealArray:
        points = self.points_in_cube(rng, vmin, vmax, (num_queries, vmin.size))
        return np.stack((points, points), axis=-2)

    @pytest.mark.nd_tree
    @pytest.mark.parametrize("vmin,vmax,num_lines,num_queries",
                             [(np.array([0.0, 0.0]), np.array([100.0, 80.0]), 100, 50)])
    def test_quad(self, vmin: RealArray, vmax: RealArray, lines: RealArray,
                  queries: RealArray, num_neighbours: int):
        tree = build_quad_tree(np.reshape(lines, lines.shape[:-2] + (-1,)), vmin, vmax)
        dist, index = tree.find_k_nearest(np.reshape(queries, queries.shape[:-2] + (-1,)),
                                          num_neighbours)
        dmatrix = self.distance(lines[:, None], queries[None])
        index_matrix = np.argsort(dmatrix, axis=0)
        dmatrix = np.take_along_axis(dmatrix, index_matrix, axis=0)

        np.testing.assert_allclose(dist, dmatrix[:num_neighbours].T)
        np.testing.assert_allclose(self.distance(lines[index.T], queries),
                                   dmatrix[:num_neighbours])

    @pytest.mark.nd_tree
    @pytest.mark.parametrize("vmin,vmax,num_lines,num_queries",
                             [(np.array([0.0, 0.0, 0.0]), np.array([100.0, 80.0, 90.0]), 200, 50)])
    def test_oct(self, vmin: RealArray, vmax: RealArray, lines: RealArray,
                 queries: RealArray, num_neighbours: int):
        tree = build_octree(np.reshape(lines, lines.shape[:-2] + (-1,)), vmin, vmax)
        dist, index = tree.find_k_nearest(np.reshape(queries, queries.shape[:-2] + (-1,)),
                                          num_neighbours)
        dmatrix = self.distance(lines[:, None], queries[None])
        index_matrix = np.argsort(dmatrix, axis=0)
        dmatrix = np.take_along_axis(dmatrix, index_matrix, axis=0)

        np.testing.assert_allclose(dist, dmatrix[:num_neighbours].T)
        np.testing.assert_allclose(self.distance(lines[index.T], queries),
                                   dmatrix[:num_neighbours])
