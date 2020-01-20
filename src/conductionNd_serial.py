"""
Copyright 2017 Ben Mather

This file is part of Conduction <https://git.dias.ie/itherc/conduction/>

Conduction is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Conduction is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Conduction.  If not, see <http://www.gnu.org/licenses/>.
"""

try: range = xrange
except: pass

import numpy as np
from scipy import sparse
from scipy.sparse import linalg

class ConductionND(object):
    """
    Implicit ND steady-state heat equation solver over a structured grid
    (Serial version)
    """

    def __init__(self, minCoord, maxCoord, res):

        dim = len(res)
        extent = np.zeros(dim*2)

        index = 0
        sizes = 1

        bbox  = list(range(dim))
        n = np.zeros(dim, dtype=np.int32)

        for i in range(0, dim):
            extent[index]   = minCoord[i]
            extent[index+1] = maxCoord[i]
            index += 2
            sizes *= res[i]
            bbox[i] = (minCoord[i], maxCoord[i])
            n[i] = res[i]

        # Setup matrix sizes
        self.sizes = (sizes, sizes)
        self.dim = dim
        self.extent = extent
        self.bbox = bbox

        nn = sizes
        self.n = n[::-1]
        self.nn = nn

        # stencil size
        self.stencil_width = 2*dim + 1

        # local numbering
        self.nodes = np.arange(0, nn, dtype=np.int32)

        # closure depends on dim
        if dim == 1:
            closure = [(0,-2),(2,0),(1,-1)]
        elif dim == 2:
            closure = [(0,-2), (1,-1), (2,0), (1,-1), (1,-1)]
        elif dim == 3:
            closure = [(0,-2), (1,-1), (1,-1), (2,0), (1,-1), (1,-1), (1,-1)]
        self.closure = self._create_closure_object(closure)


        # interior slices
        self.interior_slice = [None]*dim
        for i in range(0, dim):
            self.interior_slice[i] = slice(1, -1)


        self._initialise_mesh_variables()
        self._initialise_boundary_dictionary()
        self._initialise_COO_vectors()

        # thermal properties
        self.diffusivity  = np.zeros(nn)
        self.heat_sources = np.zeros(nn)
        self.temperature  = np.zeros(nn)

        # right hand side vector
        self.rhs = np.zeros(nn)


    def _initialise_COO_vectors(self):

        nn = self.nn
        n = self.n

        index = np.empty(n + 2, dtype=np.int32)
        index.fill(-1)
        index[self.interior_slice] = self.nodes.reshape(n)
        self.index = index

        self.rows = np.empty((self.stencil_width, nn), dtype=np.int32)
        self.cols = np.empty((self.stencil_width, nn), dtype=np.int32)
        self.vals = np.empty((self.stencil_width, nn))



    def _initialise_mesh_variables(self):

        dim = self.dim
        bbox = self.bbox
        n = self.n[::-1]

        extent = np.zeros(dim*2)

        index = 0
        for bs, be in bbox:
            extent[index]   = bs
            extent[index+1] = be
            index += 2

        self.extent = extent

        # local coordinates
        grid_coords = [None]*dim
        for i in range(0, dim):
            minI, maxI = bbox[i]
            size = n[i]
            grid_coords[i] = np.linspace(minI, maxI, size)


        coord_arrays = np.meshgrid(*grid_coords[::-1], indexing='ij')
        coords = np.empty((self.nn, dim))
        for i in range(0, dim):
            coords[:,i] = coord_arrays[::-1][i].ravel()

        self.grid_coords = grid_coords
        self.coords = coords


    def _initialise_boundary_dictionary(self):

        coords = self.coords
        grid_coords = self.grid_coords
        dim = self.dim

        minCoords = coords.min(axis=0)
        maxCoords = coords.max(axis=0)

        bbox = self.bbox
        n = self.n[::-1]

        # Setup boundary dictionary
        bc = dict()

        wall = [("minX", "maxX"), ("minY", "maxY"), ("minZ", "maxZ")]

        for i in range(0, dim):
            w0, w1 = wall[i]
            c0, c1 = bbox[i]
            m0, m1 = self.coords[:,i] == c0, self.coords[:,i] == c1
            d0 = d1 = (c1 - c0)/(n[i] - 1)

            bc[w0] = {"val": 0.0, "delta": d0, "flux": True, "mask": m0}
            bc[w1] = {"val": 0.0, "delta": d1, "flux": True, "mask": m1}

        self.bc = bc
        self.dirichlet_mask = np.zeros(self.nn, dtype=bool)


    def _create_closure_object(self, closure):

        n = self.n
        obj = [[0] * self.dim for i in range(self.stencil_width)]

        for i in range(0, self.stencil_width):
            # construct slicing object
            for j in range(0, self.dim):
                start, end = closure[i-j]
                obj[i][j] = slice(start, n[j]+end+2)

        return obj


    def update_properties(self, diffusivity, heat_sources):
        """
        Update diffusivity and heat sources
        """
        self.diffusivity = diffusivity
        self.heat_sources = heat_sources


    def boundary_condition(self, wall, val, flux=True):
        """
        Set the boundary conditions on each wall of the domain.
        By default each wall is a Neumann (flux) condition.
        If flux=True, positive val indicates a flux vector towards the centre
        of the domain.

        val can be a vector with the same number of elements as the wall
        """
        wall = str(wall)

        if wall in self.bc:
            self.bc[wall]["val"]  = val
            self.bc[wall]["flux"] = flux
            d = self.bc[wall]

            mask = d['mask']

            if flux:
                self.dirichlet_mask[mask] = False
                self.bc[wall]["val"] /= -d['delta']
            else:
                self.dirichlet_mask[mask] = True

        else:
            raise ValueError("Wall should be one of {}".format(self.bc.keys()))



    def construct_matrix(self, derivative=False):
        """
        Construct the coefficient matrix
        i.e. matrix A in Ax = b

        We vectorise the 7-point stencil for fast matrix insertion.
        An extra border of dummy values around the domain allows for automatic
        Neumann (flux) boundary creation.
        These are stomped on if there are any Dirichlet conditions.

        """

        nodes = self.nodes
        nn = self.nn
        n = self.n
        dim = self.dim

        index = self.index

        rows = self.rows
        cols = self.cols
        vals = self.vals

        dirichlet_mask = self.dirichlet_mask

        u = self.diffusivity.reshape(n)

        k = np.zeros(n + 2)
        k[self.interior_slice] = u

        for i in range(0, self.stencil_width):
            obj = self.closure[i]

            rows[i] = nodes
            cols[i] = index[obj].ravel()

            distance = np.linalg.norm(self.coords[cols[i]] - self.coords, axis=1)
            distance[distance==0] = 1e-12 # protect against dividing by zero
            delta = 1.0/(2.0*distance**2)

            vals[i] = delta*(k[obj] + u).ravel()


        # Dirichlet boundary conditions (duplicates are summed)
        cols[:,dirichlet_mask] = nodes[dirichlet_mask]
        vals[:,dirichlet_mask] = 0.0

        # zero off-grid coordinates
        vals[cols < 0] = 0.0

        # centre point
        vals[-1] = 0.0
        if derivative:
            vals[-1][dirichlet_mask] = 0.
        else:
            vals[-1][dirichlet_mask] = -1.0


        row = rows.ravel()
        col = cols.ravel()
        val = vals.ravel()


        # mask off-grid entries and sum duplicates
        mask = col >= 0
        row = row[mask]
        col = col[mask]
        val = val[mask]

        mat = sparse.coo_matrix((val, (row, col)), shape=self.sizes).tocsr()
        diag = np.ravel(mat.sum(axis=1))
        diag *= -1
        mat.setdiag(diag)

        return mat


    def construct_rhs(self):
        """
        Construct the right-hand-side vector
        i.e. vector b in Ax = b

        Boundary conditions are grabbed from the dictionary and
        summed to the rhs.
        Be careful of duplicate entries on the corners!!
        """
        
        rhs = -1.0*self.heat_sources

        for wall in self.bc:
            val  = self.bc[wall]['val']
            flux = self.bc[wall]['flux']
            mask = self.bc[wall]['mask']
            if flux:
                rhs[mask] += val
            else:
                rhs[mask] = val

        return rhs


    def solve(self, matrix=None, rhs=None):
        """
        Construct the matrix A and vector b in Ax = b
        and solve for x

        GMRES method is default
        """
        if matrix is None:
            matrix = self.construct_matrix()
        if rhs is None:
            rhs = self.construct_rhs()
        res = self.temperature

        T = linalg.spsolve(matrix, rhs)
        self.temperature = T

        return T


    def gradient(self, vector, **kwargs):

        return np.gradient(vector.reshape(self.n), *self.grid_coords[::-1], **kwargs)


    def heatflux(self):

        T = self.temperature
        k = self.diffusivity * -1
        divT = self.gradient(T)
        for i in range(0, self.dim):
            div = k*divT[i].ravel()
            divT[i] = div

        return divT
