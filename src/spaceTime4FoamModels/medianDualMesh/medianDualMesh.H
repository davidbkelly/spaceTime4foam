/*---------------------------------------------------------------------------*\
License
    This file is part of spaceTime4foam.

    spaceTime4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    spaceTime4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spaceTime4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::medianDualMesh

Description
    Median-dual geometry of the space-time triangulation for the
    vertex-centred (edge-based) finite volume method in 1D+t, following
    CLAUDE.md section 7 and Tufillaro, Williams and Nishikawa (2026), Eq. 12,
    18 and Remark 2 with d = 2.

    Topology. The triangulation is the set of faces of the given patch
    (default frontAndBack) that lie on the lower z plane z = zMin (the front
    faces). They must all be triangles and be exactly half of the patch
    faces; the other half must lie on z = zMax (the back faces). The nodes
    are the points of the front faces in the local (patch) numbering, used
    as (x, t) = (x, y) with z set to 0. The edges are the patch edges,
    internal edges first (0 to nInternalEdges() - 1), then boundary edges.

    Edge orientation. OpenFOAM does not order the two points of a patch
    edge, so every edge is stored here as (j, k) with j < k (local node
    labels), and every edge vector n_jk points from j to k, i.e. it is the
    directed area of the dual face seen from j. The vector from k is
    n_kj = -n_jk. Boundary edge i (0 to nBoundaryEdges() - 1) is edge
    nInternalEdges() + i.

    Geometry, with |T| the area of triangle T:
        V_j  = (1/3) sum_(T contains j) |T|
        n_jk = (1/3) sum_(T contains jk) n_j^T + n_B/6 (boundary edges only)
    where n_j^T is the normal of the edge of T opposite j, with length equal
    to that edge's length and pointing away from j, and n_B is the outward
    normal of the boundary edge, with length equal to the edge length.
    Every normal is oriented by an explicit dot-product test (the vector
    from the node, or from the opposite vertex, to the edge midpoint), so no
    sign depends on the face winding written by gmsh or gmshToFoam. The
    signed triangle areas must all have the same sign and none may be near
    zero, otherwise the construction stops with a fatal error.

    Front/back pairing. Every front node is paired with exactly one back
    mesh point with the same (x, y), matched by coordinates to within
    1e-9 times the diagonal of the (x, y) bounding box.

    Boundary edges. Each boundary edge is mapped to the face containing its
    two front mesh points on a patch that is not empty and is not the
    triangulated patch (e.g. tStart); there must be exactly one such face.

    Gradients. Unweighted linear least squares over the edge neighbours:
        grad(u)_j = M_j^-1 sum_k (p_k - p_j) (u_k - u_j),
        M_j = sum_k (p_k - p_j) (p_k - p_j)^T
    with M_j stored as a symmTensor with zz = 1 so that inv() applies. The
    coefficient vectors M_j^-1 (p_k - p_j) and M_k^-1 (p_j - p_k) of each
    edge are precomputed, so gradient() is a single edge loop.

    Mesh size h = sqrt(A_domain/(nTriangles/2)), where A_domain = sum |T|,
    as in spaceTimeErrors (h = 1/N on the N x N triangle meshes).

    Only serial runs are supported.

SourceFiles
    medianDualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef medianDualMesh_H
#define medianDualMesh_H

#include "fvMesh.H"
#include "edgeList.H"
#include "faceList.H"
#include "scalarField.H"
#include "vectorField.H"
#include "symmTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class medianDualMesh Declaration
\*---------------------------------------------------------------------------*/

class medianDualMesh
{
    // Private data

        //- The finite volume mesh
        const fvMesh& mesh_;

        //- Name of the triangulated patch
        const word patchName_;

        //- Index of the triangulated patch
        const label patchID_;

        //- Front triangles in local node labels, in the stored point order
        faceList triangles_;

        //- Node points (x, t, 0)
        pointField points_;

        //- Edges (j, k) with j < k, internal edges first
        edgeList edges_;

        //- Number of internal edges
        label nInternalEdges_;

        //- Edges of each triangle
        labelListList triangleEdges_;

        //- Triangles of each edge
        labelListList edgeTriangles_;

        //- Edges of each node
        labelListList nodeEdges_;

        //- Front mesh point of each node
        labelList frontMeshPoints_;

        //- Back mesh point of each node
        labelList backMeshPoints_;

        //- Mesh face (on a non-empty patch) of each boundary edge
        labelList boundaryEdgeFaces_;

        //- Patch of each boundary edge
        labelList boundaryEdgePatches_;

        //- True for nodes on a boundary edge
        boolList isBoundaryNode_;

        //- Triangle areas |T|
        scalarField triangleAreas_;

        //- Dual areas V_j
        scalarField V_;

        //- Directed-area vector n_jk of each edge, from j to k (j < k)
        vectorField n_;

        //- Outward normal n_B of each boundary edge, length = edge length
        vectorField nB_;

        //- LSQ coefficient M_j^-1 (p_k - p_j) of each edge, for node j
        vectorField lsqCoeffsStart_;

        //- LSQ coefficient M_k^-1 (p_j - p_k) of each edge, for node k
        vectorField lsqCoeffsEnd_;

        //- Domain area A_domain = sum |T|
        scalar domainArea_;

        //- Mesh size h = sqrt(A_domain/(nTriangles/2))
        scalar h_;


    // Private Member Functions

        //- Find the patch index, stopping if it does not exist
        label findPatch() const;

        //- Build the triangulation and the node and edge addressing
        void calcTopology();

        //- Pair every front node with its back mesh point
        void calcFrontBackPairs();

        //- Map every boundary edge to its face and patch
        void calcBoundaryEdgePatches();

        //- Compute the triangle areas, V, n and n_B
        void calcGeometry();

        //- Compute the LSQ coefficients
        void calcLeastSquares();

        //- No copy construct
        medianDualMesh(const medianDualMesh&) = delete;

        //- No copy assignment
        void operator=(const medianDualMesh&) = delete;


public:

    // Static Member Functions

        //- Signed areas of triangles (local labels into points), from the
        //  stored point order: 0.5 ((p1 - p0) x (p2 - p0)) . z
        static tmp<scalarField> signedTriangleAreas
        (
            const faceList& triangles,
            const pointField& points
        );

        //- The node of triangle f that is not on edge e
        static label oppositeNode(const face& f, const edge& e);

        //- Normal of the segment (a, b) with length |b - a|, oriented away
        //  from the point from: (0.5 (a + b) - from) . n > 0
        static vector normalAwayFrom
        (
            const point& a,
            const point& b,
            const point& from
        );


    // Constructors

        //- Construct from a mesh and the name of the triangulated patch
        medianDualMesh
        (
            const fvMesh& mesh,
            const word& patchName = "frontAndBack"
        );


    //- Destructor
    ~medianDualMesh() = default;


    // Member Functions

        // Access

            //- The finite volume mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Name of the triangulated patch
            const word& patchName() const
            {
                return patchName_;
            }

            //- Number of nodes
            label nNodes() const
            {
                return points_.size();
            }

            //- Number of edges
            label nEdges() const
            {
                return edges_.size();
            }

            //- Number of internal edges (edges 0 to nInternalEdges() - 1)
            label nInternalEdges() const
            {
                return nInternalEdges_;
            }

            //- Number of boundary edges
            label nBoundaryEdges() const
            {
                return edges_.size() - nInternalEdges_;
            }

            //- Number of triangles
            label nTriangles() const
            {
                return triangles_.size();
            }

            //- Node points (x, t, 0)
            const pointField& points() const
            {
                return points_;
            }

            //- Triangles in local node labels
            const faceList& triangles() const
            {
                return triangles_;
            }

            //- Edges (j, k) with j < k, internal edges first
            const edgeList& edges() const
            {
                return edges_;
            }

            //- Edges of each triangle
            const labelListList& triangleEdges() const
            {
                return triangleEdges_;
            }

            //- Triangles of each edge
            const labelListList& edgeTriangles() const
            {
                return edgeTriangles_;
            }

            //- Edges of each node
            const labelListList& nodeEdges() const
            {
                return nodeEdges_;
            }

            //- Front mesh point of each node
            const labelList& frontMeshPoints() const
            {
                return frontMeshPoints_;
            }

            //- Back mesh point of each node (same x and t)
            const labelList& backMeshPoints() const
            {
                return backMeshPoints_;
            }

            //- Mesh face of each boundary edge
            const labelList& boundaryEdgeFaces() const
            {
                return boundaryEdgeFaces_;
            }

            //- Patch index of each boundary edge
            const labelList& boundaryEdgePatches() const
            {
                return boundaryEdgePatches_;
            }

            //- Patch name of boundary edge i
            const word& boundaryEdgePatchName(const label i) const
            {
                return mesh_.boundaryMesh()[boundaryEdgePatches_[i]].name();
            }

            //- True for nodes on a boundary edge
            const boolList& isBoundaryNode() const
            {
                return isBoundaryNode_;
            }

            //- Triangle areas |T|
            const scalarField& triangleAreas() const
            {
                return triangleAreas_;
            }

            //- Dual areas V_j
            const scalarField& V() const
            {
                return V_;
            }

            //- Directed-area vector n_jk of each edge, from j to k (j < k)
            const vectorField& n() const
            {
                return n_;
            }

            //- Outward normal n_B of each boundary edge
            const vectorField& nB() const
            {
                return nB_;
            }

            //- Domain area A_domain = sum |T|
            scalar domainArea() const
            {
                return domainArea_;
            }

            //- Mesh size h = sqrt(A_domain/(nTriangles/2))
            scalar h() const
            {
                return h_;
            }


        // Evaluation

            //- LSQ gradient at the nodes of the nodal values u (z = 0)
            tmp<vectorField> gradient(const scalarField& u) const;


        // Output

            //- Print the sizes, domain area and h
            void printSummary() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
