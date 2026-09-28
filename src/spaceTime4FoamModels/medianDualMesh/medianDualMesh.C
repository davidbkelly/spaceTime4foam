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

\*---------------------------------------------------------------------------*/

#include "medianDualMesh.H"
#include "emptyPolyPatch.H"
#include "matchPoints.H"
#include "boundBox.H"
#include "HashSet.H"
#include "DynamicList.H"
#ifdef FOAMEXTEND
    #include "PrimitivePatch.H"
#else
    #include "primitivePatch.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Primitive patch of the front triangles, owning its face list. Only the
    // OpenFOAM.com form has been tested.
#ifdef FOAMEXTEND
    typedef PrimitivePatch<face, List, const pointField&> frontTrianglePatch;
#else
    typedef primitiveFacePatch frontTrianglePatch;
#endif
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::medianDualMesh::signedTriangleAreas
(
    const faceList& triangles,
    const pointField& points
)
{
    tmp<scalarField> tareas(new scalarField(triangles.size()));
    scalarField& areas = tareas.ref();

    forAll(triangles, trii)
    {
        const face& f = triangles[trii];

        if (f.size() != 3)
        {
            FatalErrorInFunction
                << "Face " << trii << " is not a triangle: " << f
                << exit(FatalError);
        }

        const vector d1 = points[f[1]] - points[f[0]];
        const vector d2 = points[f[2]] - points[f[0]];

        areas[trii] = 0.5*(d1.x()*d2.y() - d1.y()*d2.x());
    }

    return tareas;
}


Foam::label Foam::medianDualMesh::oppositeNode(const face& f, const edge& e)
{
    forAll(f, fp)
    {
        if (f[fp] != e.start() && f[fp] != e.end())
        {
            return f[fp];
        }
    }

    FatalErrorInFunction
        << "Triangle " << f << " has no node opposite edge " << e
        << exit(FatalError);

    return -1;
}


Foam::vector Foam::medianDualMesh::normalAwayFrom
(
    const point& a,
    const point& b,
    const point& from
)
{
    const vector d = b - a;

    // Rotated by -90 degrees in (x, t): |normal| = |d|
    vector normal(d.y(), -d.x(), 0);

    // Explicit orientation test, independent of the point order
    if (((0.5*(a + b) - from) & normal) < 0)
    {
        normal = -normal;
    }

    return normal;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::medianDualMesh::findPatch() const
{
    const label patchi = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Cannot find the patch " << patchName_ << ". Patches: "
            << mesh_.boundaryMesh().names() << exit(FatalError);
    }

    return patchi;
}


void Foam::medianDualMesh::calcTopology()
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchID_];
    const pointField& meshPoints = mesh_.points();

    if (pp.empty())
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " has no faces" << exit(FatalError);
    }

    // z range of the patch
    scalar zMin = GREAT;
    scalar zMax = -GREAT;

    forAll(pp.meshPoints(), pointi)
    {
        const scalar z = meshPoints[pp.meshPoints()[pointi]].z();
        zMin = min(zMin, z);
        zMax = max(zMax, z);
    }

    const scalar zTol = 1e-6*(zMax - zMin);

    if (zMax - zMin <= 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " lies on a single z plane: it must"
            << " have front (z = zMin) and back (z = zMax) faces"
            << exit(FatalError);
    }

    // Front faces (all points on z = zMin); all others must be back faces
    DynamicList<label> frontFaces(pp.size()/2);
    label nBack = 0;

    forAll(pp, facei)
    {
        const face& f = pp[facei];
        bool onFront = true;
        bool onBack = true;

        forAll(f, fp)
        {
            const scalar z = meshPoints[f[fp]].z();

            if (mag(z - zMin) > zTol)
            {
                onFront = false;
            }
            if (mag(z - zMax) > zTol)
            {
                onBack = false;
            }
        }

        if (onFront)
        {
            if (f.size() != 3)
            {
                FatalErrorInFunction
                    << "Front face " << pp.start() + facei << " of patch "
                    << patchName_ << " has " << f.size() << " points: all"
                    << " front faces must be triangles" << exit(FatalError);
            }

            frontFaces.append(pp.start() + facei);
        }
        else if (onBack)
        {
            nBack++;
        }
        else
        {
            FatalErrorInFunction
                << "Face " << pp.start() + facei << " of patch " << patchName_
                << " lies neither on z = " << zMin << " nor on z = " << zMax
                << exit(FatalError);
        }
    }

    if (2*frontFaces.size() != pp.size() || nBack != frontFaces.size())
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " has " << pp.size() << " faces, "
            << frontFaces.size() << " on the front and " << nBack
            << " on the back: the front faces must be exactly half of them"
            << exit(FatalError);
    }

    faceList frontFaceList(frontFaces.size());

    forAll(frontFaces, i)
    {
        frontFaceList[i] = mesh_.faces()[frontFaces[i]];
    }

    const frontTrianglePatch front(frontFaceList, meshPoints);

    triangles_ = front.localFaces();
    points_ = front.localPoints();
    points_.replace(vector::Z, scalar(0));
    frontMeshPoints_ = front.meshPoints();
    nInternalEdges_ = front.nInternalEdges();
    triangleEdges_ = front.faceEdges();
    edgeTriangles_ = front.edgeFaces();
    nodeEdges_ = front.pointEdges();

    // Orient every edge from the lower to the higher node label: the edges
    // of a PrimitivePatch are not ordered
    const edgeList& frontEdges = front.edges();
    edges_.setSize(frontEdges.size());

    forAll(frontEdges, edgei)
    {
        const edge& e = frontEdges[edgei];

        edges_[edgei] = edge(min(e.start(), e.end()), max(e.start(), e.end()));
    }

    // Manifold: two triangles per internal edge, one per boundary edge
    forAll(edgeTriangles_, edgei)
    {
        const label nExpected = (edgei < nInternalEdges_ ? 2 : 1);

        if (edgeTriangles_[edgei].size() != nExpected)
        {
            FatalErrorInFunction
                << "Edge " << edgei << " " << edges_[edgei] << " has "
                << edgeTriangles_[edgei].size() << " triangles, expected "
                << nExpected << exit(FatalError);
        }
    }
}


void Foam::medianDualMesh::calcFrontBackPairs()
{
    const polyPatch& pp = mesh_.boundaryMesh()[patchID_];
    const pointField& meshPoints = mesh_.points();

    // Back mesh points: the patch points that are not front points
    labelHashSet frontSet(frontMeshPoints_);
    labelHashSet backSet(pp.meshPoints().size());

    forAll(pp.meshPoints(), pointi)
    {
        const label meshPointi = pp.meshPoints()[pointi];

        if (!frontSet.found(meshPointi))
        {
            backSet.insert(meshPointi);
        }
    }

    const labelList backPoints(backSet.sortedToc());

    if (backPoints.size() != nNodes())
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " has " << nNodes()
            << " front points and " << backPoints.size() << " back points"
            << exit(FatalError);
    }

    // Match by (x, y) only
    pointField backXY(backPoints.size());

    forAll(backPoints, i)
    {
        backXY[i] = meshPoints[backPoints[i]];
        backXY[i].z() = 0;
    }

    const scalar pairTol = 1e-9*mag(boundBox(points_).span());

    labelList frontToBack;

    const bool allMatched = matchPoints
    (
        points_,
        backXY,
        scalarField(nNodes(), pairTol),
        false,
        frontToBack
    );

    if (!allMatched)
    {
        label nUnpaired = 0;

        forAll(frontToBack, nodei)
        {
            if (frontToBack[nodei] < 0)
            {
                nUnpaired++;
            }
        }

        FatalErrorInFunction
            << nUnpaired << " front points of patch " << patchName_
            << " have no back point with the same (x, y) to within "
            << pairTol << exit(FatalError);
    }

    // Every back point must be used exactly once
    boolList used(backPoints.size(), false);
    backMeshPoints_.setSize(nNodes());

    forAll(frontToBack, nodei)
    {
        const label backi = frontToBack[nodei];

        if (used[backi])
        {
            FatalErrorInFunction
                << "Back point " << backPoints[backi] << " is paired with"
                << " more than one front point" << exit(FatalError);
        }

        used[backi] = true;
        backMeshPoints_[nodei] = backPoints[backi];
    }
}


void Foam::medianDualMesh::calcBoundaryEdgePatches()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const labelListList& pointFaces = mesh_.pointFaces();
    const faceList& faces = mesh_.faces();

    boundaryEdgeFaces_.setSize(nBoundaryEdges());
    boundaryEdgePatches_.setSize(nBoundaryEdges());
    isBoundaryNode_.setSize(nNodes());
    isBoundaryNode_ = false;

    for (label i = 0; i < nBoundaryEdges(); i++)
    {
        const edge& e = edges_[nInternalEdges_ + i];
        const label meshPointA = frontMeshPoints_[e.start()];
        const label meshPointB = frontMeshPoints_[e.end()];

        // The side face through the edge on a non-empty patch
        label nFound = 0;
        label foundFace = -1;

        forAll(pointFaces[meshPointA], pfi)
        {
            const label facei = pointFaces[meshPointA][pfi];

            if (mesh_.isInternalFace(facei))
            {
                continue;
            }

            const label patchi = bMesh.whichPatch(facei);

            if (patchi == patchID_ || isA<emptyPolyPatch>(bMesh[patchi]))
            {
                continue;
            }

            if (faces[facei].found(meshPointB))
            {
                nFound++;
                foundFace = facei;
            }
        }

        if (nFound != 1)
        {
            FatalErrorInFunction
                << "Boundary edge " << nInternalEdges_ + i << " (mesh points "
                << meshPointA << " and " << meshPointB << ") is on "
                << nFound << " faces of non-empty patches, expected 1"
                << exit(FatalError);
        }

        boundaryEdgeFaces_[i] = foundFace;
        boundaryEdgePatches_[i] = bMesh.whichPatch(foundFace);
        isBoundaryNode_[e.start()] = true;
        isBoundaryNode_[e.end()] = true;
    }
}


void Foam::medianDualMesh::calcGeometry()
{
    // Signed areas: one sign, none near zero
    const tmp<scalarField> tsignedAreas =
        signedTriangleAreas(triangles_, points_);
    const scalarField& signedAreas = tsignedAreas();

    const scalar meanArea = sum(mag(signedAreas))/nTriangles();
    const scalar areaTol = 1e-8*meanArea;
    label nPositive = 0;
    label nNegative = 0;

    forAll(signedAreas, trii)
    {
        if (mag(signedAreas[trii]) <= areaTol)
        {
            FatalErrorInFunction
                << "Triangle " << trii << " " << triangles_[trii]
                << " has area " << signedAreas[trii] << " (mean triangle area "
                << meanArea << ")" << exit(FatalError);
        }

        if (signedAreas[trii] > 0)
        {
            nPositive++;
        }
        else
        {
            nNegative++;
        }
    }

    if (nPositive > 0 && nNegative > 0)
    {
        FatalErrorInFunction
            << "The triangles of patch " << patchName_ << " are not"
            << " consistently oriented: " << nPositive << " have a positive"
            << " and " << nNegative << " a negative signed area (a folded"
            << " or inverted mesh)" << exit(FatalError);
    }

    triangleAreas_ = mag(signedAreas);
    domainArea_ = sum(triangleAreas_);
    h_ = Foam::sqrt(domainArea_/(0.5*nTriangles()));

    // Dual areas
    V_.setSize(nNodes());
    V_ = 0;

    forAll(triangles_, trii)
    {
        const face& f = triangles_[trii];

        forAll(f, fp)
        {
            V_[f[fp]] += triangleAreas_[trii]/3.0;
        }
    }

    // Interior contribution (1/3) n_j^T for every edge (j, k) of every
    // triangle, seen from the lower label j
    n_.setSize(nEdges());
    n_ = Zero;

    forAll(triangles_, trii)
    {
        const face& f = triangles_[trii];
        const labelList& fEdges = triangleEdges_[trii];

        forAll(fEdges, fei)
        {
            const label edgei = fEdges[fei];
            const edge& e = edges_[edgei];
            const label j = e.start();
            const label k = e.end();
            const label l = oppositeNode(f, e);

            // n_j^T: normal of the edge (k, l) opposite j, away from j
            n_[edgei] +=
                normalAwayFrom(points_[k], points_[l], points_[j])/3.0;
        }
    }

    // Boundary correction n_B/6, with n_B pointing away from the opposite
    // node of the edge's triangle, i.e. out of the domain
    nB_.setSize(nBoundaryEdges());

    for (label i = 0; i < nBoundaryEdges(); i++)
    {
        const label edgei = nInternalEdges_ + i;
        const edge& e = edges_[edgei];
        const face& f = triangles_[edgeTriangles_[edgei][0]];
        const label l = oppositeNode(f, e);

        nB_[i] =
            normalAwayFrom(points_[e.start()], points_[e.end()], points_[l]);
        n_[edgei] += nB_[i]/6.0;
    }
}


void Foam::medianDualMesh::calcLeastSquares()
{
    // M_j = sum_k (p_k - p_j) (p_k - p_j)^T
    symmTensorField M(nNodes(), Zero);

    forAll(edges_, edgei)
    {
        const edge& e = edges_[edgei];
        const vector d = points_[e.end()] - points_[e.start()];

        M[e.start()] += sqr(d);
        M[e.end()] += sqr(d);
    }

    symmTensorField invM(nNodes());

    forAll(M, nodei)
    {
        symmTensor& Mj = M[nodei];

        // The (x, t) block is 2 x 2: zz = 1 makes the tensor invertible
        Mj.zz() = 1;

        const scalar det2D = Mj.xx()*Mj.yy() - sqr(Mj.xy());
        const scalar scale = sqr(0.5*(Mj.xx() + Mj.yy()));

        if (det2D <= 1e-12*scale)
        {
            FatalErrorInFunction
                << "The least-squares matrix of node " << nodei << " at "
                << points_[nodei] << " is singular: its edge neighbours are"
                << " collinear" << exit(FatalError);
        }

        invM[nodei] = inv(Mj);
    }

    lsqCoeffsStart_.setSize(nEdges());
    lsqCoeffsEnd_.setSize(nEdges());

    forAll(edges_, edgei)
    {
        const edge& e = edges_[edgei];
        const vector d = points_[e.end()] - points_[e.start()];

        lsqCoeffsStart_[edgei] = invM[e.start()] & d;
        lsqCoeffsEnd_[edgei] = invM[e.end()] & (-d);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::medianDualMesh::medianDualMesh
(
    const fvMesh& mesh,
    const word& patchName
)
:
    mesh_(mesh),
    patchName_(patchName),
    patchID_(findPatch()),
    triangles_(),
    points_(),
    edges_(),
    nInternalEdges_(0),
    triangleEdges_(),
    edgeTriangles_(),
    nodeEdges_(),
    frontMeshPoints_(),
    backMeshPoints_(),
    boundaryEdgeFaces_(),
    boundaryEdgePatches_(),
    isBoundaryNode_(),
    triangleAreas_(),
    V_(),
    n_(),
    nB_(),
    lsqCoeffsStart_(),
    lsqCoeffsEnd_(),
    domainArea_(0),
    h_(0)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "medianDualMesh supports serial runs only" << exit(FatalError);
    }

    calcTopology();
    calcFrontBackPairs();
    calcBoundaryEdgePatches();
    calcGeometry();
    calcLeastSquares();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::medianDualMesh::gradient
(
    const scalarField& u
) const
{
    if (u.size() != nNodes())
    {
        FatalErrorInFunction
            << "Field size " << u.size() << " is not the number of nodes "
            << nNodes() << exit(FatalError);
    }

    tmp<vectorField> tgrad(new vectorField(nNodes(), Zero));
    vectorField& grad = tgrad.ref();

    forAll(edges_, edgei)
    {
        const label j = edges_[edgei].start();
        const label k = edges_[edgei].end();
        const scalar du = u[k] - u[j];

        grad[j] += lsqCoeffsStart_[edgei]*du;
        grad[k] -= lsqCoeffsEnd_[edgei]*du;
    }

    return tgrad;
}


void Foam::medianDualMesh::printSummary() const
{
    label nBoundaryNodes = 0;

    forAll(isBoundaryNode_, nodei)
    {
        if (isBoundaryNode_[nodei])
        {
            nBoundaryNodes++;
        }
    }

    Info<< "Median-dual mesh of patch " << patchName_ << " (z = zMin faces):"
        << nl
        << "    nodes:      " << nNodes() << " (" << nNodes() - nBoundaryNodes
        << " interior, " << nBoundaryNodes << " boundary)" << nl
        << "    edges:      " << nEdges() << " (" << nInternalEdges_
        << " internal, " << nBoundaryEdges() << " boundary)" << nl
        << "    triangles:  " << nTriangles() << nl
        << "    domain area " << domainArea_ << nl
        << "    h = sqrt(A_domain/(nTriangles/2)) = " << h_ << endl;
}


// ************************************************************************* //
