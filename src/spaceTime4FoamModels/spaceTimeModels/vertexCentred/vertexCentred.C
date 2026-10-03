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

#include "vertexCentred.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace spaceTimeModels
{
    defineTypeNameAndDebug(vertexCentred, 0);
    addToRunTimeSelectionTable(spaceTimeModel, vertexCentred, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::spaceTimeModels::vertexCentred::readOption
(
    const word& key,
    const wordList& options
) const
{
    const word value(coeffs_.get<word>(key));

    if (!options.found(value))
    {
        FatalIOErrorInFunction(coeffs_)
            << "Unknown " << key << " " << value << ". Valid options: "
            << options << exit(FatalIOError);
    }

    return value;
}


void Foam::spaceTimeModels::vertexCentred::calcUnknowns()
{
    const boolList& isBoundaryNode = dual_.isBoundaryNode();

    isUnknown_.setSize(dual_.nNodes());
    nUnknowns_ = 0;

    forAll(isUnknown_, nodei)
    {
        isUnknown_[nodei] = !(strongBoundary_ && isBoundaryNode[nodei]);

        if (isUnknown_[nodei])
        {
            nUnknowns_++;
        }
    }

    if (nUnknowns_ == 0)
    {
        FatalErrorInFunction
            << "The mesh has no unknown nodes" << exit(FatalError);
    }
}


void Foam::spaceTimeModels::vertexCentred::calcBoundaryData()
{
    uExact_ = analytical().value(dual_.points());

    // Inflow where A . n_B < 0, from the sign alone (not the patch name)
    const vectorField& nB = dual_.nB();

    isInflowEdge_.setSize(dual_.nBoundaryEdges());

    forAll(nB, i)
    {
        isInflowEdge_[i] = ((A() & nB[i]) < 0);
    }
}


void Foam::spaceTimeModels::vertexCentred::calcSource()
{
    sourceV_.setSize(dual_.nNodes());
    sourceV_ = 0;

    if (analytical().hasSource())
    {
        sourceV_ = analytical().source(dual_.points())*dual_.V();
    }
}


void Foam::spaceTimeModels::vertexCentred::calcPseudoTimeStep()
{
    const edgeList& edges = dual_.edges();
    const vectorField& n = dual_.n();
    const vectorField& nB = dual_.nB();
    const label nInternalEdges = dual_.nInternalEdges();

    // D_j = sum_k 0.5 |n_jk . A| (Eq. 103)
    scalarField D(dual_.nNodes(), 0);

    forAll(edges, edgei)
    {
        const scalar w = 0.5*mag(n[edgei] & A());

        D[edges[edgei].start()] += w;
        D[edges[edgei].end()] += w;
    }

    // Weak closure (deviation from the paper, CLAUDE.md section 7): the
    // boundary flux of area |n_B|/2 at each node of a boundary edge adds
    // 0.5 |n_B . A|/2. This changes the convergence path only.
    if (!strongBoundary_)
    {
        forAll(nB, i)
        {
            const edge& e = edges[nInternalEdges + i];
            const scalar w = 0.5*mag(nB[i] & A())/2.0;

            D[e.start()] += w;
            D[e.end()] += w;
        }
    }

    dTau_.setSize(dual_.nNodes());
    dTau_ = 0;

    forAll(D, nodei)
    {
        if (!isUnknown_[nodei])
        {
            continue;
        }

        if (D[nodei] <= SMALL*dual_.h())
        {
            FatalErrorInFunction
                << "Node " << nodei << " at " << dual_.points()[nodei]
                << " has no flux through its dual cell (D = " << D[nodei]
                << "): the pseudo-time step is undefined" << exit(FatalError);
        }

        dTau_[nodei] = CFL_*dual_.V()[nodei]/D[nodei];
    }
}


void Foam::spaceTimeModels::vertexCentred::initialise()
{
    u_.setSize(dual_.nNodes());

    if (initialisation_ == "analytical")
    {
        u_ = uExact_;
    }
    else
    {
        u_ = 0;

        // Strong boundary treatment: the boundary nodes hold the exact data
        forAll(u_, nodei)
        {
            if (!isUnknown_[nodei])
            {
                u_[nodei] = uExact_[nodei];
            }
        }
    }
}


void Foam::spaceTimeModels::vertexCentred::addSource(scalarField& res) const
{
    res -= sourceV_;
}


void Foam::spaceTimeModels::vertexCentred::addWeakBoundaryFluxes
(
    const scalarField& u,
    scalarField& res
) const
{
    const edgeList& edges = dual_.edges();
    const vectorField& nB = dual_.nB();
    const label nInternalEdges = dual_.nInternalEdges();

    forAll(nB, i)
    {
        const edge& e = edges[nInternalEdges + i];
        const scalar magNB = mag(nB[i]);
        const scalar lambda = (nB[i]/magNB) & A();

        // Each node of the edge gets half of the edge
        const scalar area = 0.5*magNB;

        forAll(e, ei)
        {
            const label m = e[ei];

            // Boundary state: analytical data on inflow, the node's own
            // value on outflow
            scalar ub = u[m];

            if (isInflowEdge_[i])
            {
                ub = uExact_[m];
            }

            const scalar phi =
                0.5*lambda*(u[m] + ub) - 0.5*mag(lambda)*(ub - u[m]);

            res[m] += phi*area;
        }
    }
}


void Foam::spaceTimeModels::vertexCentred::updatePointField()
{
    scalarField& uP = uPoint_.primitiveFieldRef();
    const labelList& front = dual_.frontMeshPoints();
    const labelList& back = dual_.backMeshPoints();

    forAll(u_, nodei)
    {
        uP[front[nodei]] = u_[nodei];
        uP[back[nodei]] = u_[nodei];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeModels::vertexCentred::vertexCentred(Time& runTime)
:
    spaceTimeModel(typeName, runTime),
    coeffs_(spaceTimeProperties().subDict(typeName + "Coeffs")),
    linearReconstruction_
    (
        readOption("reconstruction", wordList({"none", "linear"}))
     == "linear"
    ),
    strongBoundary_
    (
        readOption("boundaryTreatment", wordList({"weak", "strong"}))
     == "strong"
    ),
    CFL_(coeffs_.get<scalar>("CFL")),
    initialisation_
    (
        readOption("initialisation", wordList({"zero", "analytical"}))
    ),
    dual_(mesh()),
    isUnknown_(),
    nUnknowns_(0),
    uExact_(),
    sourceV_(),
    isInflowEdge_(),
    dTau_(),
    u_(),
    res_(),
    uPoint_
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh()),
        dimensionedScalar(dimless, Zero)
    ),
    residual_(GREAT)
{
    if (CFL_ <= 0)
    {
        FatalIOErrorInFunction(coeffs_)
            << "CFL must be positive, not " << CFL_ << exit(FatalIOError);
    }

    // Every mesh point must carry a nodal value
    if (2*dual_.nNodes() != mesh().nPoints())
    {
        FatalErrorInFunction
            << "The mesh has " << mesh().nPoints() << " points, but the"
            << " front patch has " << dual_.nNodes() << " nodes: the mesh"
            << " must be one cell thick (front and back points only)"
            << exit(FatalError);
    }

    dual_.printSummary();

    calcUnknowns();
    calcBoundaryData();
    calcSource();
    calcPseudoTimeStep();
    initialise();

    res_.setSize(dual_.nNodes());
    calcResidual(u_, res_);
    residual_ = rmsResidual(res_);

    updatePointField();

    label nInflowEdges = 0;

    forAll(isInflowEdge_, i)
    {
        if (isInflowEdge_[i])
        {
            nInflowEdges++;
        }
    }

    Info<< nl << "vertexCentred settings:" << nl
        << "    reconstruction     "
        << (linearReconstruction_ ? "linear (VC-2)" : "none (VC-1)") << nl
        << "    boundaryTreatment  "
        << (strongBoundary_ ? "strong" : "weak") << nl
        << "    CFL                " << CFL_ << nl
        << "    initialisation     " << initialisation_ << nl
        << "    source term        "
        << (analytical().hasSource() ? "yes" : "no") << nl
        << "    unknown nodes      " << nUnknowns_ << " of "
        << dual_.nNodes() << nl
        << "    inflow edges       " << nInflowEdges << " of "
        << dual_.nBoundaryEdges() << " boundary edges (A . n_B < 0)" << nl
        << "Initial RMS residual of Res_j/V_j over the unknowns = "
        << residual_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spaceTimeModels::vertexCentred::~vertexCentred()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spaceTimeModels::vertexCentred::calcResidual
(
    const scalarField& u,
    scalarField& res
) const
{
    if (u.size() != dual_.nNodes())
    {
        FatalErrorInFunction
            << "Field size " << u.size() << " is not the number of nodes "
            << dual_.nNodes() << exit(FatalError);
    }

    const pointField& p = dual_.points();
    const edgeList& edges = dual_.edges();
    const vectorField& n = dual_.n();

    res.setSize(dual_.nNodes());
    res = 0;

    // LSQ gradient for the linear reconstruction (VC-2)
    vectorField gradU(dual_.nNodes(), Zero);

    if (linearReconstruction_)
    {
        gradU = dual_.gradient(u);
    }

    // Edge loop: Phi_jk with n_jk from j to k
    forAll(edges, edgei)
    {
        const label j = edges[edgei].start();
        const label k = edges[edgei].end();

        const scalar magN = mag(n[edgei]);
        const scalar lambda = (n[edgei]/magN) & A();
        const vector d = p[k] - p[j];

        scalar uL = u[j];
        scalar uR = u[k];

        if (linearReconstruction_)
        {
            uL += 0.5*(gradU[j] & d);
            uR -= 0.5*(gradU[k] & d);
        }

        // Upwind flux. Eq. 97 of Tufillaro et al. is printed with
        // -0.5 |nHat . a| (uL - uR), which is downwind: the sign here is
        // the upwind one (CLAUDE.md section 7)
        const scalar phi =
            0.5*lambda*(uL + uR) - 0.5*mag(lambda)*(uR - uL);

        res[j] += phi*magN;
        res[k] -= phi*magN;
    }

    if (!strongBoundary_)
    {
        addWeakBoundaryFluxes(u, res);
    }

    addSource(res);

    // Strong boundary treatment: the boundary nodes are not unknowns
    forAll(res, nodei)
    {
        if (!isUnknown_[nodei])
        {
            res[nodei] = 0;
        }
    }
}


Foam::scalar Foam::spaceTimeModels::vertexCentred::rmsResidual
(
    const scalarField& res
) const
{
    const scalarField& V = dual_.V();
    scalar sumSqr = 0;

    forAll(res, nodei)
    {
        if (isUnknown_[nodei])
        {
            sumSqr += sqr(res[nodei]/V[nodei]);
        }
    }

    return Foam::sqrt(sumSqr/nUnknowns_);
}


bool Foam::spaceTimeModels::vertexCentred::evolve()
{
    const scalarField& V = dual_.V();

    // Stage 1: u* = u^n - (dTau/V) Res(u^n); res_ holds Res(u^n)
    scalarField uStar(u_);

    forAll(uStar, nodei)
    {
        if (isUnknown_[nodei])
        {
            uStar[nodei] -= dTau_[nodei]/V[nodei]*res_[nodei];
        }
    }

    // Stage 2: u^(n+1) = 0.5 (u* + u^n) - 0.5 (dTau/V) Res(u*)
    scalarField resStar(dual_.nNodes());
    calcResidual(uStar, resStar);

    forAll(u_, nodei)
    {
        if (isUnknown_[nodei])
        {
            u_[nodei] =
                0.5*(uStar[nodei] + u_[nodei])
              - 0.5*dTau_[nodei]/V[nodei]*resStar[nodei];
        }
    }

    // Residual of the new solution: the convergence measure, and Res(u^n)
    // of the next step
    calcResidual(u_, res_);
    residual_ = rmsResidual(res_);

    updatePointField();

    Info<< "vertexCentred: RMS residual of Res_j/V_j = " << residual_ << endl;

    if (!std::isfinite(residual_))
    {
        FatalErrorInFunction
            << "The pseudo-time iterations diverged: the RMS residual is "
            << residual_ << exit(FatalError);
    }

    return residual_ <= convergenceTolerance();
}


void Foam::spaceTimeModels::vertexCentred::writeFields()
{
    updatePointField();

    spaceTimeModel::writeFields();
}


// ************************************************************************* //
