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
    Foam::spaceTimeModels::vertexCentred

Description
    Vertex-centred (node-centred, edge-based) finite volume (VCFV)
    space-time discretisation of

        div_st(A u) = f,    A = (a, 1, 0)

    on the median-dual mesh of the triangulated front patch
    (medianDualMesh), following CLAUDE.md section 7 and Tufillaro, Williams
    and Nishikawa (2026), Section 5.1, Eq. 92-103. The unknowns are the
    nodal values u_j.

    Residual (one edge loop over the edges (j, k), j < k, with the
    directed-area vector n_jk from j to k and nHat = n_jk/|n_jk|):

        uL = u_j + 0.5 grad(u)_j . (x_k - x_j)      (VC-2, Eq. 93)
        uR = u_k - 0.5 grad(u)_k . (x_k - x_j)      (VC-2, Eq. 94)
        (VC-1: uL = u_j, uR = u_k)
        Phi_jk = 0.5 (nHat . A)(uL + uR) - 0.5 |nHat . A| (uR - uL)
        Res_j += Phi_jk |n_jk|,   Res_k -= Phi_jk |n_jk|
        Res_j -= f_j V_j        (only if the analytical solution has one)

    grad(u) is the unweighted least-squares gradient of medianDualMesh.
    The upwind term has the sign above (upwind); Eq. 97 of the paper is
    printed with -0.5 |nHat . a| (uL - uR), which is the downwind sign.

    The source f_j = analyticalSolution::source(x_j) (S = A & grad(u) from
    the exact gradient) is added by its own function, so that a later
    model problem changes the source, or adds terms to the edge loop,
    without restructuring the residual.

    Boundary treatment (boundaryTreatment):
    - weak: every node is an unknown. Each boundary edge with outward
      normal n_B (|n_B| = edge length) contributes to each of its two nodes
      m, with area |n_B|/2 and nHat = n_B/|n_B|, the upwind flux
          Phi = 0.5 (nHat . A)(u_m + u_b) - 0.5 |nHat . A| (u_b - u_m)
      between the nodal value u_m and the boundary state u_b. The edge is
      inflow if A . n_B < 0, and then u_b is the analytical value at node
      m; otherwise (outflow, A . n_B >= 0) u_b = u_m and Phi = (nHat . A) u_m.
      Inflow and outflow are decided per edge from the sign of A . n_B, not
      from the patch names. A corner node, e.g. (x, t) = (1, 0) between
      tStart (inflow) and xRight (outflow), gets one contribution of each
      kind: the analytical value through its inflow edge and its own value
      through its outflow edge.
    - strong: the boundary nodes (every node of a boundary edge) are held
      at the analytical value, set at initialisation and never changed, and
      their residual is set to zero. Only the interior nodes are unknowns.

    No sign depends on the face winding: n_jk and n_B are oriented by
    medianDualMesh with explicit dot-product tests.

    Pseudo-time: two-stage Runge-Kutta, Eq. 101-103,
        u*      = u^n - (dTau/V) Res(u^n)
        u^(n+1) = 0.5 (u* + u^n) - 0.5 (dTau/V) Res(u*)
    The paper's Eq. 102 does not show the argument of Res; it is evaluated
    at u*, as in the standard two-stage (Heun, SSP) Runge-Kutta scheme.
    The local pseudo-time step is
        dTau_j = CFL V_j / D_j,   D_j = sum_k 0.5 |n_jk . A|
    over the edge neighbours k of j. With weak boundary treatment a
    boundary node also adds 0.5 |n_B . A| / 2 for each of its boundary
    edges. This is a deviation from the paper (which uses strong Dirichlet
    data and has no boundary fluxes); it changes only the convergence path,
    not the converged solution.

    Iterations: one Runge-Kutta step per call to evolve(), i.e. per
    runTime iteration. After each step the residual of the new solution is
    computed and
        residual() = sqrt(sum_U (Res_j/V_j)^2 / nU)
    over the set U of the nU unknown nodes (all nodes for weak, interior
    nodes for strong); this is the "final residual" written by
    spaceTime4Foam, and it is the residual of the solution that is written.
    Converged when residual() <= convergenceTolerance. A non-finite
    residual is a fatal error (diverged).

    Output: u is written as a pointScalarField; the front and back mesh
    points of each node get the same nodal value. Its patch fields are
    calculated, and empty on frontAndBack. Every mesh point must be a front
    or a back point (one-layer mesh).

    Settings in constant/spaceTimeProperties (all required):
    \verbatim
    spaceTimeModel       vertexCentred;
    a                    1;
    convergenceTolerance 1e-10;
    analyticalSolution { type travellingSine; k 1; }

    vertexCentredCoeffs
    {
        reconstruction     linear;    // none (VC-1) or linear (VC-2)
        boundaryTreatment  weak;      // weak or strong
        CFL                0.5;
        initialisation     zero;      // zero or analytical
    }
    \endverbatim

SourceFiles
    vertexCentred.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentred_H
#define vertexCentred_H

#include "spaceTimeModel.H"
#include "medianDualMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace spaceTimeModels
{

/*---------------------------------------------------------------------------*\
                        Class vertexCentred Declaration
\*---------------------------------------------------------------------------*/

class vertexCentred
:
    public spaceTimeModel
{
    // Private data

        //- The vertexCentredCoeffs sub-dictionary
        const dictionary coeffs_;

        //- Linear reconstruction (VC-2) if true, none (VC-1) if false
        const bool linearReconstruction_;

        //- Strong (true) or weak (false) boundary treatment
        const bool strongBoundary_;

        //- CFL number of the pseudo-time step
        const scalar CFL_;

        //- Initialisation: zero or analytical
        const word initialisation_;

        //- Median-dual geometry of the front triangles
        const medianDualMesh dual_;

        //- True for the nodes whose value is solved for
        boolList isUnknown_;

        //- Number of unknown nodes
        label nUnknowns_;

        //- Analytical value at every node
        scalarField uExact_;

        //- Source term f_j V_j at every node (zero if there is no source)
        scalarField sourceV_;

        //- True for inflow boundary edges (A . n_B < 0)
        boolList isInflowEdge_;

        //- Local pseudo-time step dTau_j
        scalarField dTau_;

        //- Nodal solution u_j
        scalarField u_;

        //- Residual Res_j of the current solution
        scalarField res_;

        //- Output field: u at the mesh points
        pointScalarField uPoint_;

        //- RMS of Res_j/V_j over the unknowns, for the current solution
        scalar residual_;


    // Private Member Functions

        //- Read a word from coeffs_ and check it is one of the options
        word readOption(const word& key, const wordList& options) const;

        //- Set isUnknown_ and nUnknowns_
        void calcUnknowns();

        //- Set uExact_ and isInflowEdge_
        void calcBoundaryData();

        //- Set sourceV_ = f_j V_j
        void calcSource();

        //- Set dTau_
        void calcPseudoTimeStep();

        //- Set the initial u_
        void initialise();

        //- Subtract the source f_j V_j from the residual
        void addSource(scalarField& res) const;

        //- Add the weak boundary fluxes to the residual
        void addWeakBoundaryFluxes
        (
            const scalarField& u,
            scalarField& res
        ) const;

        //- Copy u_ to the front and back points of uPoint_
        void updatePointField();

        //- No copy construct
        vertexCentred(const vertexCentred&) = delete;

        //- No copy assignment
        void operator=(const vertexCentred&) = delete;


public:

    //- Runtime type information
    TypeName("vertexCentred");


    // Constructors

        //- Construct from Time
        vertexCentred(Time& runTime);


    //- Destructor
    virtual ~vertexCentred();


    // Member Functions

        // Access

            //- Return the median-dual geometry
            const medianDualMesh& dual() const
            {
                return dual_;
            }

            //- Return the nodal solution
            const scalarField& u() const
            {
                return u_;
            }

            //- Return the output point field
            const pointScalarField& uPoint() const
            {
                return uPoint_;
            }

            //- True for linear reconstruction (VC-2)
            bool linearReconstruction() const
            {
                return linearReconstruction_;
            }

            //- True for strong boundary treatment
            bool strongBoundary() const
            {
                return strongBoundary_;
            }

            //- True for the nodes whose value is solved for
            const boolList& isUnknown() const
            {
                return isUnknown_;
            }

            //- Number of unknown nodes
            label nUnknowns() const
            {
                return nUnknowns_;
            }

            //- Local pseudo-time step dTau_j
            const scalarField& dTau() const
            {
                return dTau_;
            }


        // Evaluation

            //- Compute the residual Res_j of the nodal values u (see the
            //  class description), including the source and the boundary
            //  treatment
            void calcResidual(const scalarField& u, scalarField& res) const;

            //- RMS of Res_j/V_j over the unknown nodes
            scalar rmsResidual(const scalarField& res) const;


        // Evolution

            //- Perform one two-stage Runge-Kutta step; return true once
            //  converged
            virtual bool evolve();

            //- Return the RMS residual of the current solution
            virtual scalar residual() const
            {
                return residual_;
            }


        // Write

            //- Write the fields at the current time (iteration)
            virtual void writeFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace spaceTimeModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
