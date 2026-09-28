//-----------------------------------------------------------------------------
// License
//     This file is part of spaceTime4foam, licensed under GNU General Public
//     License <http://www.gnu.org/licenses/>.
//
// File
//     unstructuredSquare.geo
//
// Description
//     Unstructured triangle mesh of the space-time unit square (x, t) in
//     [0, 1] x [0, 1] for the medianDualIdentities test: gmsh 2D algorithm 6
//     (Frontal-Delaunay) with a uniform characteristic length lc (default
//     0.08, a few hundred triangles), extruded one layer in z (thickness 1)
//     into prisms, with the same patch names as the tutorial mesh (tStart,
//     tEnd, xLeft, xRight, frontAndBack). The mesh is deterministic: the
//     random seed of the mesher is fixed and one thread is used.
//
//     Usage:
//         gmsh -3 -format msh22 -o mesh.msh unstructuredSquare.geo
//
//-----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// Characteristic length
If (!Exists(lc))
    lc = 0.08;
EndIf

// Deterministic meshing
General.NumThreads = 1;
Mesh.RandomSeed = 1;
Mesh.Algorithm = 6;
Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Extrusion thickness in z
thickness = 1.0;

// Tolerance for the bounding-box patch selection
eps = 1e-6;

Rectangle(1) = {0, 0, 0, 1, 1};

Extrude {0, 0, thickness} { Surface{1}; Layers{1}; Recombine; }

tStart[] = Surface In BoundingBox{-eps, -eps, -eps, 1 + eps, eps, thickness + eps};
tEnd[] = Surface In BoundingBox{-eps, 1 - eps, -eps, 1 + eps, 1 + eps, thickness + eps};
xLeft[] = Surface In BoundingBox{-eps, -eps, -eps, eps, 1 + eps, thickness + eps};
xRight[] = Surface In BoundingBox{1 - eps, -eps, -eps, 1 + eps, 1 + eps, thickness + eps};
front[] = Surface In BoundingBox{-eps, -eps, -eps, 1 + eps, 1 + eps, eps};
back[] = Surface In BoundingBox{-eps, -eps, thickness - eps, 1 + eps, 1 + eps, thickness + eps};

Physical Surface("tStart") = {tStart[]};
Physical Surface("tEnd") = {tEnd[]};
Physical Surface("xLeft") = {xLeft[]};
Physical Surface("xRight") = {xRight[]};
Physical Surface("frontAndBack") = {front[], back[]};
Physical Volume("spaceTime") = {Volume{:}};

//-----------------------------------------------------------------------------
