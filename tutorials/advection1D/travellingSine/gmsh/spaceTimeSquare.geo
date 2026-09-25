//-----------------------------------------------------------------------------
// License
//     This file is part of spaceTime4foam, licensed under GNU General Public
//     License <http://www.gnu.org/licenses/>.
//
// File
//     spaceTimeSquare.geo
//
// Description
//     Space-time unit square: mesh x is physical x in [0, 1] and mesh y is
//     time t in [0, 1]. The square is split into N x N squares, each cut into
//     two right triangles by a diagonal in ONE direction (transfinite "Right"
//     arrangement: every diagonal joins (x_i, t_j) and (x_i+1, t_j+1)).
//     The triangles are extruded one layer in z, so every cell is a prism.
//
//     The z thickness is 1 (nondimensional), so a cell volume equals its
//     triangle area.
//
//     The OpenCASCADE kernel is used because it places the transfinite
//     nodes at exactly i/N; the built-in kernel gives errors of about 1e-12.
//     Patches are selected by bounding box, so they do not depend on the
//     order of the surfaces returned by Extrude.
//
//     Usage (N defaults to 8):
//         gmsh -3 -setnumber N 16 -format msh22 -o mesh.msh spaceTimeSquare.geo
//
//-----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// Number of squares in x and in t
If (!Exists(N))
    N = 8;
EndIf

// Extrusion thickness in z
thickness = 1.0;

// Tolerance for the bounding-box patch selection
eps = 1e-6;

// Corner points: (x, t, z)
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};    // t = 0
Line(2) = {2, 3};    // x = 1
Line(3) = {3, 4};    // t = 1
Line(4) = {4, 1};    // x = 0

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// N + 1 nodes on each side gives N x N squares
Transfinite Curve{1, 2, 3, 4} = N + 1;

// "Right" puts all diagonals in the same direction
Transfinite Surface{1} = {1, 2, 3, 4} Right;

// Extrude one layer. Recombine turns the extruded triangles into prisms and
// the extruded lines into quads; the triangles themselves are kept.
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
