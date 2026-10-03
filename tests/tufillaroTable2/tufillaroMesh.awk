#------------------------------------------------------------------------------
# License
#     This file is part of spaceTime4foam, licensed under GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# File
#     tufillaroMesh.awk
#
# Description
#     Write mesh n of the 2D mesh family of Tufillaro, Williams and Nishikawa
#     (2026), Section 6.1, Table 1 and Fig. 8, as a one-layer extruded gmsh
#     MSH 2.2 file for gmshToFoam:
#       - mesh 0: the unit square (x, y) in [0, 1]^2 split by both diagonals
#         into 4 triangles, with 5 nodes (the corners and the centre);
#       - mesh n + 1: every triangle of mesh n split into 4 at its edge
#         midpoints (the 3 corner triangles and the middle one).
#     Mesh n has 4^(n + 1) triangles, 2^n edges on each side and h = 1/2^n.
#     For n >= 2 this is not the same as crossing each square of a grid by
#     both diagonals: the node sets agree, the connectivity does not.
#     All coordinates are dyadic, so they are exact in binary.
#
#     y plays the role of t. The triangles are extruded one layer in z
#     (thickness 1) into prisms. Patches (physical names): tStart (y = 0),
#     tEnd (y = 1), xLeft (x = 0), xRight (x = 1) and frontAndBack (the
#     triangles at z = 0 and z = 1). Every triangle is counter-clockwise in
#     (x, y), so every prism has its bottom triangle at z = 0 ordered
#     counter-clockwise seen from z > 0.
#
#     Usage:
#         awk -v n=3 -f tufillaroMesh.awk > mesh3.msh
#
#------------------------------------------------------------------------------

# midpoint(a, b): the node at the midpoint of nodes a and b, created once
function midpoint(a, b,    key)
{
    key = (a < b) ? a "," b : b "," a

    if (!(key in mid))
    {
        nNodes++
        X[nNodes] = 0.5*(X[a] + X[b])
        Y[nNodes] = 0.5*(Y[a] + Y[b])
        mid[key] = nNodes
    }

    return mid[key]
}

# addTriangle(a, b, c): append a (counter-clockwise) triangle to the new list
function addTriangle(a, b, c)
{
    nNew++
    NA[nNew] = a
    NB[nNew] = b
    NC[nNew] = c
}

BEGIN {
    if (n == "" || n !~ /^[0-9]+$/)
    {
        print "tufillaroMesh.awk: set the mesh number, -v n=<0, 1, ...>" \
            > "/dev/stderr"
        exit 1
    }

    # Mesh 0: corners 1-4 and centre 5
    nNodes = 5
    X[1] = 0; Y[1] = 0
    X[2] = 1; Y[2] = 0
    X[3] = 1; Y[3] = 1
    X[4] = 0; Y[4] = 1
    X[5] = 0.5; Y[5] = 0.5

    nTri = 4
    TA[1] = 1; TB[1] = 2; TC[1] = 5
    TA[2] = 2; TB[2] = 3; TC[2] = 5
    TA[3] = 3; TB[3] = 4; TC[3] = 5
    TA[4] = 4; TB[4] = 1; TC[4] = 5

    # Refinement: every triangle into 4 at its edge midpoints
    for (level = 1; level <= n; level++)
    {
        nNew = 0
        delete mid

        for (t = 1; t <= nTri; t++)
        {
            a = TA[t]; b = TB[t]; c = TC[t]
            mab = midpoint(a, b)
            mbc = midpoint(b, c)
            mca = midpoint(c, a)

            addTriangle(a, mab, mca)
            addTriangle(mab, b, mbc)
            addTriangle(mca, mbc, c)
            addTriangle(mab, mbc, mca)
        }

        nTri = nNew
        for (t = 1; t <= nTri; t++)
        {
            TA[t] = NA[t]; TB[t] = NB[t]; TC[t] = NC[t]
        }
    }

    # Boundary edges: edges of one triangle only, kept in the triangle's
    # (counter-clockwise) direction
    for (t = 1; t <= nTri; t++)
    {
        v[1] = TA[t]; v[2] = TB[t]; v[3] = TC[t]
        for (i = 1; i <= 3; i++)
        {
            a = v[i]
            b = v[i % 3 + 1]
            key = (a < b) ? a "," b : b "," a
            edgeCount[key]++
        }
    }

    # Second pass in triangle order, so the output order does not depend on
    # the order of awk's associative arrays
    nBoundary = 0
    for (t = 1; t <= nTri; t++)
    {
        v[1] = TA[t]; v[2] = TB[t]; v[3] = TC[t]
        for (i = 1; i <= 3; i++)
        {
            a = v[i]
            b = v[i % 3 + 1]
            key = (a < b) ? a "," b : b "," a
            if (edgeCount[key] != 1)
            {
                continue
            }

            # Patch from the edge midpoint (exact dyadic coordinates)
            xm = 0.5*(X[a] + X[b])
            ym = 0.5*(Y[a] + Y[b])
            if (ym == 0) { patch = 1 }
            else if (ym == 1) { patch = 2 }
            else if (xm == 0) { patch = 3 }
            else if (xm == 1) { patch = 4 }
            else
            {
                print "tufillaroMesh.awk: boundary edge " a " " b \
                    " is not on the square" > "/dev/stderr"
                exit 1
            }

            nBoundary++
            BA[nBoundary] = a
            BB[nBoundary] = b
            BP[nBoundary] = patch
        }
    }

    # MSH 2.2: back node i + nNodes is front node i at z = 1
    print "$MeshFormat"
    print "2.2 0 8"
    print "$EndMeshFormat"
    print "$PhysicalNames"
    print 6
    print "2 1 \"tStart\""
    print "2 2 \"tEnd\""
    print "2 3 \"xLeft\""
    print "2 4 \"xRight\""
    print "2 5 \"frontAndBack\""
    print "3 6 \"spaceTime\""
    print "$EndPhysicalNames"

    print "$Nodes"
    print 2*nNodes
    for (i = 1; i <= nNodes; i++)
    {
        printf "%d %.17g %.17g 0\n", i, X[i], Y[i]
    }
    for (i = 1; i <= nNodes; i++)
    {
        printf "%d %.17g %.17g 1\n", i + nNodes, X[i], Y[i]
    }
    print "$EndNodes"

    # Elements: elm-number elm-type 2 physical elementary node-list
    print "$Elements"
    print 3*nTri + nBoundary
    id = 0
    for (t = 1; t <= nTri; t++)
    {
        id++
        printf "%d 6 2 6 1 %d %d %d %d %d %d\n", id, TA[t], TB[t], TC[t], \
            TA[t] + nNodes, TB[t] + nNodes, TC[t] + nNodes
    }
    for (t = 1; t <= nTri; t++)
    {
        id++
        printf "%d 2 2 5 2 %d %d %d\n", id, TA[t], TC[t], TB[t]
        id++
        printf "%d 2 2 5 2 %d %d %d\n", id, TA[t] + nNodes, TB[t] + nNodes, \
            TC[t] + nNodes
    }
    for (i = 1; i <= nBoundary; i++)
    {
        id++
        printf "%d 3 2 %d %d %d %d %d %d\n", id, BP[i], 2 + BP[i], BA[i], \
            BB[i], BB[i] + nNodes, BA[i] + nNodes
    }
    print "$EndElements"
}

#------------------------------------------------------------------------------
