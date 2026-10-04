#------------------------------------------------------------------------------
# License
#     This file is part of spaceTime4foam, licensed under GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# Script
#     plotConvergence.gp
#
# Description
#     gnuplot script for the S5 benchmark figures. Every figure is drawn
#     from the tables written by makeTables (<scheme>_<family>.dat and
#     orders_<scheme>_<family>.dat) only. Figures whose tables are missing
#     (for a partial sweep) are skipped.
#
#     Every scheme keeps the same colour and point type in every figure.
#
#     Usage: gnuplot -e "resultsDir='results'" plotConvergence.gp
#            (figures in <resultsDir>/figures)
#
#------------------------------------------------------------------------------

if (!exists("resultsDir")) { resultsDir = "results" }
figureDir = resultsDir."/figures"
system("mkdir -p '".figureDir."'")

schemes = "CC-1 CC-2 CC-2-LS CC-2-EX CC-2-LS-EX VC-1 VC-2"
colours = "#2a78d6 #eb6834 #1baf7a #eda100 #e87ba4 #008300 #4a3aa7"
# Filled circle, square, triangle, inverted triangle, diamond for CCFV;
# open circle and open square for VCFV
pointTypes = "7 5 9 11 13 6 4"

# Neutral colours: text, axes, grid and reference lines
textColour = "#222222"
axisColour = "#444444"
gridColour = "#dddddd"
refColour = "#777777"

# Table columns (see the header of every <scheme>_<family>.dat)
colN = 1; colH = 2; colUnknowns = 4
colL2 = 6; colLinf = 7; colL2Int = 9; colL2Bnd = 12; colL2T = 15
colIter = 24; colWall = 27
# Orders table: N_coarse N_fine p_L1_all p_L2_all p_Linf_all p_L1_int
# p_L2_int ...
colOrdersNFine = 2; colOrdersL2Int = 7

table(i, family) = sprintf("%s/%s_%s.dat", resultsDir, word(schemes, i), \
    family)
ordersTable(i, family) = sprintf("%s/orders_%s_%s.dat", resultsDir, \
    word(schemes, i), family)
colour(i) = word(colours, i)
pointType(i) = word(pointTypes, i) + 0

# available(family, indices): the indices (1-7, in scheme order) of the
# given schemes that have a table for the family
available(family, indices) = system(sprintf( \
    "for i in %s; do s=$(echo '%s' | cut -d ' ' -f $i);" \
    . " [ -f '%s/'$s'_%s.dat' ] && printf '%%s ' $i; done", \
    indices, schemes, resultsDir, family))


familyTitle(family) = (family eq "leftStructured") ? \
    "structured left-diagonal triangles" : \
    (family eq "leftPerturbed") ? \
    "perturbed left-diagonal triangles (seed 12345)" : family

# Note added to the titles of the perturbed figures
probeNote(family) = (family eq "leftPerturbed") ? \
    "\nCC-2-LS-EX omitted: unstable here (see ccLsExPerturbedProbe.dat)" : ""

set terminal pngcairo size 900,650 enhanced font "Helvetica,14" \
    background rgb "#ffffff"
set encoding utf8
set datafile missing "-"
set samples 2

set border 3 lc rgb axisColour lw 1
set tics nomirror out textcolor rgb textColour
set grid xtics ytics lc rgb gridColour lw 1 dt 1
# The legend is outside the plot, so it never hides data
set key outside right top vertical Left reverse
set key textcolor rgb textColour spacing 1.3 samplen 2.5 noopaque nobox
set title textcolor rgb textColour
set xlabel textcolor rgb textColour
set ylabel textcolor rgb textColour
set style line 100 lc rgb refColour lw 1 dt (6, 4)

hTics = 'set xtics ("1/256" 1./256, "1/128" 1./128, "1/64" 1./64,' \
    . ' "1/32" 1./32, "1/16" 1./16, "1/8" 1./8)'
nTics = 'set xtics ("8" 8, "16" 16, "32" 32, "64" 64, "128" 128,' \
    . ' "256" 256)'

# Series style of scheme s (an index 1-7, as a string)
seriesStyle = 'lw 2 pt pointType(s + 0) ps 1.6 lc rgb colour(s + 0)'

# minFine: smallest value of column col on the finest mesh (N = nMax) over
# the schemes in list (sets yRef)
minFine = 'yRef = 1e30; do for [sc in list] { ' \
    . 'stats table(sc + 0, family) using (column(colN) == nMax ? ' \
    . 'column(col) : NaN) nooutput; ' \
    . 'if (STATS_min < yRef) { yRef = STATS_min } }'

# Reference slope p through (hMin, 0.3 yRef), below all the data on the
# finest mesh, drawn for hMin <= h <= 1/8 and labelled at its coarse end
refSlope(x, p) = 0.3*yRef*(x/hMin)**p
refLabels = 'unset label; ' \
    . 'set label 1 sprintf("slope %g", p1) at 1./8, refSlope(1./8, p1) ' \
    . 'offset 0.6, 0 left textcolor rgb refColour; ' \
    . 'set label 2 sprintf("slope %g", p2) at 1./8, refSlope(1./8, p2) ' \
    . 'offset 0.6, 0 left textcolor rgb refColour'
refLines = '[x = hMin:1./8] "+" using 1:(refSlope($1, p1)) ' \
    . 'with lines ls 100 notitle, ' \
    . '[x = hMin:1./8] "+" using 1:(refSlope($1, p2)) ' \
    . 'with lines ls 100 notitle'

# h axis: the finest to the coarsest mesh, with room for the slope labels
hAxis = 'set logscale xy; eval hTics; set xrange [hMin/1.3:1./8*2.2]; ' \
    . 'set format x "%g"; set format y "10^{%L}"; set xlabel "h = 1/N"'

mainFamilies = "leftStructured leftPerturbed"
allIndices = "1 2 3 4 5 6 7"


do for [family in mainFamilies] {

    list = available(family, allIndices)
    if (strlen(list) == 0) { continue }

    # Finest mesh of the family (from the first table)
    stats table(word(list, 1) + 0, family) using colN nooutput
    nMax = STATS_max
    hMin = 1.0/nMax
    p1 = 1
    p2 = 2

    # 1. L2 (all) against h
    col = colL2
    eval minFine
    eval refLabels
    eval hAxis
    set output sprintf("%s/L2_vs_h_%s.png", figureDir, family)
    set title sprintf("L2 error over all unknowns against h,\n%s%s", \
        familyTitle(family), probeNote(family))
    set ylabel "L2 error (all cells/nodes)"
    plot for [s in list] table(s + 0, family) \
            using colH:colL2 with linespoints @seriesStyle \
            title word(schemes, s + 0), \
        @refLines

    # 3. Linf (all) against h
    col = colLinf
    eval minFine
    eval refLabels
    set output sprintf("%s/Linf_vs_h_%s.png", figureDir, family)
    set title sprintf("Linf error over all unknowns against h,\n%s%s", \
        familyTitle(family), probeNote(family))
    set ylabel "Linf error (all cells/nodes)"
    plot for [s in list] table(s + 0, family) \
            using colH:colLinf with linespoints @seriesStyle \
            title word(schemes, s + 0), \
        @refLines

    # 5. Primary t = T L2 against h
    col = colL2T
    eval minFine
    eval refLabels
    set output sprintf("%s/L2T_vs_h_%s.png", figureDir, family)
    set title sprintf("L2 error at t = T (written tEnd value for CCFV," \
        . " nodal for VCFV) against h,\n%s%s", familyTitle(family), \
        probeNote(family))
    set ylabel "L2 error at t = T (primary value)"
    plot for [s in list] table(s + 0, family) \
            using colH:colL2T with linespoints @seriesStyle \
            title word(schemes, s + 0), \
        @refLines

    # 4. Interior against boundary L2 for CC-2, CC-2-EX and VC-2
    listSubset = available(family, "2 4 7")
    if (strlen(listSubset) > 0) {
        listAll = list
        list = listSubset
        col = colL2Int
        eval minFine
        yInt = yRef
        col = colL2Bnd
        eval minFine
        if (yInt < yRef) { yRef = yInt }
        list = listAll
        p1 = 1.5
        eval refLabels
        set output sprintf("%s/interiorBoundary_vs_h_%s.png", figureDir, \
            family)
        set title sprintf("L2 error over interior (solid) and boundary" \
            . " (dashed) unknowns against h,\n%s", familyTitle(family))
        set ylabel "L2 error (interior or boundary cells/nodes)"
        plot for [s in listSubset] table(s + 0, family) \
                using colH:colL2Int with linespoints @seriesStyle \
                title word(schemes, s + 0)." interior", \
            for [s in listSubset] table(s + 0, family) \
                using colH:colL2Bnd with linespoints @seriesStyle \
                dt (8, 5) title word(schemes, s + 0)." boundary", \
            @refLines
        p1 = 1
    }

    # 2. L2 (all) against the number of unknowns
    unset label
    set output sprintf("%s/L2_vs_unknowns_%s.png", figureDir, family)
    set title sprintf("L2 error over all unknowns against the number of" \
        . " unknowns,\n%s%s", familyTitle(family), probeNote(family))
    set xtics autofreq
    set xrange [*:*]
    set format x "10^{%L}"
    set xlabel "number of unknowns (cells: 2 N^2; nodes: (N + 1)^2)"
    set ylabel "L2 error (all cells/nodes)"
    plot for [s in list] table(s + 0, family) \
            using colUnknowns:colL2 with linespoints @seriesStyle \
            title word(schemes, s + 0)

    # 6a. Iterations against N
    set output sprintf("%s/iterations_vs_N_%s.png", figureDir, family)
    set title sprintf("Iterations to convergence (tolerance 1e-10)" \
        . " against N,\n%s%s", familyTitle(family), probeNote(family))
    eval nTics
    set xrange [8/1.3:nMax*1.3]
    set format x "%g"
    set format y "%g"
    set xlabel "N (N x N squares, h = 1/N)"
    set ylabel "iterations"
    plot for [s in list] table(s + 0, family) \
            using colN:colIter with linespoints @seriesStyle \
            title word(schemes, s + 0)

    # 6b. Wall time against N
    set output sprintf("%s/wallTime_vs_N_%s.png", figureDir, family)
    set title sprintf("Wall time of the solver loop against N,\n%s%s", \
        familyTitle(family), probeNote(family))
    set format y "10^{%L}"
    set ylabel "wall time of the solver loop [s]"
    plot for [s in list] table(s + 0, family) \
            using colN:colWall with linespoints @seriesStyle \
            title word(schemes, s + 0)

    # 6c. L2 (all) against wall time (efficiency)
    set output sprintf("%s/L2_vs_wallTime_%s.png", figureDir, family)
    set title sprintf("L2 error over all unknowns against wall time" \
        . " (efficiency),\n%s%s", familyTitle(family), probeNote(family))
    set xtics autofreq
    set xrange [*:*]
    set format x "10^{%L}"
    set xlabel "wall time of the solver loop [s]"
    set ylabel "L2 error (all cells/nodes)"
    plot for [s in list] table(s + 0, family) \
            using colWall:colL2 with linespoints @seriesStyle \
            title word(schemes, s + 0)
    set format x "%g"
}


# 7. Interior L2 orders of CC-2 and VC-2, structured against perturbed

listS = available("leftStructured", "2 7")
listP = available("leftPerturbed", "2 7")
if (strlen(listS) > 0 && strlen(listP) > 0) {
    unset label
    unset logscale y
    set logscale x
    set output sprintf("%s/interiorOrders_structuredVsPerturbed.png", \
        figureDir)
    set title "Observed order of the interior L2 error against N" \
        . " (pair N/2 to N),\nCC-2 and VC-2, structured (solid) and" \
        . " perturbed (dashed) left-diagonal triangles"
    eval nTics
    set xrange [16/1.3:256*1.3]
    set yrange [0:3.5]
    set ytics 0.5
    set format y "%.1f"
    set xlabel "N of the finer mesh (pair N/2 to N)"
    set ylabel "observed order p of the interior L2 error"
    set arrow 1 from graph 0, first 2 to graph 1, first 2 nohead ls 100
    set arrow 2 from graph 0, first 3 to graph 1, first 3 nohead ls 100
    set label 1 "order 2" at graph 0.12, first 2 offset 0, -0.7 \
        textcolor rgb refColour
    set label 2 "order 3" at graph 0.12, first 3 offset 0, 0.6 \
        textcolor rgb refColour
    plot for [s in listS] ordersTable(s + 0, "leftStructured") \
            using colOrdersNFine:colOrdersL2Int with linespoints \
            @seriesStyle title word(schemes, s + 0)." structured", \
        for [s in listP] ordersTable(s + 0, "leftPerturbed") \
            using colOrdersNFine:colOrdersL2Int with linespoints \
            @seriesStyle dt (8, 5) title word(schemes, s + 0)." perturbed"
    unset arrow
    unset label
    set autoscale y
    set ytics autofreq
}


# 8. Right-diagonal (aligned) family: an illustration only

family = "rightStructured"
list = available(family, allIndices)
if (strlen(list) > 0) {
    stats table(word(list, 1) + 0, family) using colN nooutput
    nMax = STATS_max
    hMin = 1.0/nMax
    p1 = 1
    p2 = 2
    col = colL2
    eval minFine
    eval refLabels
    eval hAxis
    set output sprintf("%s/L2_vs_h_rightStructured_illustration.png", \
        figureDir)
    set title "Characteristic-alignment ILLUSTRATION (not a main result):" \
        . "\nL2 error over all unknowns against h, right diagonals" \
        . " parallel to A = (1, 1)"
    set ylabel "L2 error (all cells/nodes)"
    plot for [s in list] table(s + 0, family) \
            using colH:colL2 with linespoints @seriesStyle \
            title word(schemes, s + 0), \
        @refLines
}

#------------------------------------------------------------------------------
