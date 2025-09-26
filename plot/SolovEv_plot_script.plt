# Define filename

outputfile = sprintf("%s.png", filename)

# Set terminal and output
set terminal pngcairo size 1000,800 font "Arial,14"
set output outputfile
set xlabel "R"
set ylabel "Z"

# 2D top-down view (map)
set view map

# Enable pm3d for filled contour
set pm3d at b

# Enable contour lines at the base
set contour base

# Set number of contour levels
set cntrparam levels 35
set pm3d at b flush begin clip



# Don't draw 3D surface mesh
unset surface

# Hide key
unset key

# Set aspect ratio

set size ratio kappa

# Plot filled contour + contour lines overlay
splot filename with pm3d, filename with lines lc rgb "black" lw 2

# Finish output (close file)
set output
