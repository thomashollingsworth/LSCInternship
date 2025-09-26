

names = "density x-velocity y-velocity z-velocity pressure B_x B_y B_z psi"

do for [col=3:11] {
    set terminal pngcairo size 800,600
    set pm3d
    name = word(names, col-2)

    set output sprintf("%s_%s.png", filename,name)
    set xlabel "x"
    set ylabel "y"
    unset surface
    set view map
 

    
    splot filename using 1:2:col notitle
}




