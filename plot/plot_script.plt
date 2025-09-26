

names = "density x-velocity y-velocity pressure"

do for [col=3:6] {
    set terminal pngcairo size 800,600
    set pm3d
    name = word(names, col-2)

    set output sprintf("%s_%s.png", filename,name)
    set xlabel "x"
    set ylabel "y"
    set zlabel sprintf("%s", name)
    unset colorbox
    set view 60, 30  # 3D view angles
    splot filename using 1:2:col notitle
}




