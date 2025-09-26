

names = "density v_x v_y v_z pressure B_x B_y B_z"

do for [col=2:9] {
    set terminal pngcairo size 800,600

    name = word(names, col-1)

    set output sprintf("%s_%s.png", filename,name)
    set xlabel "x"
    set ylabel sprintf("%s", name)
    
    plot filename using 1:col with lines notitle
}




