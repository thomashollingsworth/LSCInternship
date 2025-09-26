

names = "density velocity pressure"

do for [col=2:4] {
    set terminal pngcairo size 800,600

    name = word(names, col-1)

    set output sprintf("%s_%s.png", filename,name)
    set xlabel "x"
    set ylabel sprintf("%s", name)
    
    plot filename using 1:col with lines notitle
}




