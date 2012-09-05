"""Module to help writing and plotting datafiles"""
import subprocess
import os

imagetype = "png"
#imagetype = "eps"


def plotwitherrorbarsnames(basename, data_errors, shifts, autoscale=True):
    """ plots data with errorbars, data_errors must be a dict with the
    name given by the key and the values are a tuple with (data,errors)"""

    if not os.path.exists(os.path.dirname(basename)):
        os.makedirs(os.path.dirname(basename))

    writedataerrorsfile(basename, data_errors, shifts)
    plotfile = open(basename + ".plt", "w")

    names = data_errors.keys()

    if(autoscale):
        maxval = 0.0
        for _, data_error in data_errors.iteritems():
            maxval = max(max(data_error[0].values()), maxval)
        maxrange = maxval + maxval * 0.2

        minval = maxval
        for k, data_error in data_errors.iteritems():
            minval = min(min(data_error[0].values()), maxval)
        minrange = minval - minval * 0.2

    else:
        maxrange = 2.0
        minrange = 0.0

    plotfile.write("set yrange[0:%f] \n" % maxrange)
    #plotfile.write("set yrange[%f:%f] \n" % (minrange, maxrange))
    plotfile.write("set xrange[0:%d] \n" % len(shifts))
    if imagetype is "eps":
        plotfile.write("set terminal postscript eps color enhanced \n")
    else:
        plotfile.write("set terminal png \n")
    plotfile.write("set ylabel \"%s\" \n" % basename)
    plotfile.write("set output \"%s.%s\" \n" % (basename, imagetype))
    plotfile.write("set style line 1 lt 1 lc 1 pt 9 ps 2  lw 2 \n")
    plotfile.write("set style line 2 lt 1 lc 3 pt 7 ps 2  lw 2 \n")
    for index, name in enumerate(names):
        plotfile.write("%s" % "plot" if index is 0 else ", ")
        plotfile.write(" \"%s.out\" using 1:%d:%d with "
                       "yerrorbars linestyle %d title \"%s\" "
                       % (basename, index * 2 + 2, index * 2 + 3, index + 1, name))
    plotfile.write("\n")
#    subprocess.check_call(["gnuplot ", " %s.plt" % basename])
    plotfile.flush()
    plotfile.close()
    subprocess.call(["sync"])
    returncode = subprocess.call(["gnuplot", basename + ".plt"])
    if returncode:
        raise OSError("gnuplot failed to plot" + basename + ".plt")
    else:
        print("plotted " + basename + ".plt" + " successfully")


def writedataerrorsfile(basename, data_errors, shifts):
    """Write a data file so it may be plotted"""
    writefile = open(basename + ".out", "w")

    names = data_errors.keys()
    writefile.write("# time \t \"%s \" \n" % "\t".join(names))

    for shift in shifts:
        writefile.write(str(shift))
        for _, values in data_errors.iteritems():
            writefile.write(",\t")
            data, errors = values
            writefile.write(" {0}".format(float(data[shift]), "7.3g"))
            writefile.write(",\t")
            writefile.write(" {0}".format(errors[shift], "7.3g"))
        writefile.write("\n")
    writefile.flush()
    writefile.close()
    print("wrote " + basename + ".out ")


def writedatafile(basename, data, indexes):
    """Write a data file so it may be plotted"""
    writefile = open(basename + ".out", "w")

    names = data.keys()
    writefile.write("# time \t \"%s \" \n" % "\t".join(names))

    for index in indexes:
        writefile.write(str(index))
        for _, values in data.iteritems():
            writefile.write(",\t")
            writefile.write(" {0}".format(values[index], "7.3g"))
        writefile.write("\n")
    writefile.flush()
    writefile.close()
    print("wrote " + basename + ".out ")


def plotwithnames(basename, datas, shifts):
    """ plots data with errorbars, datas must be a dict with the
    name given by the key and the values are a tuple with (data, errors)"""

    writedatafile(basename, datas, shifts)
    plotfile = open(basename + ".plt", "w")

    names = datas.keys()

    maxval = 0.0
    for _, data in datas.iteritems():
        maxval = max(max(data), maxval)
    maxrange = maxval + maxval * 0.1

    minval = maxval
    for k, data in datas.iteritems():
        minval = min(min(data), minval)
    minrange = minval - minval * 0.1

#    plotfile.write("set yrange[0:%f] \n" % maxrange)
    plotfile.write("set yrange[%f:%f] \n" % (minval, maxval))
    plotfile.write("set xrange[0:%d] \n" % len(shifts))
    plotfile.write("set nokey\n")
    if imagetype is "eps":
        plotfile.write("set terminal postscript eps color enhanced \n")
    else:
        plotfile.write("set terminal png \n")
    plotfile.write("set ylabel \"%s\" \n" % basename)
    plotfile.write("set output \"%s.%s\" \n" % (basename, imagetype))
    plotfile.write("set style line 1 lt 1 lc 1 pt 9 ps 2  lw 2 \n")
    plotfile.write("set style line 2 lt 1 lc 3 pt 7 ps 2  lw 2 \n")
    for index, name in enumerate(names):
        plotfile.write("%s" % "plot" if index is 0 else ", ")
        plotfile.write(" \"%s.out\" using 1:%d with"
                       "linespoints linestyle %d title \"%s\" "
                       % (basename, index + 2, index + 1, name))

    plotfile.write("\n")
#    subprocess.check_call(["gnuplot ", " %s.plt" % basename])
    plotfile.flush()
    plotfile.close()
    subprocess.call(["sync"])
    returncode = subprocess.call(["gnuplot", basename + ".plt"])
    if returncode:
        raise OSError("gnuplot failed to plot" + basename + ".plt")
    else:
        print("plotted " + basename + ".plt" + " successfully")
