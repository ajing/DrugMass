#
#  Read ML file
#

import pymzml
import Tkinter, tkFileDialog

DEBUG = False

def InRange(spec, rtrange):
    try:
        rt_time = spec["scan time"]
        if not rtrange is None:
            if rt_time < rtrange[0]:
                return False
            elif rt_time > rtrange[1]:
                return False
    except:
        return False
    return rt_time

def HighestPeaks(peaklist):
    max_item= (0,0)
    for eachitem in peaklist:
        if eachitem[1] > max_item[1]:
            max_item = eachitem
    return max_item

def SelectPeaks(peaks, mzRange):
    return [ (mz,i) for mz, i in peaks if mzRange[0] <= mz <= mzRange[1] ]

def PlotRange(run, rt_time):
    p = pymzml.plot.Factory()
    n = 0
    for spec in run:
        if spec['ms level'] != 1:
            continue
        print spec["scan time"], abs(spec["scan time"] - rt_time)
        if abs(spec["scan time"] - rt_time) > 0.003:
            continue
        n = n + 1
        print "plot graph..."
        p.newPlot()
        p.add(spec.peaks, color=(200,00,00), style='circles')
        p.add(spec.centroidedPeaks, color=(00,00,00), style='sticks')
        p.add(spec.reprofiledPeaks, color=(00,255,00), style='circles')
        p.save( filename="output/plotAspect_%s.xhtml" %(n))

def Test():
    filename = "E165ug.mzML"
    run = pymzml.run.Reader(filename)
    print run[1701].extremeValues('i')
    print run[1700].extremeValues('i')
    print run[1702].extremeValues('i')
    print run[1702].peaks

if __name__ == "__main__":
    root = Tkinter.Tk()
    root.withdraw()
    inputfile = tkFileDialog.askopenfilename()
    run = pymzml.run.Reader(inputfile, noiseThreshold = 100)
    PlotRange(run, 5.681)
    #PlotRange(run, 6.38)
