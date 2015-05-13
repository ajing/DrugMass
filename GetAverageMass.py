'''
    Get the average mass for 3 specs
'''

import Tkinter, tkFileDialog
from ExtractChrom import ExtractSpec
from DrugMass import GetMaxPeakInSpecs, HighestPeaks, SelectPeaks
import pymzml

def GetAverage(specs):
    for spec in specs:
        print spec
    return sum(specs) / len(specs)

def WithinTime(spec, rt_range):
    return spec['scan time'] > rt_range[0] and spec['scan time'] < rt_range[1]

def GetAveSpecForRT(mass_list, run, exspec, rt_time, t_interval = 0.01):
    # get the average spec for a special retention time
    # t interval 0.01 is 3 spectrum
    all_spec = []
    for spectrum in run:
        if spectrum['ms level'] == 1:
            if WithinTime(spectrum, [rt_time - t_interval * 1.0 / 2, rt_time + t_interval * 1.0 / 2]):
                all_spec.append(spectrum.deRef())
    spec_sum = None
    for spec in all_spec:
        if spec_sum is None:
            spec_sum = spec
        else:
            spec_sum += spec
    print "the number of spectrum selected: ", len(all_spec)
    spec_ave = spec_sum / len(all_spec)
    return spec_ave

def FindPeaksInSpec(mass_list, spec, tolerance = 0.2):
    mz_peak_list = []
    max_int_dict = dict()
    for eachmz in mass_list:
        eachrange = (eachmz - tolerance, eachmz + tolerance)
        try:
            mz, intensity = HighestPeaks(SelectPeaks(spec.peaks, eachrange))
        except Exception as e:
            print e.message
            continue
        max_int_dict[eachmz] = {"max_int": intensity, "max_mz": mz}
    if max_int_dict and all([value["max_int"] > 0 for key, value in max_int_dict.iteritems()]):
        sum_intensity = sum([value["max_int"] for key, value in max_int_dict.iteritems()])
        mz_peak_list = [ [max_int_dict[x]["max_int"], max_int_dict[x]["max_mz"]] for x in max_int_dict.keys()]
        print max_int_dict
        print mz_peak_list

def main():
    root = Tkinter.Tk()
    root.withdraw()
    #ms_file = tkFileDialog.askopenfilename()
    ms_file   = "./Data/CCG224144MIDSample5minMS2.mzML"
    #ms_file   = "./Data/CCG224144MIDSample5min.mzML"
    #mass_list = [423.3, 268.2]
    mass_list = [439, 421, 312.2, 252, 170.8]
    rt_time = 5.681
    average_interval = 0.5
    mz_tolerance     = 0.3
    exspec = ExtractSpec(ms_file)
    run = pymzml.run.Reader(ms_file, noiseThreshold = 100)
    spec_ave = GetAveSpecForRT(mass_list, run, exspec, rt_time, average_interval)
    FindPeaksInSpec(mass_list, spec_ave, mz_tolerance)

if __name__ == "__main__":
    main()
