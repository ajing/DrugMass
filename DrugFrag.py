'''
    Get drug fragment mass
'''

from ExtractChrom import ExtractSpec
from DrugMass import HighestPeaks, SelectPeaks
import pymzml
import Tkinter, tkFileDialog

def DrugFragMass(masslist, i):
    '''
        The function sort masslist first and then add mass for each element before index i
        for example if masslist = [439, 421, 312.2, 252, 170.8], i = 2
        result is [455, 437, 312.2, 252, 170.8]
    '''
    mass_list_mod = sorted(list(masslist), reverse = True)
    for j in range(i + 1):
        mass_list_mod[j] = mass_list_mod[j] + 16
    return mass_list_mod

def GetSumIntensityInOneSpec(mz_list, one_spec, tolerance = 0.2):
    '''
        Get the sum intensity for one spectrum for a list of mz in a collection of specs
    '''
    print_list = []
    max_int_dict = dict()
    spec = one_spec.spec
    for eachmz in mz_list:
        eachrange = (eachmz - tolerance, eachmz + tolerance)
        try:
            mz, intensity = HighestPeaks(SelectPeaks(spec["peaks"], eachrange))
        except Exception as e:
            #print e.message
            continue
        max_int_dict[eachmz] = {"max_int": intensity, "max_mz": mz, "max_time": spec["scan time"], "max_id": spec["id"]}
    sum_intensity = sum([value["max_int"] for key, value in max_int_dict.iteritems()])
    return sum_intensity

def SpecSumIntensity4MassList(mzlist, rt_time, exspec):
    specs = exspec.extractWithTime(rt_time)
    if not specs:
        # if specs empty
        tol   = 0.01
        specs = exspec.extractWithTimeRange(rt_time - tol, rt_time + tol)
        spec  = specs[len(specs) / 2]
    else:
        spec  = specs[0]
    return GetSumIntensityInOneSpec(mzlist, spec)

def main():
    #root = Tkinter.Tk()
    #root.withdraw()
    #ms_file = tkFileDialog.askopenfilename()
    ms_file   = "./Data/CCG224144MIDSample5minMS2.mzML"
    ms_file   = "./Data/CCG224144MIDSample5min.mzML"
    ms_file   = "./Data/5minMRM_Biotrans.mzML"
    mass_list = [423, 405, 296, 268, 171]  #  only for parent drug
    #mass_list = [439, 421, 312.2, 252, 170.8]
    exspec = ExtractSpec(ms_file)
    run = pymzml.run.Reader(ms_file, noiseThreshold = 100)
    #intensity_1 = SpecSumIntensity4MassList(mass_list, 5.83, exspec)
    #print 5.83, intensity_1
    for i in range(len(mass_list)):
        mass_list_mod = DrugFragMass(mass_list, i)
        #intensity_1 = SpecSumIntensity4MassList(mass_list_mod, 5.83, exspec)
        #print 5.83, mass_list_mod, intensity_1
        intensity_1 = SpecSumIntensity4MassList(mass_list_mod, 5.681, exspec)
        print 5.681, mass_list_mod, intensity_1

if __name__ == "__main__":
    #mass_list = [439, 421, 312.2, 252, 170.8]
    #print DrugFragMass(mass_list, 2)
    mass_list = [423, 405, 296, 268, 171]
    main()

