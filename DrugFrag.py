'''
    Get drug fragment mass
'''

from ExtractChrom import ExtractSpec
from DrugMass import HighestPeaks, SelectPeaks, DRange
import pymzml
import Tkinter, tkFileDialog
import pprint
import string

from csv import DictReader

MOD_FILE = "mod_list.txt"

def DrugFragMass(masslist, i, mod_value):
    '''
        The function sort masslist first and then add mass for each element before index i
        for example if masslist = [439, 421, 312.2, 252, 170.8], i = 2
        result is [455, 437, 312.2, 252, 170.8]
    '''
    mass_list_mod = sorted(list(masslist), reverse = True)
    for j in range(i + 1):
        #mass_list_mod[j] = mass_list_mod[j] + 16
        #mass_list_mod[j] = mass_list_mod[j] + mod_value
        mass_list_mod[j] = mass_list_mod[j] + mod_value
    return mass_list_mod

def GetSumIntensityInOneSpec(mz_list, one_spec, tolerance = 0.11):
    '''
        Get the sum intensity for one spectrum for a list of mz in a collection of specs
    '''
    print_list = []
    max_int_dict = dict()
    spec = one_spec.spec
    mz, max_int = HighestPeaks(spec["peaks"])
    for eachmz in mz_list:
        eachrange = (eachmz - tolerance, eachmz + tolerance)
        try:
            mz, intensity = HighestPeaks(SelectPeaks(spec["peaks"], eachrange))
        except Exception as e:
            #print e.message
            continue
        max_int_dict[eachmz] = {"max_int": intensity, "max_mz": mz, "max_time": spec["scan time"], "max_id": spec["id"], "max_abundance": intensity / max_int}
    #pprint.pprint(max_int_dict)
    sum_intensity = sum([value["max_int"] for key, value in max_int_dict.iteritems()])
    sum_abundance = sum([value["max_abundance"] for key, value in max_int_dict.iteritems()])
    return sum_intensity, sum_abundance, spec["scan time"]

def SpecSumIntensity4MassList(mzlist, rt_time, exspec):
    specs = exspec.extractWithTime(rt_time)
    #print "current time:", rt_time
    if not specs:
        # if specs empty
        tol   = 0.01
        while True:
            specs = exspec.extractWithTimeRange(rt_time - tol, rt_time + tol)
            if len(specs) > 0:
                break
            else:
                tol = tol + 0.01
        spec  = specs[len(specs) / 2]
    else:
        spec  = specs[0]
    return GetSumIntensityInOneSpec(mzlist, spec)

def MetCandidate(i, n, mod_part):
    uppercase = string.ascii_uppercase
    mod_part_p= "(" + mod_part + ")"
    molecule  = uppercase[:n]
    molecule  = list(molecule)
    molecule[i] = molecule[i] + mod_part_p
    return "-".join(molecule)

def main():
    #root = Tkinter.Tk()
    #root.withdraw()
    #ms_file = tkFileDialog.askopenfilename()
    ms_file   = "./Data/CCG224144MIDSample5minMS2.mzML"
    ms_file   = "./Data/CCG224144MIDSample5min.mzML"
    ms_file   = "./Data/5minMRM_Biotrans.mzML"
    mass_list = [423, 405, 296, 268, 171]  #  only for parent drug
    exspec = ExtractSpec(ms_file)
    run = pymzml.run.Reader(ms_file, noiseThreshold = 100)

    mod_dict  = ParseModShift(MOD_FILE)
    print mod_dict
    mass_list = [423, 296, 268, 171]  #  only for parent drug

    for mod in mod_dict:
        mass_shift = float(mod_dict[mod]["MassShift"])
        time_set   = float(mod_dict[mod]["RetentionTime"])
        RunForEachMod(mass_list, mod, mass_shift, exspec, time_set)

def ParseModShift(filename):
    mod_dict = dict()
    for line in DictReader(open(filename), delimiter = " "):
        mod_dict[line["Modification"]] = line
    return mod_dict

def RunForEachMod(mass_list, mod_type, mod_val, exspec, time_set):
    abund_dict = dict()
    # initialize abund_dict
    for i in range(len(mass_list)):
        mass_list_mod = DrugFragMass(mass_list, i, mod_val)
        abund_dict[tuple(mass_list_mod)] = {"abundance": 0, "rt_time": 0, "intensity": 0}
    for rt_time in DRange(time_set - 0.5, time_set + 0.5, exspec.interval):
        ## ignore the spec out of rtrange
        #if (not rtrange is None) and (rt_time < min(rtrange) or rt_time > max(rtrange)):
        #    continue
        for i in range(len(mass_list)):
            mass_list_mod = DrugFragMass(mass_list, i, mod_val)
            try:
                sum_intensity, sum_abundance, scan_time = SpecSumIntensity4MassList(mass_list_mod, rt_time, exspec)
            except Exception, e:
                print e
            if sum_intensity > abund_dict[tuple(mass_list_mod)]["intensity"]:
            #if abundance_1 > abund_dict[tuple(mass_list_mod)]["abundance"]:
                abund_dict[tuple(mass_list_mod)] = {"abundance": sum_abundance, "meta_cand": MetCandidate(i, len(mass_list), mod_type), "rt_time": scan_time, "intensity": sum_intensity}
    pprint.pprint(abund_dict)


def single_time(rt_time, modtype):
    #ms_file   = "./Data/CCG224144MIDSample5minMS2.mzML"
    #ms_file   = "./Data/CCG224144MIDSample5min.mzML"
    ms_file   = "./Data/5minMRM_Biotrans.mzML"

#    root = Tkinter.Tk()
#    root.withdraw()
#    ms_file = tkFileDialog.askopenfilename()
    # initialize abund_dict
    abund_dict = dict()
    for i in range(len(mass_list)):
        mass_list_mod = DrugFragMass(mass_list, i, mod_val)
        try:
            intensity_1, abundance_1, scan_time = SpecSumIntensity4MassList(mass_list_mod, rt_time, exspec)
        except Exception, e:
            print e
        abund_dict[tuple(mass_list_mod)] = {"abundance": abundance_1, "meta_cand": MetCandidate(i, len(mass_list), modtype), "rt_time": scan_time, "intensity": intensity_1}
    #print abund_dict
    pprint.pprint(abund_dict)



    exspec = ExtractSpec(ms_file)
    run = pymzml.run.Reader(ms_file, noiseThreshold = 100)



if __name__ == "__main__":
    main()
    #single_time(5.681)
    #modtype = "-2H"
    #single_time(5.83, modtype)
    #single_time(7.373)

