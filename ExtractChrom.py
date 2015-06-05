'''
    Save spectrum info to a data structure, so it's easy to retrieve spect for a special time
'''

import pymzml

class SpecBasic:
    def __init__(self, rt, index, mslevel = 1):
        # here the program assume the retention time is a float number
        self._rt   = float(rt)
        self._idx  = index
        mslevel    = 1
        self._spec = dict()

    @property
    def rtime(self):
        return self._rt

    @property
    def index(self):
        return self._idx

    @property
    def spec(self):
        return self._spec

    @spec.setter
    def spec(self, spec):
        self._spec = spec

    def __str__(self):
        return "retention time: %s, index: %s" %(self._rt, self._idx)


class SpecDict(dict):
    def __init__(self):
        self._dict = dict()

    def __getitem__(self, time):
        if float(time).is_integer():
            return self._dict[int(time)]
        else:
            matchlist = []
            dec_len = len(str(time).split(".")[-1])
            for specbasic in self._dict[int(time)]:
                if (specbasic.rtime - time) * 10 ** dec_len // 1 == 0:
                    matchlist.append(specbasic)
            return matchlist

    def __setitem__(self, time, specbasic):
        if not isinstance(specbasic, SpecBasic):
            raise Exception("SpectDict only accept SpecBasic")
        try:
            self._dict[int(time)].append(specbasic)
        except:
            self._dict[int(time)] = [specbasic]

    def __str__(self):
        return "number of spec: %s" % (len(self._dict.keys()))

    def getRange(self, s_time, f_time):
        matchlist = []
        for int_time in range(int(s_time), int(f_time) + 1):
            if not int_time in self._dict:
                continue
            for specbasic in self._dict[int_time]:
                if specbasic.rtime > s_time and specbasic.rtime < f_time:
                    matchlist.append(specbasic)
        return matchlist

class ExtractSpec:
    def __init__(self, filename):
        self.start_time  = 0
        self.end_time    = 0
        self.interval    = 0
        self.run, self.specdict = self.setup(filename)
        #print self.specdict

    def setup(self, filename):
        # set up basic data structure
        run = pymzml.run.Reader(filename)
        specdict = SpecDict()
        #self.start_time = run[1]["scan time"]
        #self.interval   = run[2]["scan time"] - run[1]["scan time"]
        # first spec flag
        flag_f = 1
        for spectrum in run:
            if spectrum['ms level'] == 2:
                if flag_f:
                    self.start_time = spectrum["scan time"]
                    flag_f = 0
                #print max(spectrum.i), min(spectrum.i)
                print spectrum['scan time']
                specbasic = SpecBasic(spectrum['scan time'], spectrum['id'])
                specbasic.spec["peaks"] = spectrum.peaks
                specbasic.spec["scan time"] = spectrum['scan time']
                specbasic.spec["id"] = spectrum['id']
                specdict[spectrum['scan time']] = specbasic
                if spectrum["scan time"] > self.end_time:
                    self.interval = spectrum["scan time"] - self.end_time
                    self.end_time = spectrum["scan time"]
        return run, specdict

    def extractWithTime(self, time):
        return self.specdict[time]

    def extractWithTimeRange(self, s_time, f_time):
        return self.specdict.getRange(s_time, f_time)


def ExtractIonChrom(run):
    timeDependentIntensities = []
    #MASS_2_FOLLOW = 1402
    MASS_2_FOLLOW = 1403
    print "mass 2 follow:", MASS_2_FOLLOW
    for spectrum in run:
        print spectrum['id']
    #for spectrum in run:
    #    if spectrum['ms level'] == 1:
    #        matchList = spectrum.hasPeak(MASS_2_FOLLOW)
    #        if matchList != []:
    #            print "matched"
    #            for mz,I in matchList:
    #                if I > 100:
    #                    timeDependentIntensities.append( [ spectrum['scan time'], I , mz ])
    for rt, i, mz in timeDependentIntensities:
        print('{0:5.3f} {1:13.4f}       {2:10}'.format( rt, i, mz ))



def ExtractTest():
    exspec = ExtractSpec("./4tRNA1_102009.mzML")
    # extract spectrums for specific time
    print "for time 1"
    specs  = exspec.extractWithTime(1)
    for spec in specs:
        print spec
    print "for time 1.1"
    specs  = exspec.extractWithTime(1.1)
    for spec in specs:
        print spec
    print "for time 1.1 to 1.5"
    specs  = exspec.extractWithTimeRange(1.1, 1.5)
    for spec in specs:
        print spec

if __name__ == "__main__":
    ExtractTest()
