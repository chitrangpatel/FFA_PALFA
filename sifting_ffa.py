#!/usr/bin/env python
import sys
import re
import os
import copy
import numpy as Num
import matplotlib.pyplot as plt

"""
Slightly modified version of sifting.py (PRESTO). 
Modifications that have been applied are commented in capital letters.
"""

fund_re = re.compile("^\d")
harms_re = re.compile("^[ ]\d")
DM_re = re.compile("DM(\d+\.\d{2})")
# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
r_err = 1.1
# Longest period candidates to consider (s)
long_period = 15.0
# Shortest period candidates to consider (s)
short_period = 0.0005
# Ignore candidates with a sigma (from incoherent power summation) less than this
sigma_threshold = 6.0
# Ignore candidates with a coherent power less than this
c_pow_threshold = 100.0
# Ignore any candidates where at least one harmonic does not exceed this power
harm_pow_cutoff = 8.0

# If the birds file works well, the following shouldn't
# be needed at all...
#                (ms, err)
known_birds_p = []
#                (Hz, err)
known_birds_f = []
from presto import candidate_sigma

def remove_harmonics(candlist, *args, **kwargs):
    """Remove harmonics. The candlist is modified
        **in-place**.

        Note: This function is defined to maintain support
            for old code. It simply calls the 
            'remove_harmonics' method of candlist.

        Inputs:
            ** All arguments are passed onto the
            'remove_harmonics' method of candlist.

        Output:
            candlist: The modified candidate list.
    """
    candlist.remove_harmonics(*args, **kwargs)
    return candlist


def cmp_freq(self, other):
    return cmp(self.r, other.r)

def cmp_sigma(self, other):
    retval = -cmp(self.snr, other.snr)
    if retval==0:
        return -cmp(self.snr, other.snr)
    else:
        return retval


def cmp_snr(self, other):
    retval = -cmp(self.snr, other.snr)
    if retval==0:
        return -cmp(self.snr, other.snr)
    else:
        return retval


def cmp_dms(self, other):
    return cmp(float(self[0]), float(other[0]))

def write_candlist(candlist, *args, **kwargs):
    candlist.to_file(*args, **kwargs)


def candlist_from_candfile(filename, trackbad=False, trackdupes=False):
    print filename
    candfile = open(filename, 'r')
# CHANGED FEW THINGS 
    # First identify the length of the observation searched
    for line in candfile:
        if line.startswith(" Number of bins in the time series"):
            numsamp = int(line.split()[-1])
        if line.startswith(" Width of each time series bin (sec)"):
            dt = float(line.split()[-1])
	if line.startswith(" Dispersion measure (cm-3 pc)"):
	    dm = float(line.split()[-1])
    print "DM of the cands:", dm
    tobs = numsamp * dt
    # Go back to the start of the file to read the candidates
    candfile.seek(0)

    cands = []
    candnums = []
    current_goodcandnum = 0
    for line in candfile:
        # Identify the candidates in the top of the file
        if fund_re.match(line):
# CHANGED HERE FOR FFA
            split_line = line.split()
	    candnum = int(split_line[0])
            p = float(split_line[1])       # Spin period in sec
	    f = 1.0/p    		   # Spin freq in hz
	    snr = float(split_line[2])
	    width = float(split_line[3])
	    bin = tobs*f
	    

            # Add it to the candidates list
	    if DM_re.search(filename) == None: DMstr = '9999999.9999'
            else: DMstr = DM_re.search(filename).groups()[0]
            cands.append(Candidate(candnum,p, snr, width,bin, dm,DMstr, filename, tobs))
            continue
     
    candfile.close()
    return Candlist(cands, trackbad=trackbad, trackdupes=trackdupes)


def read_candidates(filenms, prelim_reject=True, track=False):
    """Read in accelsearch candidates from the test ACCEL files.
        Return a Candlist object of Candidate instances.

        Inputs:
            filenms: A list of files to read candidates from.
            prelim_reject: If True, perform preliminary rejection of
                candidates. (Default: True)
            track: If True, keep track of bad/duplicate candidates.
                (Default: False)

    """
    candlist = Candlist(trackbad=track, trackdupes=track)
# CHANGED THE READING HERE 
    fi = open(filenms,'r')
    fns=fi.readlines()
    fns=[item.strip() for item in fns]
    numfiles = len(fns)
    ii = 0
    if fns:
        print "\nReading candidates from %d files...." % len(filenms)
        for filenm in fns:
            curr_candlist = candlist_from_candfile(filenm, trackbad=track, trackdupes=track)
            if prelim_reject:
                curr_candlist.default_rejection()
            candlist.extend(curr_candlist)
            sys.stdout.write(" Read %d of %d files (%d cands)\r" % (ii+1, numfiles, len(candlist)))
            sys.stdout.flush()
	    ii+=1
        print "\nDone"
    else:
        print "Error:  There are no candidate files to read!"
    return candlist

def remove_duplicate_candidates(candlist, *args, **kwargs):
    """Remove duplicate candidates. The candlist is modified
        **in-place**.

        Note: This function is defined to maintain support
            for old code. It simply calls the 
            'remove_duplicate_candidates' method of candlist.

        Inputs:
            ** All arguments are passed onto the
            'remove_duplicate_candidates' method of candlist.

        Output:
            candlist: The modified candidate list.
    """
    candlist.remove_duplicate_candidates(*args, **kwargs)
    return candlist


def cmp_freq(self, other):
    return cmp(self.r, other.r)

def print_sift_globals():
    print "r_err =", r_err
    print "short_period =", short_period     
    print "long_period =", long_period     
    print "sigma_threshold =", sigma_threshold 
    print "c_pow_threshold =", c_pow_threshold 
    print "harm_pow_cutoff =", harm_pow_cutoff 
    print "known_birds_p =", known_birds_p   
    print "known_birds_f =", known_birds_f


class Candidate(object):
# CHANGED sigma -> SNR
# REMOVED POWER STUFF AND HARM 
    def __init__(self,candnum, p, snr, dt ,bin, dm, DMstr, filename, T):
        self.path, self.filename = os.path.split(filename)
	self.candnum = int(candnum)
	self.p = p
        self.snr = snr
	self.dt = dt
        self.f = 1.0/p
        self.T = T
	self.r = bin
        self.DMstr = DMstr
        self.DM = float(dm)
        self.hits = []
        self.note = ""

    def add_as_hit(self, other):
        self.hits.extend(other.hits)

    def __str__(self):
        cand = self.filename + ':	' + `self.candnum`
	print cand
        return "%-65s   %7.2f  %6.2f   %s   %7.5f   "% (cand, self.p*1000, self.snr, self.DM, self.dt*1000)
# REMOVED HARMS_TO_SNR


class Candlist(object):
    def __init__(self, cands=None, trackbad=False, trackdupes=False):
        if cands is None:
            self.cands = []
        else:
            self.cands = cands
        self.trackbad = trackbad # Should we keep track of bad candidates
        self.trackdupes = trackdupes # Should we keep track of duplicates
        # Set default badlists
# REMOVED HARMONIC RELATED STUFF
        self.badlists = {'knownbirds': [], \
                         'longperiod': [], \
                         'shortperiod': [], \
                         'threshold': [], \
                         'dmproblem': []}
        self.duplicates = []

    def __iter__(self):
        return iter(self.cands)

    def __getitem__(self, key):
        return self.cands[key]

    def __delitem__(self, key):
        del(self.cands[key])

    def sort(self, *args, **kwargs):
        self.cands.sort(*args, **kwargs)

# CHANGED USEFREQS DEFAULT
    def plot_summary(self, usefreqs=False):
        """Produce a plot summarizing the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Get all candidates and sort by sigma
        allcands = self.get_all_cands()
        sigmas = Num.array([c.snr for c in allcands])
        isort = sigmas.argsort()
        sigmas = sigmas[isort]
        if usefreqs:
            xdata = Num.array([c.f for c in allcands])[isort]
            xlabel = "Freq (Hz)"
            xscale = "log"
        else:
            xdata = Num.array([c.p for c in allcands])[isort]
            xlabel = "Period (s)"
            xscale = "loglin"
        dms = Num.array([c.DM for c in allcands])[isort]

        # Plot the all candidates 
# REMOVED HARM STUFF 
        scatt = plt.scatter(xdata, dms, s=sigma_to_size(sigmas), marker='o', alpha=0.7, zorder=-1) 
        plt.set_cmap("Spectral") 
  
        # Add colorbar 
        fmtr = matplotlib.ticker.FuncFormatter(lambda x, pos: "%d" % 2**x)
        cax = plt.axes((0.18, 0.06, 0.67, 0.035))
        cb = plt.colorbar(scatt, cax=cax, ticks=(0,1,2,3,4), format=fmtr, \
                            orientation="horizontal")
        cb.set_label("Num harmonics summed") 
        
        plt.axes(ax) # Set scatter plot's axes as current
        plt.xscale(xscale)
        plt.xlabel(xlabel)
        mindm = Num.min(dms)
        maxdm = Num.max(dms)
        dmrange = Num.ptp(dms)
        plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)
        plt.ylabel(r"DM (pc cm$^{-3}$)") 
        if not usefreqs:
            plt.gca().xaxis.set_ticks(Num.concatenate((\
                                        Num.logspace(-4,0,4, endpoint=False), \
                                        Num.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(Num.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(xdata)/5.0), \
                        min(long_period+0.5, max(xdata)+0.5))
        ax.format_coord = lambda x,y: "x=%g, y=%g" % (x,y)
        return fig


    def remove_harmonics(self, verbosity=1):
        """Remove the candidates that are lower significance harmonics
            of other candidates from the candlist.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Note:  should probably put the harmonics into the fundamental as hits (use sets)
        numremoved = 0
        self.cands.sort(cmp_snr)
        f_err = r_err/self.cands[0].T
        if verbosity >= 1:
            print "\nSearching for duplicate harmonics..."
        ii = 0
        while 1:
            fundcand = self.cands[ii]
            jj = len(self.cands) - 1
            zapj = 0
            while 1:
                harmcand = self.cands[jj]
                if zapj:  print "Hey!"
                for factor in Num.arange(1.0, 17.0):
                    if Num.fabs(fundcand.f - harmcand.f*factor) < f_err*factor:
                        zapj = 1
                        harmstr = "1/%dth" % factor
                    elif Num.fabs(fundcand.f - harmcand.f/factor) < f_err/factor:
                        zapj = 1
                        if factor==2.0:
                            harmstr = "%dnd" % factor
                        else:
                            harmstr = "%dth" % factor
                    if zapj:
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%s) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        harmstr, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        break
                # Check a few other common ratios
                for numer,denom in zip([3.0, 5.0, 2.0, 4.0, 5.0, \
                                        3.0, 5.0, 2.0, 3.0, 4.0],
                                       [2.0, 2.0, 3.0, 3.0, 3.0, \
                                        4.0, 4.0, 5.0, 5.0, 5.0]):
                    factor = numer/denom
                    if Num.fabs(fundcand.f-harmcand.f*factor) < f_err*factor:
                        if verbosity >= 2:
                            print "Removing %s:%d (%.2f Hz) because it is " \
                                    "a harmonic (%d/%dth) of %s:%d (%.2f Hz)" % \
                                    (harmcand.filename, \
                                        harmcand.candnum, \
                                        harmcand.f, \
                                        denom, \
                                        numer, \
                                        fundcand.filename, \
                                        fundcand.candnum, \
                                        fundcand.f)
                        harmstr = "%d/%dth" % (denom, numer)
                        zapj = 1
                        break
                if zapj:
                    harmcand.note = "This candidate (P=%.4f s, DM=%.2f) is " \
                                    "a harmonic (%s) of %s:%d " \
                                    "(P=%.4f s, DM=%.2f)." % \
                                (harmcand.p, harmcand.DM, harmstr, \
                                    fundcand.filename, fundcand.candnum, \
                                    fundcand.p, fundcand.DM)
                    numremoved += 1
                    self.mark_as_bad(jj, 'harmonic')
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (harmcand.filename, harmcand.candnum, jj)
                        print "    %s" % harmcand.note
                    zapj = 0
                jj -= 1
                if jj == ii:
                    break
            ii += 1
            if ii >= len(self.cands) - 1:
                break
        if verbosity >= 1:
            print "Removed a total of %d harmonics.\n" % numremoved

    def plot_goodcands(self, usefreqs=False):
        """Produce a plot highlighting good candidates as selected by
            the sifiting performed.

            Input:
                usefreqs: If True, the horizontal axis will use
                    frequency. If False, use period.
            
            Output:
                fig: A matplotlib figure instance.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,8)) 
        ax = plt.axes((0.08, 0.18, 0.87, 0.80)) 
        
        # Plot candidates
        labels = []
        candlists = []
        for key in self.badlists:
            labels.append(key.title())
            candlists.append(self.badlists[key])
        candlists.append(self.cands)
        labels.append('Good cands')
        colours = ['#FF0000', '#800000', '#008000', '#00FF00', \
                    '#00FFFF', '#0000FF', '#FF00FF', '#800080', 'r']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
        zorders = [-2, -2, -2, -2, -2, -2, -2, -2, 0]
        sizes = [10, 10, 10, 10, 10, 10, 10, 10, 50]
        fixedsizes = [1, 1, 1, 1, 1, 1, 1, 1, 0]
        lws = [1,1,1,1,1,1,1,1,1,1]
        ecs = ['none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'k']
        alphas = [1,1,1,1,1,1,1,1,0.7]
        handles = []
        for cands, colour, marker, zorder, size, fixedsize, lw, alpha, ec in \
                zip(candlists, colours, markers, zorders, sizes, fixedsizes, lws, alphas, ecs):
            sigmas = []
            dms = []
            xdata = []
            for c in cands:
                sigmas.extend([h.sigma for h in c.hits])
                dms.extend([h[0] for h in c.hits])
                if usefreqs:
                    xval = c.f
                else:
                    xval = c.p
                xdata.extend([xval]*len(c.hits))
            sigmas = Num.array(sigmas)
            dms = Num.array(dms)
            xdata = Num.array(xdata)

            isort = sigmas.argsort()
            sigmas = sigmas[isort]
            dms = dms[isort]
            xdata = xdata[isort]
            if usefreqs:
                xlabel = "Freq (Hz)"
                xscale = "log"
            else:
                xlabel = "Period (s)"
                xscale = "loglin"
            
            # Plot the candidates
            if fixedsize:
                plt.scatter(xdata, dms, s=size, lw=lw, edgecolors=ec, \
                            c=colour, marker=marker, alpha=alpha, zorder=zorder)
            else:
                plt.scatter(xdata, dms, s=sigma_to_size(sigmas), lw=lw, edgecolors=ec, \
                            c=colour, marker=marker, alpha=alpha, zorder=zorder)
            handles.append(plt.scatter([], [], s=size, c=colour, \
                                    marker=marker, alpha=0.7))

        fig.legend(handles, labels, 'lower center', \
                        prop={'size':'x-small'}, ncol=4)

        plt.xscale(xscale) 
        plt.xlabel(xlabel) 
        mindm = Num.min(dms)
        maxdm = Num.max(dms)
        dmrange = Num.ptp(dms)
        plt.ylim(mindm-0.1*dmrange, maxdm+0.1*dmrange)
        plt.ylabel(r"DM (pc cm$^{-3}$)")
        if not usefreqs:
            plt.gca().xaxis.set_ticks(Num.concatenate((\
                                        Num.logspace(-4,0,4, endpoint=False), \
                                        Num.linspace(1,15,8))))
            plt.gca().xaxis.set_ticks(Num.logspace(-4,0,40), minor=True)
            plt.gca().xaxis.set_ticklabels([r"10$^{-4}$", r"10$^{-3}$", \
                        r"10$^{-2}$", r"10$^{-1}$", "1", "3", "5", "7", \
                        "9", "11", "13", "15"])
            plt.xlim(max(short_period/5.0, min(xdata)/5.0), \
                        min(long_period+0.5, max(xdata)+0.5))
        return fig

    def mark_as_bad(self, icand, badlistname):
        cand = self.cands.pop(icand)
        if self.trackbad:
            badlist = self.badlists.setdefault(badlistname, [])
            badlist.append(cand)

    def mark_as_duplicate(self, icand):
        cand = self.cands.pop(icand)
        if self.trackdupes:
            self.duplicates.append(cand)

    def get_all_cands(self):
        cands = self.get_all_goodcands()
        return self.get_all_goodcands() + self.get_all_badcands()

    def get_all_goodcands(self):
        return self.cands + self.duplicates

    def get_all_badcands(self):
        cands = []
        for key in self.badlists.keys():
            cands += self.badlists[key]
        return cands


    
# REMOVED THE FUNCTIONS REJECT_LONG/SHORT PERIODS

    def reject_knownbirds(self, known_birds_f=[], known_birds_p=[]):
        """Find and remove candidates conincident with known birds.

            Inputs:
                known_birds_f: A list of tuples containing bad frequencies
                    and widths. The tuples should contain
                        (<bad freq (Hz)>, <one-sided width (Hz)>)
                    (Default: Globally defined "known_birds_f")
                known_birds_p: A list of tuples containing bad peridocities
                    and widths. The tuples should contain
                        (<bad freq (ms)>, <one-sided width (ms)>)
                    (Default: Globally defined "known_birds_p")

            Outputs:
                None
        """
        if known_birds_f is None:
            known_birds_f = globals()['known_birds_f']
        if known_birds_p is None:
            known_birds_p = globals()['known_birds_p']
        for ii in reversed(range(len(self.cands))):
            cand = self.cands[ii]
            known_bird = 0
            for bird, err in known_birds_f:
                if (Num.fabs(cand.f-bird) < err):
                    known_bird = 1
                    cand.note = "Freq (%.2f Hz) is within %g Hz " \
                                    "of a known birdie centred at %.2f Hz" % \
                                    (cand.f, err, bird)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue
            for bird, err in known_birds_p:
                if (Num.fabs(cand.p*1000.0-bird) < err):
                    known_bird = 1
                    cand.note = "Period (%.2f ms) is within %g ms " \
                                    "of a known birdie centred at %.2f ms" % \
                                    (cand.f*1000, err, bird)
                    break
            if known_bird:
                self.mark_as_bad(ii, 'knownbirds')
                continue

# REMOVED REJECT_TRESHOLD (that is done in the FFA)
# REMOVED THE 'REJECT_' HARMONIC RELATED 

# DEFAULT REJECTION: REMOVED DELETED REJECT_FUNC        
    def default_rejection(self):
        """Run all rejection methonds with default arguments.

            Inputs:
                None

            Outputs:
                None
        """
        self.reject_knownbirds()

    def remove_duplicate_candidates(self, verbosity=1):
        """Remove lower-significance 'duplicate' (i.e. same period)
            candidates from a list of candidates.  For the highest
            significance candidate, include a list of the DMs (and SNRs)
            of all the other detections.

            Inputs:
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        if verbosity >= 1:
            print "  Sorting the %d candidates by frequency..." % \
                        self.get_numcands()
        self.cands.sort(cmp_freq)
        if verbosity >= 1:
            print "  Searching for dupes..."
        ii = 0
        # Find any match
        while ii < self.get_numcands():
            jj = ii + 1
            if jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                # Find others that match
                jj += 1
                while jj < self.get_numcands() and \
                        Num.fabs(self.cands[ii].r-self.cands[jj].r) < r_err:
                    jj += 1
                matches = self.cands[ii:jj]
	        print cmp_snr
                matches.sort(cmp_snr)
                bestindex = self.cands.index(matches[0])
                #sigmas = [c.sigma for c in matches]
                #bestindex = Num.argmax(sigmas)+ii
                # flag the duplicates
                bestcand = self.cands[bestindex]
                # Add other matching cands as hit of highest-sigma cand
                for matchind in reversed(range(ii, jj)):
                    if matchind == bestindex:
                        # The current candidate is the highest-sigma cand
                        # Don't remove it
                        continue
                    match = self.cands[matchind]
                    bestcand.add_as_hit(match)
                    match.note = "This candidate is a duplicate of %s:%d" % \
                                (bestcand.filename, bestcand.candnum)
                    self.mark_as_duplicate(matchind)
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (match.filename, match.candnum, matchind)
                        print "    %s" % match.note
                # If the best candidate isn't at the same freq
                # as ii, then it's possible even more hits should
                # be added. So we don't increment the index
                # (note that the best cand has moved into position ii).
            else:
                ii += 1 # No candidates to be added as hits, move on
        if verbosity >= 1:
            print "Found %d candidates.\n" % self.get_numcands()
        self.cands.sort(cmp_snr)


    def remove_DM_problems(self, numdms, dmlist, low_DM_cutoff, verbosity=1):
        """Remove the candidates where any of the following are true:
            1) The number of hits is < numdms
            2) The highest S/N candidate occurs below a DM of low_DM_cutoff
            3) The minimum difference in DM indices between the hits is > 1

            Inputs:
                numdms: The minimum number of hits for a good candidate.
                dmlist: List of DMs.
                low_DM_cutoff: The lowest DM possible for a good candidate.
                verbosity: Verbosity level. (Default: 1)

            Ouputs:
                None
        """
        # Create a dictionary where the key is the dmstr 
        # and the values are the index
        dmdict = {}
        dms = Num.unique([float(dm) for dm in dmlist])
        dmstrs = ['%.2f'%dm for dm in dms]
        dmdict = dict(zip(dmstrs, range(len(dms))))
        numremoved = 0
        num_toofew = 0
        num_toolow = 0
        num_gaps = 0
        self.cands.sort(cmp_snr)
        for ii in reversed(range(len(self.cands))):
            currcand = self.cands[ii]
            # Remove all the candidates without enough DM hits
            if len(currcand.hits) < numdms:
                numremoved += 1
                num_toofew += 1
                currcand.note = "Candidate has only %d DM hits. This is less " \
                                "than minimum for 'good' cands (%d hits)" % \
                                (len(currcand.hits), numdms)
                self.mark_as_bad(ii, 'dmproblem')
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            # Remove all the candidates where the max sigma DM is 
            # less than the cutoff DM
            # Recall - A hit is a 3-tuple: (DM, SNR, sigma)
            imax = Num.argmax(Num.array([hit[2] for hit in currcand.hits]))
            hitdm, hitsnr, hitsigma = currcand.hits[imax]
            if float(hitdm) <= low_DM_cutoff:
                numremoved += 1
                num_toolow += 1
                currcand.note = "Hit with max sigma (%g) has dm (%.2f) " \
                                "<= low DM cutoff (%.2f) " % \
                                    (hitsigma, hitdm, low_DM_cutoff)
                self.mark_as_bad(ii, 'dmproblem')
                if verbosity >= 2:
                    print "Removing %s:%d (index: %d)" % \
                            (currcand.filename, currcand.candnum, ii)
                    print "    %s" % currcand.note
                continue

            # Remove all the candidates where there are no hits at consecutive DMs
            if len(currcand.hits) > 1:
                currcand.hits.sort(cmp_dms)
                dm_indices = Num.asarray([dmdict["%.2f"%currcand.hits[jj][0]]
                                          for jj in range(len(currcand.hits))])
                min_dmind_diff = min(dm_indices[1:] - dm_indices[:-1])
                if min_dmind_diff > 1:
                    numremoved += 1
                    num_gaps += 1
                    currcand.note = "DM list of hits has gaps (i.e. " \
                                    "consecutive DMs don't have hits)."
                    self.mark_as_bad(ii, 'dmproblem')
                    if verbosity >= 2:
                        print "Removing %s:%d (index: %d)" % \
                                (currcand.filename, currcand.candnum, ii)
                        print "    %s" % currcand.note
                    continue

        if verbosity >= 1:
            print "Removed %d candidates with DM problems.\n" % numremoved
        if verbosity >= 2:
            print "  # with too few hits:", num_toofew
            print "  # with peak SNR too low:", num_toolow
            print "  # with gaps in DM hits:", num_gaps

    def print_cand_summary(self, summaryfilenm=None):
        """Write a summary of all candidates to file (or stdout).

            Input:
                summaryfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if summaryfilenm is None:
            summaryfile = sys.stdout
        elif summaryfilenm in [sys.stdout, sys.stderr]:
            summaryfile = summaryfilenm
        else:
            summaryfile = open(summaryfilenm, "w")
        summaryfile.write("   Candlist contains %d 'good' candidates\n" %len(self.cands))
        summaryfile.write("      # Known RFI rejects:           %d\n" %len(self.badlists['knownbirds']))
# REMOVED HARMONIC/TRESHOLD
        summaryfile.write("      # Duplicate candidates:        %d\n" %len(self.duplicates))
        summaryfile.write("      # Harmonic candidates:         %d\n" %len(self.badlists['harmonic']))
        summaryfile.write("      # Candidates with DM problems: %d\n" %len(self.badlists['dmproblem']))
        if summaryfilenm not in [None, sys.stdout, sys.stderr]:
            summaryfile.close()
  
    def write_cand_report(self, reportfilenm=None):
        """Write a report of all bad candidates to file (or stdout).

            Input:
                reportfilenm: Name of file to write to. If None write to stdout.
                    (Default: write to stdout).

            Outputs:
                None
        """
        if reportfilenm is None:
            reportfile = sys.stdout
        else:
            reportfile = open(reportfilenm, "w")
# CHANGED THE TITLES
        reportfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                       "SNR".center(8) +  "P(ms)".center(14) +
                       "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
        badcands = self.get_all_badcands()
        for badcand in badcands:
            reportfile.write("%s (%d)\n" % (str(badcand), len(badcand.hits)))
            reportfile.write("    Note: %s\n\n" % badcand.note)
        if reportfilenm is not None:
            reportfile.close()
        

    def __add__(self, other):
        copy_of_self = copy.deepcopy(self)
        copy_of_self.extend(other)
        return copy_of_self

    def get_numcands(self):
        """Get the number of good candidates (i.e. len(self.cands)).

            Inputs:
                None

            Outputs:
                None
        """
        return len(self)

    def __len__(self):
        # return the number of good candidates
        return len(self.cands)

    def extend(self, other):
        """Extend Candlist with another. This combines
            the candidates, as well as the lists of bad cands.
        
            Inputs:
                other: A second Candlist object to extend from.

            Outputs:
                None - the original Candlist object is extended in place.
        """
        self.cands.extend(other.cands)
        self.duplicates.extend(other.duplicates)
        for key in other.badlists:
            bad = self.badlists.setdefault(key, [])
            bad.extend(other.badlists[key])

    def to_file(self, candfilenm=None):
        """Write Candlist to file (or stdout).
            
            Input:
                candfilenm: Name of file to write to. If None,
                    write to stdout. (Default: write to stdout).

            Outputs:
                None
        """
        if candfilenm is None:
            candfile = sys.stdout
        else:
            candfile = open(candfilenm, "w")
        candfile.write("#" + "file:".center(20)+"candnum".center(30) + "			P(ms)".center(9) +
                       "SNR".center(12) + "DM".center(2) +
                       "dt (ms)".center(18) + "numhits" + "\n")
        for goodcand in self.cands:
            candfile.write("%s (%d)\n" % (str(goodcand), len(goodcand.hits)))
            if (len(goodcand.hits) > 1):
                goodcand.hits.sort(cmp_dms)
                for hit in goodcand.hits:
                    numstars = int(hit[2]/3.0)
                    candfile.write("  DM=%6.2f SNR=%5.2f   "%hit + \
                                    numstars*'*' + '\n')
        if candfilenm is not None:
            candfile.close()
