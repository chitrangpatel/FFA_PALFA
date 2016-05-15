import ffa_tools as ft
import ffa_stages as fs
import sifting_ffa as sf
import sys

"""
Run this script once you have applied the FFA on all de-dispersed time series.
The candidate file inputed should include candidates from multiple DMs (if it applies)
The argument filenm is the _cands.ffa file for a particular beam.
"""

filenm = sys.argv[1]
sf.siftDM(filenm)
ft.format_final_cands(filenm)
print "Done with ",filenm
