#! /usr/bin/env python3
'''
MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses
'''

# standard imports
from datetime import datetime
import sys

# useful constants
VERSION = '0.0.1'
global QUIET; QUIET = False
global LOGFILE; LOGFILE = None

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to log (prefixed by current time)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    if not QUIET:
        print(tmp, file=sys.stderr, end=end); sys.stderr.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# main program execution
def main():
    print_log("=== MultiVirusConsensus (MVC) v%s ===" % VERSION)
    pass # TODO

# run tool
if __name__ == "__main__":
    main()
