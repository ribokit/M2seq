"""
utils.py

All utility functions pertaining to the main m2seq.py process are stored here
"""

import os
import time

# https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def timeStamp():
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return '%i:%i:%i, %i/%i'%(hour, minute, second, month, day)


def make_dir(path):
    try:
        os.mkdir(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def get_sequence(sequencefile):
    sequencefile_lines = sequencefile.readlines()

    sequence = sequencefile_lines[1].strip().upper()
    # check for Us
    for e in sequence:
        if e == 'U':
            raise ValueError("your sequence in " + sequencefile + " has Us not Ts please fix!" )

    return sequence