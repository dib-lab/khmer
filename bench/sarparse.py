#! /usr/bin/env python
import datetime
from collections import namedtuple
import time
import gzip

def parse_pidstat(filename):
    """
    Parse the 'pidstat' output.

    12:08:07 AM   UID       PID  minflt/s  majflt/s     VSZ    RSS   %MEM  Command
    """
    parsed = []

    typ = namedtuple('PIDStats', 'time UID PID minflt_s majflt_s VSZ RSS '
                                 'Percent_MEM Command')

    for n, line in enumerate(gzip.open(filename, "rt")):
        if n < 3:                       # skip headers
            continue
        if line.startswith('Average'):  # skip footers
            continue

        line = line.strip().split()

        assert len(line) == 10, len(line)
        line2 = [line[0] + ' ' + line[1], line[2], line[3]]
        line2.extend(( float(x) for x in line[4:-1]))
        line2.extend([line[-1]])

        t = typ(*line2)
        parsed.append(t)

    return parsed

def parse_sartime(t):
    """
    convert a 'sar' timestamp into a hour/minute/second:

    06:21:38 PM
    """
    t = t.split(' ')[0]
    hour, minute, second = t.split(':')
    return int(hour), int(minute), int(second)

def get_sar_start_time(sar_data, timelog_timestamp):
    """
    Since 'sar' timestamp output is braindead and has neither date nor
    timezone, register it to our timelog timestamp.  Assume same
    day/hour for start. BORK BORK BORK @CTB.
    """
    sar_time = sar_data[0][0]

    hour, minute, second = parse_sartime(sar_time)

    d = datetime.datetime(timelog_timestamp.year,
                          timelog_timestamp.month,
                          timelog_timestamp.day,
                          timelog_timestamp.hour,
                          int(minute),
                          int(second))

    return d

def make_timediff(sar_data):
    """
    Calculate 'sar' sampling frequency in seconds. Must be less than 1 hr.

    Assume the first two times are immediately adjacent (don't use disk
    output!)
    """
    t1 = parse_sartime(sar_data[0][0])
    t2 = parse_sartime(sar_data[1][0])

    assert t1[0] == t2[0]
    secdiff = (t2[1] - t1[1]) * 60 + t2[2] - t1[2]

    return secdiff

def fixtime(sar_data, start, secdiff):
    "Fix the hh::mm::ss timestamps output by 'sar' to full datetimes."
    delta = datetime.timedelta(0, secdiff)

    currentime = start

    sar_data2 = []
    for x in sar_data:
        sar_data2.append(x._replace(time=currentime))
        currentime += delta

    return sar_data2

def make_time(x, start=None):
    "Convert datetimes into seconds with time.mktime, optionally - start."
    sub = 0
    if start:
        sub = time.mktime(start.timetuple())
    return time.mktime(x.timetuple()) - sub
