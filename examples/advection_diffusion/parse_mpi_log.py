# coding=utf-8
"""
.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""
from __future__ import print_function

import os
import datetime
import pandas as pd
import numpy as np
import re


def parse_log(file):
    TIME_MATCH = re.compile(r'^(?P<time_h>[0-9]+):(?P<time_m>[0-9]+):(?P<time_s>[0-9]+),(?P<time_ms>[0-9]+)\s.*')
    LINE_MATCH = re.compile(r'.*\sstep:\s+(?P<step>[0-9]+)\s+iter:\s+(?P<iter>[0-9]+)\s+n1:\s+(?P<n1>[0-9]+)\s+n2:\s+(?P<n2>[0-9]+)\s+residual:\s+(?P<residual>[0-9\.e-]+)\s+err:\s+(?P<error>[0-9\.e-]+).*')
    data_type = np.dtype([('step', np.uint32), ('iter', np.uint32), ('n1', np.uint32), ('n2', np.uint32),
                          ('residual', np.float64), ('error', np.float64), ('time', np.float64)])

    file_date = datetime.datetime.fromtimestamp(os.stat(file).st_mtime).date()
    print("file:      %s" % file)
    print("  created: %s" % file_date)
    matched = []
    first_datetime = None
    with open(file, mode='r') as f:
        for line in f:
            data_match = LINE_MATCH.match(line)
            if data_match:
                time_match = TIME_MATCH.match(line)
                line_datetime = None
                if time_match:
                    g = time_match.groupdict()
                    line_time = datetime.time(int(g['time_h']), int(g['time_m']), int(g['time_s']), int(g['time_ms']))
                    line_datetime = datetime.datetime.combine(file_date, line_time)
                    if not first_datetime:
                        first_datetime = line_datetime
                g = data_match.groupdict()
                time_diff = line_datetime - first_datetime
                matched.append(np.array((int(g['step']), int(g['iter']), int(g['n1']), int(g['n2']), float(g['residual']), float(g['error']), time_diff.total_seconds()), dtype=data_type))
    print("%d records" % len(matched))
    matched = np.asarray(matched)
    return pd.DataFrame({
        'time': matched['time'],
        'step': matched['step'],
        'iter': matched['iter'],
        'n1': matched['n1'],
        'n2': matched['n2'],
        'residual': matched['residual'],
        'error': matched['error']
    })

