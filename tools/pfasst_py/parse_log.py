# coding=utf-8
"""
.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""
import datetime
import os
import re
import pickle
import pathlib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class DataPattern(object):
    """
    regex pattern provider for log data

    RE will get applied to each log line
    """

    RE = None


class AdvecDiffDataPattern(DataPattern):
    RE = re.compile(r'.*MPI\s+(?P<rank>[0-9]+).*step:\s+(?P<step>[0-9]+)\s+iter:\s+(?P<iter>[0-9]+)\s+n1:\s+(?P<n1>[0-9]+)\s+n2:\s+(?P<n2>[0-9]+)\s+residual:\s+(?P<residual>[0-9\.e+-]+)\s+err:\s+(?P<error>[0-9\.e+-]+).*')



class BorisDataPattern(DataPattern):
    RE = re.compile(r'.*MPI\s+(?P<rank>[0-9]+).*step:\s+(?P<step>[0-9]+)\s+iter:\s+(?P<iter>[0-9]+)\s\((?P<level>[a-z]+)\s*\)\s+residual:\s+(?P<residual>[0-9\.e+-]+)\s+energy\sdrift:\s+(?P<drift>[0-9\.e+-]+)\s+total\senergy:\s+(?P<energy>[0-9\.e+-]+).*')


class LogParser(object):
    TIME_MATCH = re.compile(r'^(?P<time_h>[0-9]+):(?P<time_m>[0-9]+):(?P<time_s>[0-9]+),(?P<time_ms>[0-9]+)\s.*')

    def __init__(self, file, data_pattern=None):
        assert(issubclass(data_pattern, DataPattern))
        assert(os.access(str(file), os.R_OK))
        self._file = file
        self._file_date = self._get_file_date()
        self._data_pattern = data_pattern
        self._first_timestamp = None
        self._parsed = False
        self._dataframe = None

    def get_dataframe(self):
        """
        gets data as a pandas.DataFrame
        """
        if not self._parsed:
            self._parse()
        return self._dataframe

    def _parse(self):
        self._first_timestamp = None
        data_dicts = list()
        with open(str(self._file), mode='r') as f:
            for line in f:
                dict = self._log_line_to_dict(line)
                if dict is not None:
                    timestamp = self._get_line_timestamp(line)
                    if self._first_timestamp is None:
                        self._first_timestamp = timestamp
                    dict['time'] = (timestamp - self._first_timestamp)
                    data_dicts.append(dict)
        self._dataframe = pd.DataFrame(data_dicts).convert_objects(convert_numeric=True)
        self._parsed = True

    def _log_line_to_dict(self, line):
        """
        extracts data from file based on data_pattern as a dict

        :param line:
        :return:
        """
        assert(isinstance(line, str))
        data_match = self._data_pattern.RE.match(line)
        if data_match:
            return data_match.groupdict()
        else:
            return None

    def _get_line_timestamp(self, line):
        """
        extracts timestamp from given log line

        :param line:
        :return:
        """
        assert(isinstance(line, str))
        assert(self._file_date is not None)

        time_matched = LogParser.TIME_MATCH.match(line)
        if time_matched:
            g = time_matched.groupdict()
            line_time = datetime.time(int(g['time_h']), int(g['time_m']), int(g['time_s']), int(g['time_ms']))
            return datetime.datetime.combine(self._file_date, line_time)
        else:
            return None

    def _get_file_date(self):
        """
        gets modification date from given file

        :return:
        """
        return datetime.datetime.fromtimestamp(os.stat(self._file).st_mtime).date()

    def __str__(self):
        _str  = "file: %s\n" % self._file
        _str += "  created: %s\n" % self._file_date
        _str += "  parsed: %s\n" % self._parsed
        if self._parsed:
            _str += "    entries: %d\n" % len(self._dataframe)
            _str += "    duration: %s\n" % self._dataframe.tail(1)['time'].iget(0)
        return _str


class LogParserSet(object):
    """
    gathers and parses rank-specific logs of a single run
    """

    def __init__(self, dir, log_prefix, max_ranks, pattern):
        """
        prepares the run parser

        :param dir: the path where all the logs are
        :param log_prefix: the log file prefix as specified with `--log_prefix`
        :param max_ranks: how many ranks were used
        :param pattern: regex pattern provider
        """
        self._dir = dir
        self._log_prefix = log_prefix
        self._max_ranks = max_ranks
        self._parser = list()
        self._pattern = pattern
        self._all_data = None

    def load_and_parse(self):
        """
        greates a LogParser for each log file and parses it
        """
        for rank in range(0, self._max_ranks):
            file = self._log_prefix + "%04d" % rank + ".log"
            self._parser.append(LogParser(self._dir + "/" + file, self._pattern))

    def get_dataframe(self):
        """
        a single pandas.DataFrame with data from all logs
        """
        self._concat_data()
        return self._all_data

    def get_duration(self):
        """
        gets last timestamp reading from last processor

        this assumes that the rank with the highest rank number is the one who terminates last
        """
        self._concat_data()
        return self._all_data[self._all_data['rank'] == self._all_data['rank'].max()].time.max()

    def _concat_data(self):
        if self._all_data is None:
            self._all_data = pd.concat([parser.get_dataframe() for parser in self._parser], ignore_index=True)


class ScaleStudy(object):
    """
    utilities for a scaling study
    """
    RE_RUNJOB_NP = re.compile(r'^runjob.*--np\s(?P<np>[0-9]+).*')
    RE_RUNJOB_NPN = re.compile(r'^runjob.*--ranks-per-node\s(?P<npn>[0-9]+).*')

    def __init__(self, path, pattern, attach=None):
        self._path = pathlib.Path(path).resolve()
        assert(self._path.exists())
        self._pattern = pattern
        self._attach = attach
        self._runs = list()

    def gather(self):
        print("gathering in base path '%s'" % self._path)
        for run_dir in self._path.iterdir():
            if run_dir.is_dir():
                print("  - %s" % run_dir.name)
                run_log = list(run_dir.glob('run.log'))[0].resolve()
                if isinstance(run_log, pathlib.Path) and run_log.is_file():
                    d = {
                        'dir': run_dir.name
                    }
                    if self._attach is not None:
                        d.update(self._attach)
                        for k, v in d.items():
                            if k != 'dir':
                                print("      %10s: %s" % (k, v))
                    with open(run_log.__str__(), mode='r') as f:
                        for line in f:
                            np_match = self.RE_RUNJOB_NP.match(line)
                            if np_match is not None:
                                d['np'] = int(np_match.groupdict()['np'])
                                print("      np        : %d" % d['np'])
                            npn_match = self.RE_RUNJOB_NPN.match(line)
                            if npn_match is not None:
                                d['npn'] = int(npn_match.groupdict()['npn'])
                                print("      npn       : %d" % d['npn'])
                    if d['np'] is not None and isinstance(d['np'], int):
                        d['parser_set'] = LogParserSet(run_dir.__str__(), 'mpi-rank-', d['np'], self._pattern)
                        d['parser_set'].load_and_parse()
                        d['duration'] = d['parser_set'].get_duration()
                        print("      duration  : %s" % d['duration'])
                    else:
                        print("    WARNING: could not parse logs")
                    self._runs.append(d)
        return self

    def as_dataframe(self):
        return pd.DataFrame(self._runs).drop('parser_set', 1)

    def plot_np_vs_duration(self):
        df = self.as_dataframe().drop('dir', 1).sort('np')
        print(df)
        nruns = tuple(df['np'].as_matrix().tolist())
        ind = np.arange(len(nruns))
        fig, ax = plt.subplots()
        ax.plot(ind, df['duration'])
        ax.set_title('num particles %d' % df['nparticles'].max())
        ax.set_ylabel('duration (ns)')
        ax.set_yscale('log')
        ax.set_xticks(ind)
        ax.set_xlabel('num procs')
        ax.set_xticklabels(nruns)

    def pickle(self, name=None):
        if name is None:
            name = self._path.resolve().stem
        with open(str(self._path / name) + '.pickle', mode='wb') as f:
            print("pickling into: '%s'" % f.name)
            pickle.dump(self, f)

    @classmethod
    def from_pickle(self, file):
        obj = None
        with open(file, mode='rb') as f:
            print("load pickle from: '%s'" % f.name)
            obj = pickle.load(f)
        return obj
