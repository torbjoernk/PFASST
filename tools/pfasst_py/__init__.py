# coding=utf-8
"""
.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""
import os
import pathlib

PFASST_BASE_DIR = str(pathlib.Path(os.path.realpath(__file__)).resolve().parent.parent.parent)
PFASST_RUN_DIR = PFASST_BASE_DIR + '/run'

from .parse_log import *
