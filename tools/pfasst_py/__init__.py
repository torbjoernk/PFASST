# coding=utf-8
"""
.. moduleauthor:: Torbjörn Klatt <t.klatt@fz-juelich.de>
"""
import os

PFASST_BASE_DIR = os.path.abspath(os.getcwd() + "../../../")
PFASST_RUN_DIR = os.path.abspath(PFASST_BASE_DIR + "/run")

from .parse_log import *
