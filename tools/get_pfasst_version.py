# encoding: utf-8

from __future__ import print_function

from sys import version_info

if version_info[0] == 2 and version_info[1] < 7:
    raise SystemExit("Insufficient Python Interpreter Version (%s < 2.7)" % (version_info,))

import re
import subprocess

from os.path import dirname, abspath, join

# get git version
try:
    version = subprocess.check_output(['git', 'describe', '--dirty']).decode().strip()
except subprocess.CalledProcessError:
    # `git describe` will fail in case the Git tags are not fetched
    #  this happens when cloning the repo only from a fork pruned of tags
    print("Could not determine PFASST version. "
          "Inspect the output of 'git describe --dirty' for errors.")
    version = 'unknown'

# read in site_config.hpp
base = dirname(dirname(abspath(__file__)))
site_config = abspath(join(base, 'include', 'pfasst', 'site_config.hpp'))

with open(site_config, 'r') as f:
    config = f.read()

# reset version
new_config = re.sub('VERSION = "[^"]*"', 'VERSION = "{version}"'.format(version=version), config)

# write site_config if it has changed
if config != new_config:
    with open(site_config, 'w') as f:
        f.write(new_config)
    print("PFASST++ version set to: %s" % (version,))
else:
    print("PFASST++ version did not change (still %s)" % (version,))
