#!/usr/bin/env python

"""Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


"""A little utility script to sort the merged output text files.

These files have tab-separated lines that start with the model name and then have
numeric fields, often in scientific notation.  The standard sort utility doesn't
recognise scientific notation and hence sorts these incorrectly.

Arguments are file names to sort; output is written to stdout.
"""

import itertools
import operator
import sys

# Utility methods
def make_key(line):
    """The sort key is defined as a tuple of the items on the line, converted to float where possible."""
    key = []
    for item in line.split():
        try:
            key.append(float(item))
        except:
            key.append(item)
    return tuple(key)

# Read the files
lines = []
for fname in sys.argv[1:]:
    lines.extend(open(fname).readlines())

# Define the sort key information
line_keys = map(make_key, lines)

# Zip lines with their keys
zipped = [(line,) + key for line,key in itertools.izip(lines, line_keys)]
# Secondary sort by timestep, largest first
zipped.sort(key=operator.itemgetter(4), reverse=True)
# Primary sort by model,solver,LT
zipped.sort(key=operator.itemgetter(1,2,3))

# Print output
for pair in zipped:
    print pair[0],
