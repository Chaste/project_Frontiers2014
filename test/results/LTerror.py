#!/usr/bin/env python
"""Copyright (c) 2005-2014, University of Oxford.
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

"""
This script analyses the error introduced by lookup tables for a range of cell models and solvers.
It takes the standard output log file produced by TestOdeSolvingTimesLiteratePaper.hpp and looks
for lines of the form:

Model <mname> solver '<sname>' [and lookup tables] MRMS error = <number>

to see how the reported error varies with lookup tables on/off.

We also check for lines of the form:

Simulating <mname> with solver '<sname>' [and lookup tables] using timestep <number>

to check whether the same timesteps were used with/without lookup tables, and lines of the form

Model <mname> solver '<sname>' [and lookup tables] took time <number>s per simulated sec

to determine the speedup from lookup tables.

CSV data is written to standard output summarising the results.
"""

import re
import sys

if len(sys.argv) < 2:
    print >> sys.stderr, "No log file name supplied!"

log = open(sys.argv[1], 'r')

# Regular expressions matching lines of interest
re_base = r" (?P<model>\S+) (with )?solver '(?P<solver>[^']+)' (?P<lt>and lookup tables )?"
error_re = re.compile("Model" + re_base + r"MRMS error = (?P<error>[0-9e.+-]+)")
timestep_re = re.compile("Simulating" + re_base + r"using timestep (?P<dt>[0-9e.+-]+)")
time_re = re.compile(r"Model" + re_base + r"took time (?P<time>[0-9e.+-]+)s per simulated sec")
count_re = re.compile(r"Model (?P<model>\S+) solver '(?P<solver>[^']+)' has tables.+?; total (?P<num>\d+) with (?P<keys>\d+) keys.")

errors = {}
timesteps = {}
times = {}
lt_counts, lt_keys = {}, {}
combos = set()
error_incrs = []
speedups = []

for line in log:
    m = error_re.match(line)
    if m:
        key = (m.group('model'), m.group('solver'), bool(m.group('lt')))
        errors[key] = m.group('error')
        combos.add((m.group('model'), m.group('solver')))
        continue
    m = timestep_re.match(line)
    if m:
        key = (m.group('model'), m.group('solver'), bool(m.group('lt')))
        timesteps[key] = m.group('dt')
        combos.add((m.group('model'), m.group('solver')))
        continue
    m = time_re.match(line)
    if m:
        key = (m.group('model'), m.group('solver'), bool(m.group('lt')))
        times[key] = m.group('time')
        combos.add((m.group('model'), m.group('solver')))
        continue
    m = count_re.match(line)
    if m:
        key = (m.group('model'), m.group('solver'))
        lt_counts[key] = m.group('num')
        lt_keys[key] = m.group('keys')
        combos.add((m.group('model'), m.group('solver')))

print >> sys.stderr, "Found", len(errors), "error lines,", len(times), "timing lines and", len(timesteps), "timestep lines"

print 'Model,Solver,Num LTs,Num LT keys,Normal Error,LT Error,Normal dt,LT dt,Timesteps Match,Error Increase (%),Speedup'
for model, solver in sorted(combos):
    line = model + ',' + solver
    if (model, solver) in lt_counts:
        line += ',%s,%s' % (lt_counts[(model, solver)], lt_keys[(model, solver)])
    else:
        line += ',,'
    err_with = errors.get((model, solver, True), '')
    err_without = errors.get((model, solver, False), '')
    line += ',' + err_without + ',' + err_with
    for lt in [False, True]:
        line += ',' + timesteps.get((model, solver, lt), '')
    if err_with and err_without:
        dt_match = timesteps.get((model, solver, False)) == timesteps.get((model, solver, True))
        error_incr = (float(err_with) - float(err_without)) / float(err_without) * 100.0
        error_incrs.append(error_incr)
        line += ',%s,%.2f' % (dt_match, error_incr)
    else:
        line += ',,'
    time_with, time_without = times.get((model, solver, True)), times.get((model, solver, False))
    if time_with and time_without:
        speedup = float(time_without) / float(time_with)
        speedups.append(speedup)
        line += ',%.2f' % speedup
    else:
        line += ','
    print line

def summarise(vector, metric):
    """Summarise a metric. The vector gives the metric value for each model/solver combination."""
    mean = sum(vector) / len(vector)
    var = sum((mean - value) ** 2 for value in vector) / (len(vector)-1)
    print >> sys.stderr, "Summary of %s: min=%.2f, max=%.2f, mean=%.2f, var=%.2f" % (metric, min(vector), max(vector), mean, var)

# Summarise the speedups
summarise(speedups, "LT speedup")

# Summarise how the errors change overall
summarise(error_incrs, "error changes")
print >> sys.stderr, "Errors over 5% MRMS threshold:"
for model, solver in sorted(combos):
    if float(errors.get((model, solver, True), 0)) > 0.05:
        print >> sys.stderr, "   ", model, solver
print >> sys.stderr, "Combinations where LT caused simulation failure:"
for model, solver in sorted(combos):
    if not (model, solver, True) in times:
        print >> sys.stderr, "   ", model, solver
