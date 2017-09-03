# AnalysisUtils

[![Build Status](https://travis-ci.org/ajkeller34/AnalysisUtils.jl.svg?branch=master)](https://travis-ci.org/ajkeller34/AnalysisUtils.jl)
[![Coverage Status](https://coveralls.io/repos/ajkeller34/AnalysisUtils.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ajkeller34/AnalysisUtils.jl?branch=master)
[![codecov.io](http://codecov.io/github/ajkeller34/AnalysisUtils.jl/coverage.svg?branch=master)](http://codecov.io/github/ajkeller34/AnalysisUtils.jl?branch=master)

Assorted analysis utilities.

Usage suggestion:

```
using Touchstone, FileIO, AxisArrays, Images

dset = loadset("BusResonator_17")
data_pre = dset[:mag, :S,2,1, :, :]
data = AxisArray(20*log10.(data_pre), axes(data_pre)[1],
    Axis{:busf}( lentof.(data_pre[Axis{:BusLengthControl}].val) ))

trackextrema(indices(data, Axis{:busf}), findextrema(data, Axis{:f})...)
```
