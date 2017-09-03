# AnalysisUtils

[![Build Status](https://travis-ci.org/ajkeller34/AnalysisUtils.jl.svg?branch=master)](https://travis-ci.org/ajkeller34/AnalysisUtils.jl)
[![Coverage Status](https://coveralls.io/repos/ajkeller34/AnalysisUtils.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ajkeller34/AnalysisUtils.jl?branch=master)
[![codecov.io](http://codecov.io/github/ajkeller34/AnalysisUtils.jl/coverage.svg?branch=master)](http://codecov.io/github/ajkeller34/AnalysisUtils.jl?branch=master)

Assorted analysis utilities.

Fitting an avoided crossing calculated using Sonnet:

```
using Touchstone, FileIO, AxisArrays, Images, AnalysisUtils
using AnalysisUtils: avoided_Cc

dset = loadset("BusResonator_17") # is a folder of .s2p files
data_pre = dset[:mag, :S,2,1, :, :]

# in the next line, `lentof` is some function to convert resonator length to frequency
data = AxisArray(20*log10.(data_pre), axes(data_pre)[1],
    Axis{:busf}( lentof.(data_pre[Axis{:BusLengthControl}].val) ))

v = trackextrema(indices(data, Axis{:busf}), findextrema(data, Axis{:f})...)

vm, vp = v[1], v[2] # these numbers may change depending on the data set
avoided_lsq(t) = sqrt(sum((avoided_Cc.(t[1], vp[2,:], t[2], +) .- vp[1,:]).^2) +
                      sum((avoided_Cc.(t[1], vm[2,:], t[2], -) .- vm[1,:]).^2))

fit = optimize(avoided_lsq, [10.4, 0.5], NelderMead())
b1, b2 = minimum(data[Axis{:busf}].val), maximum(data[Axis{:busf}].val)
let r = linspace(b1, b2, 100)
    plot!(r, ω.(fit.minimizer[1], r, fit.minimizer[2], +))
    plot!(r, ω.(fit.minimizer[1], r, fit.minimizer[2], -))
end
```
