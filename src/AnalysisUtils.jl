module AnalysisUtils
using AxisArrays, Images

"""
    findextrema(data::AxisArray, axis::Type{<:Axis}, maxima=true; smoother=Kernel.gaussian((3,1), (5,1)))
Look for extrema in a 2D matrix.

- `data`: Array containing the data.
- `axis`: The axis along which to look for extema.
- `maxima`: If true, find maxima, else find minima.
- `smoother`: A kernel for doing smoothing prior to extrema identification.

Returns `(idxs, vals)` where:

- `idxs`: a `2 x N` matrix of indices where the `N` extrema are found in a 2D matrix. The
  second row corresponds to the control parameter.
- `vals`: a `2 x N` matrix of axis values where the `N` extrema are found in a 2D matrix.
  The second row corresponds to the control parameter.
"""
function findextrema(data::AxisArray, axis::Type{<:Axis}, maxima=true; smoother=Kernel.gaussian((3,1), (5,1)))
    # Do some smoothing using Images.jl just to avoid worrying about noise
    smdata = imfilter(data, smoother)

    # Peak finding along axis 1 using Images.jl.
    # Ignore maxima at boundaries for probe frequency axis, but the boundaries along the bus
    # frequency axis are fine. findlocalmaxima returns `CartesianIndex` objects.
    f = ifelse(maxima, findlocalmaxima, findlocalminima)
    cartesianindices = f(smdata, [axisdim(data, axis)], (x->!isa(x,axis)).(axes(data)))

    # Now we want to group the local maxima we found.
    # idxs: a 2xN matrix of (probe, bus) indexes
    # vals: a 2xN matrix of (probe, bus) frequencies
    idxs = reshape(reinterpret(Int, cartesianindices), 2, :)
    vals = getindex.(axes(data), idxs)

    return (idxs, vals)
end

"""
    trackextrema(idxrange, idxs, vals, follow_trajectory = true)

Often the result of a measurement or simulation will yield a signal that has some
number of peaks (or generically, extrema) that move around in a continuous way as a
function of a control parameter. The peaks may change in number as a function of the control
parameter: peaks can vanish or appear spontaneously.

In either measurement or simulation, the control parameter is changed discretely, leading
to a set of extrema for each value of the control parameter. `trackextrema` provides an
automated way to cluster these discrete observations into groups across all values of the
control parameter.

This seems to perform better for this purpose than clustering algorithms like k-means,
hierarchical clustering, etc., which do not exploit our key assumptions: that samples being
clustered are low noise, discretized samples of continuous lines.

- `idxrange`: the range of indices on the control parameter axis to use for tracking.
- `idxs`: a `2 x N` matrix of indices where the `N` extrema are found in a 2D matrix. The
  second row corresponds to the control parameter.
- `vals`: a `2 x N` matrix of axis values where the `N` extrema are found in a 2D matrix.
  The second row corresponds to the control parameter.
- `follow_trajectory`: optionally use the peak change in peak position to estimate where
  the next peak position should be (linear interpolation). This estimate is then used in
  the distance matrix calculation, instead of the previous peak position, to decide which
  cluster the next peak should belong to. This is useful when you have intersecting lines.
"""
function trackextrema(idxrange, idxs, vals, follow_trajectory = true)
    state = start(idxrange)
    prevvals = vals[:, findin(idxs[2,:], Base.first(idxrange))]
    nprevvals = size(prevvals,2)
    valdict = Dict{Int, Matrix{Float64}}()  # keys are "cluster indexes"
    first = true
    while !done(idxrange, state)
        previdx, state = next(idxrange, state)
        previdx == last(idxrange) && break
        nextidx = previdx+1

        # These are the peaks in the next column to process
        nextvals = vals[:, findin(idxs[2,:], nextidx)]
        nnextvals = size(nextvals,2)

        # Values to use in distance matrix calculation..
        prevdist = prevvals
        if !first && follow_trajectory
            prevdist = copy(prevvals)
            dif = nextvals[2,1] - prevvals[2,1]
            for i in 1:npairs
                prevdist[1:2, deltaidxs[i]] += deltas[:,i] * dif/deltas[2,i]
            end
        end

        # Create a distance matrix. It will always have # cols >= # rows,
        # the iswap function lets us be sure of that and handle the indices easily.
        iswap = nprevvals <= nnextvals ? (x,y)->(x,y) : (x,y)->(y,x)
        distmat = Matrix{Float64}(iswap(nprevvals, nnextvals)...)

        # We will attempt to pair up peaks from one column to the next.
        npairs  = min(size(distmat)...)
        ngroups = max(size(distmat)...)

        # Calculate a distance matrix
        for i in 1:nprevvals, j in 1:nnextvals
            distmat[iswap(i,j)...] = norm(-(iswap(prevdist[1:2,i] , nextvals[1:2,j])...))
        end

        # Pairs will be 2 x npairs, containing the indices from prevvals and nextvals.
        pairs = Matrix{Int}(2, npairs)
        for i in 1:npairs
            pairs[:, i] = [iswap(i, findmin(distmat[i,:])[2])...]
        end

        # First time through the loop, valdict is empty, we can just store these using the
        # pair index as the cluster index.
        nextvalclusters = Vector{Int}(nnextvals)
        if first
            first = false
            for i in 1:npairs
                valdict[i] = [get(valdict, i, Matrix{Float64}(2,0)) prevvals[1:2, pairs[1,i]]]
                valdict[i] = [get(valdict, i, Matrix{Float64}(2,0)) nextvals[:, pairs[2,i]]]
                nextvalclusters[pairs[2,i]] = i
            end
            nextcidx = npairs + 1
        else
            for i in 1:npairs
                cidx = prevvals[3, pairs[1,i]]  # "cluster index"
                existing = get(valdict, cidx, Matrix{Float64}(2,0))
                valdict[cidx] = [existing nextvals[:, pairs[2,i]]]
                nextvalclusters[pairs[2,i]] = cidx
            end
        end

        if nprevvals < nnextvals
            unpaired = setdiff(1:nnextvals, pairs[2,:])
            for i in unpaired
                existing = get(valdict, nextcidx, Matrix{Float64}(2,0))
                valdict[nextcidx] = [existing nextvals[:, i]]
                nextvalclusters[i] = nextcidx
                nextcidx += 1
            end
        end

        if follow_trajectory
            deltas = Matrix{Float64}(2, npairs)
            deltaidxs = Vector{Int}(npairs)
            for i in 1:npairs
                deltas[:, i] = nextvals[:, pairs[2,i]] - prevvals[1:2, pairs[1,i]]
                deltaidxs[i] = pairs[2,i]
            end
        end

        # prepare for next pass
        previdx = nextidx
        nprevvals = nnextvals

        # retain the cluster index for the next round
        prevvals = [nextvals; reshape(nextvalclusters, 1, nnextvals)]
    end

    return valdict
end

end # module
