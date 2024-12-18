using ImageMorphology
using Statistics

# Tabulated T₁ and T₂ values for the gel tubes at 3 Tesla at different temperatures (Kelvin)
const groundtruth_3T = (
    T₁ = Dict(
        #              1   2   3   4   5   6   7   8   9   10  11   12   13   14   15   16   17   18
        292 => Float64[200,299,296,463,450,444,604,596,754,745, 903,1448, 966,1034,1160,1276,1262,1415],
        296 => Float64[223,334,331,516,502,496,674,666,841,831,1007,1615,1078,1153,1293,1422,1407,1576],
        300 => Float64[249,372,368,574,559,552,750,741,935,924,1120,1796,1199,1282,1437,1581,1563,1750]
    ),
    T₂ = Dict(
        292 => Float64[52 ,73 ,113,53 ,94 ,154,95 ,136,116,157,137 ,390 ,224 ,167 ,214 ,204 ,184 ,174 ],
        296 => Float64[50 ,70 ,111,49 ,89 ,151,89 ,129,109,149,128 ,373 ,212 ,156 ,201 ,190 ,171 ,161 ],
        300 => Float64[48 ,67 ,110,46 ,85 ,148,84 ,124,103,142,121 ,359 ,203 ,148 ,191 ,180 ,161 ,151 ]
        )
)

groundtruth_3T.T₁[294] = (groundtruth_3T.T₁[292] + groundtruth_3T.T₁[296])/2
groundtruth_3T.T₂[294] = (groundtruth_3T.T₂[292] + groundtruth_3T.T₂[296])/2

"""
    function add_mean_and_std_per_tube(output::NamedTuple; tubesize = 65)

Given an MR-STAT reconstruction of gel tube phantoms (output of MR-STAT reconstruction is stored as a NamedTuple),
this function:
    1. Automatically tries to detect the several different gel tubes using the function _mask_per_tube;
    2. Computes mean values and standard deviations per tube;
    3. Adds the mean values and standard deviations to the NamedTuple (example: output.T₁.mean)
"""
function add_mean_and_std_per_region(output::NamedTuple, regions; plot_regions=true)

    # Mask maps based on MR-STAT proton density
    ρ = abs.(output.qmaps[1].ρ);

    # Check mask per tube
    if plot_regions

        figure()

        subplot(121)
            # overlay of ROIs on proton density map
            ρ_max = maximum(ρ);
            for i in eachindex(regions)
                ρ[regions[i]] .= 2*ρ_max;
            end
            imshow(ρ)
            title("Mask per tube")

        subplot(122)
            # numbered
            x = zeros(size(ρ)...)
            for (i,r) in enumerate(regions)
                x[r] .= i
            end
            imshow(x)
            colorbar()
    end

    # Compute mean and std values for each iteration of each reconstruction
    max_iter = length(output.qmaps)

    T₁_stats = (
        mean = Int.(round.(1000*[mean(output.qmaps[iter].T₁[region]) for region ∈ regions, iter ∈ 1:max_iter])),
        std  = Int.(round.(1000*[std(output.qmaps[iter].T₁[region]) for region ∈ regions, iter ∈ 1:max_iter]))
    )

    T₂_stats = (
        mean = Int.(round.(1000*[mean(output.qmaps[iter].T₂[region]) for region ∈ regions, iter ∈ 1:max_iter])),
        std  = Int.(round.(1000*[std(output.qmaps[iter].T₂[region]) for region ∈ regions, iter ∈ 1:max_iter]))
    )

    ρ_stats = (
        mean = Complex{Int}.(round.(1000*[mean(output.qmaps[iter].ρ[region]) for region ∈ regions, iter ∈ 1:max_iter])),
        std  = Int.(round.(1000*[std(output.qmaps[iter].ρ[region]) for region ∈ regions, iter ∈ 1:max_iter]))
    )

    ## sort based on T1 value of last iteration
    permutation = sortperm(T₁_stats.mean[:,end])

    T₁_stats.mean[:,:] = T₁_stats.mean[permutation,:]
    T₁_stats.std[:,:]  = T₁_stats.std[permutation,:]

    T₂_stats.mean[:,:] = T₂_stats.mean[permutation,:]
    T₂_stats.std[:,:]  = T₂_stats.std[permutation,:]

    ρ_stats.mean[:,:] = ρ_stats.mean[permutation,:]
    ρ_stats.std[:,:]  = ρ_stats.std[permutation,:]

    regions = regions[permutation]

    output = merge(output, [:T₁ => T₁_stats, :T₂ => T₂_stats, :ρ => ρ_stats])

    return output, regions
end

function add_mean_and_std_per_region(output::NamedTuple; tubesize=65,  plot_regions=true)

    ρ = abs.(output.qmaps[1].ρ);
    mask = ρ .> 0.05 * maximum(ρ);

    # Generate mask per tube
    regions = _mask_per_tube(mask, tubesize);

    output = add_mean_and_std_per_region(output, regions[1:end], plot_regions=plot_regions)

    return output
end

"""
    _mask_per_tube(mask::BitArray{2}, radius)

Automatically create a mask per individual tube by looking for connecting components,
computing their centers, and then selecting points within a certain distance from each center.
"""
function _mask_per_tube(mask::BitArray{2}, radius; background_water = false)

    # First, use ImageMorphology.label_components to label connected components of mask
    labeled_mask = ImageMorphology.label_components(mask)
    # Make a mask for each labeled component
    ntubes = maximum(labeled_mask)
    start_idx = background_water ? 2 : 1;
    tubes = [findall(x -> x == i, labeled_mask) for i = start_idx:ntubes]

    # Note that these masks could be used to compute mean values per tube,
    # but the masks are likely too large and include annoying boundaries
    # Therefore, we extract the centers and make masks with smaller radius

    # Extract tube centers
    x = tube -> round(Int, median(p[1] for p in tube))
    y = tube -> round(Int, median(p[2] for p in tube))

    tube_centers = [CartesianIndex(x(tube), y(tube)) for tube in tubes]

    # Distance function on CartesianIndices
    d(x::CartesianIndex, y::CartesianIndex) = √( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    d(x::Array{<:CartesianIndex}, y::CartesianIndex) = d.(x, (y,))

    # tube = Vector{CartesianIndices{2}}(undef, ntubes)
    allidx = CartesianIndices(mask)
    tubes = [allidx[ d.((center,), allidx) .< √radius ] for center ∈ tube_centers];

    return tubes
end


"""
    plot_tube_recons(reconstructions, iterations, descriptions; figtitle="MR-STAT", black_theme=false)

Given a Vector of MR-STAT reconstruction outputs (to which add_mean_and_std_per_tube has been applied),
plot the parameter maps as well as bar plots of mean and standard deviations per tube. Iterations is a vector with iteration numbers for each of the reconstructions. Description is a Vector of strings, one for each reconstruction.
"""
function plot_tube_recons(reconstructions, iterations, descriptions; figtitle="MR-STAT", black_theme=false, reference=false)

    # select iteration for each of the different reconstructions
    nrecons = length(reconstructions)

    has_reference = reference !== false

    # for barplots
    ntubes = size(reconstructions[1].T₁.mean,1)
    width = 0.15
    nrecons = length(iterations)
    totalwidth = (nrecons + 1 + has_reference) * width
    barcenters = LinRange(-totalwidth/2 + width/2, totalwidth/2 - width/2, nrecons + 1 + has_reference)
    # barcolors = ["#D36060" ,"#669BBC", "#264653", "#FFD275", "#2A9D8F", "#E85D75", "#E1CE7A"]
    # barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E1CE7A"]

    temp = 294
    # tube_idx = [1, 8,3,14,15,6,10,7,12,13,11,16,17]
    # tube_idx = [1,11,7,17,8,15,6,13,10,12,16,14,3]
    tube_idx = [11,13,7,9,17,12,8,18,16,13,6,3]

    # choose clims and cmaps for each reconstruction parameter
    clims = Dict(:T₁ => (0.0, 2.0), :T₂ => (0.0, 0.5), :ρ => (0.0, 2.0))
    cmaps = Dict(:T₁ => "inferno", :T₂ => "parula",   :ρ => "gray")

    if black_theme
        clrbg  = "black"
        clrtxt = "white"
        clrerr = "lightgray"
        barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E1CE7A"]
    else
        clrbg = "white"
        clrtxt = "black"
        clrerr = "black"
        barcolors = ["#E2798E", "#669BBC", "#861D32", "#264653", "#669BBC","#2A9D8F", "#E85D75", "#E1CE7A"]
    end

    fig = figure(facecolor=clrbg)

    parcnt = 1
    ncols = nrecons+1

    # make plots, one figure for each reconstruction parameter
    for par in [:T₁, :T₂]

        if parcnt == 1
            qmaps = [reconstructions[i].qmaps[iterations[i]].T₁     for i in 1:nrecons]
            means = [reconstructions[i].T₁.mean[:,iterations[i]]    for i in 1:nrecons]
            stds =  [reconstructions[i].T₁.std[:,iterations[i]]     for i in 1:nrecons]
        elseif parcnt == 2
            qmaps = [reconstructions[i].qmaps[iterations[i]].T₂     for i in 1:nrecons]
            means = [reconstructions[i].T₂.mean[:,iterations[i]]    for i in 1:nrecons]
            stds =  [reconstructions[i].T₂.std[:,iterations[i]]     for i in 1:nrecons]
        end

        for i in 1:nrecons

            # plot reconstructed map
            ax = subplot(2,ncols, i + (parcnt-1) * ncols);
                imshow(qmaps[i], cmap = cmaps[par], clim = clims[par]);
                # colorbar();
                title("$(descriptions[i]), Iteration $(iterations[i]-1) $par", color=clrtxt)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.spines["top"].set_visible(false)
                ax.spines["right"].set_visible(false)
                ax.spines["bottom"].set_visible(false)
                ax.spines["left"].set_visible(false)
                cb = colorbar()
                # cb = colorbar()
                cb.outline.set_edgecolor(clrtxt)
                cb.ax.yaxis.set_tick_params(color=clrtxt)
                plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color=clrtxt)

            ax = subplot(2,ncols, ncols + (parcnt-1)*ncols);


                # @show means[i]
                    bar( (1:ntubes) .+ barcenters[i], means[i], width = width, yerr = stds[i], color=barcolors[i], ecolor=clrerr, label=descriptions[i])
                # end

            #     # bar( (1:ntubes) .+ barcenters[end], groundtruth.T₂[294], width = width, color=barcolors[end], label="2D Mix")

            #     ax.xaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
            #     ax.yaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
            #     ax.tick_params(labelbottom=true, labeltop=false, labelleft=true, labelright=false, color="white")

            #     xlabel("# tube", color=clrtxt);
            #     ylabel("$(par) (s)", color=clrtxt);
            #     title("Mean $(par) values per tube", color=clrtxt);

            #     ylim(clims[par])

            #     ax.set_facecolor(clrbg)
            #     ax.spines["left"].set_color(clrtxt)
            #     ax.spines["bottom"].set_color(clrtxt)
            #     ax.tick_params(axis="x", colors=clrtxt)
            #     ax.tick_params(axis="y", colors=clrtxt)
            #     legend()
        end

        parcnt +=1
    end

    if has_reference

        ax = subplot(2,ncols, ncols);

        bar( (1:ntubes) .+ barcenters[end-1], vec(reference.T₁.mean), width = width, yerr = vec(reference.T₁.std), color=barcolors[end], ecolor=clrerr, label="Reference")

        ylim((0, maximum(reference.T₁.mean) * 1.1))

        legend(["EPG1", "EPG2", "Gold Standard"])

        ax = subplot(2,ncols, 2ncols);

        bar( (1:ntubes) .+ barcenters[end-1], vec(reference.T₂.mean), width = width, yerr = vec(reference.T₂.std), color=barcolors[end], ecolor=clrerr, label="Reference")

        ylim((0, maximum(reference.T₂.mean) * 1.1))

        legend(["EPG1", "EPG2", "Gold Standard"])

    end

    # ax = subplot(2,ncols, ncols);
    # bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth.T₁[temp][tube_idx], width = width, label="Ground-truth")
    # legend()
    # ax = subplot(2,ncols, 2*ncols);
    # bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth.T₂[temp][tube_idx], width = width, label="Ground-truth")
    # legend()

    fig.suptitle(figtitle, fontsize=16, color=clrtxt)
end

# Load goldstandard
function load_goldstandard_recons(path_T1, path_T2, new_x, new_y)

    # T₁
    T₁_data = MAT.matread(path_T1)["raw_data"] .|> ComplexF64;
    inversion_times = MAT.matread(path_T1)["TI"] |> vec;

    T₁_mask = abs.(T₁_data[:,:,1]) .> 0.1*maximum(abs.(T₁_data[:,:,1]));

    T₁, PD_T₁ = qMRI.goldstandard_T₁(T₁_data, T₁_mask, inversion_times);
    T₁[T₁ .> 7000] .= 0;
    T₁ ./= 1000;
    T₁ = reverse(T₁, dims=1);
    PD_T₁ = reverse(PD_T₁, dims=1);
    figure(); imshow(T₁, clim=(0,  5.0), cmap=:inferno)


    # T₂
    T₂_data = MAT.matread(path_T2)["raw_data"] .|> ComplexF64;
    echo_times = MAT.matread(path_T2)["TE"] |> vec;

    T₂_mask = abs.(T₂_data[:,:,1]) .> 0.1*maximum(abs.(T₂_data[:,:,1]));

    T₂, PD_T₂ = qMRI.goldstandard_T₂(T₂_data, T₂_mask, echo_times);
    T₂[T₂ .> 2000] .= 0;
    T₂ ./= 1000;
    T₂ = reverse(T₂, dims=1);
    PD_T₂ = reverse(PD_T₂, dims=1);
    figure(); imshow(T₂, clim=(0,  0.350), cmap=:parula)

    figure();
    subplot(121); imshow(PD_T₁)
    subplot(122); imshow(PD_T₂)

    PD = PD_T₂ # T2 is coil-combined, T1 is not

    # interpolate to higher resolution

    itp = interpolate(T₁, BSpline(Constant())) # "nearest"
    # interpolate slice to desired size
    sx,sy = size(T₁);
    x = LinRange(1,sx,new_x);
    y = LinRange(1,sy,new_y);
    T₁ = itp(x,y);

    # interpolate to higher resolution
    itp = interpolate(T₂, BSpline(Constant())) # "nearest"
    # interpolate slice to desired size
    sx,sy = size(T₂)
    x = LinRange(1,sx,new_x)
    y = LinRange(1,sy,new_y)
    T₂ = itp(x,y)

    # interpolate to higher resolution
    itp = interpolate(abs.(PD), BSpline(Constant())) # "nearest"
    # interpolate slice to desired size
    sx,sy = size(PD)
    x = LinRange(1,sx,new_x)
    y = LinRange(1,sy,new_y)
    PD = itp(x,y)
    PD = PD ./ maximum(PD);

    figure();
    subplot(131); imshow(T₁, clim=(0, 5.00), cmap=:inferno)
    subplot(132); imshow(T₂, clim=(0, 0.35), cmap=:parula)
    subplot(133); imshow(PD, clim=(0, 0.35), cmap=:gray)

    # remove some crap from PD
    PD[:,1:15] .= 0
    PD[:,end-9:end] .= 0

    mask = abs.(PD) .> 0.025 * maximum(abs.(PD));
    figure(); subplot(121); imshow(mask)

    PD[.!mask] .= 0

    # Generate mask per tube
    gs_tubes = mask_per_tube(mask, 60);

    # Check mask per tube
    figure()
    tmp = zero(abs.(T₂))
    # pd_max = maximum(pd);
    for i = eachindex(gs_tubes)
        tmp[gs_tubes[i]] .= i;
    end

    subplot(122); imshow(tmp);

    gs = (qmaps = [StructArray(T₁=T₁, T₂=T₂, ρ=complex.(PD))],);
    gs = add_mean_and_std_per_tube(gs);

end

function bias_precision(reconstructions::Vector{<:NamedTuple}, goldstandard, descriptions)

    GS = goldstandard

    bias = map(reconstructions) do r
        (   T₁ = 100 * abs.(r.T₁.mean .- GS.T₁.mean) ./ GS.T₁.mean,
            T₂ = 100 * abs.(r.T₂.mean .- GS.T₂.mean) ./ GS.T₂.mean  )
    end

    standard_deviation = map(reconstructions) do r
        (   T₁ = 100 * r.T₁.std ./ r.T₁.mean,
            T₂ = 100 * r.T₂.std ./ r.T₂.mean    )
    end

    mean_bias = [vec(mean(vcat(bias[i].T₁,bias[i].T₂),dims=1)) for i in 1:length(reconstructions)];
    mean_std = [vec(mean(vcat(standard_deviation[i].T₁,standard_deviation[i].T₂),dims=1)) for i in 1:length(reconstructions)];

    _,cartesian_index = findmin(map((x,y)->norm((x,y)), mean_bias[1], mean_std[1]))
    _,radial_index = findmin(map((x,y)->norm((x,y)), mean_bias[2], mean_std[2]))

    figure()

        for (i,r) in enumerate(reconstructions)
            plot(mean_std[i], mean_bias[i], marker="o")
        end
        xlabel("Mean standard deviation (%)")
        ylabel("Mean bias (%)")
        title("Mean bias & standard deviation")
        xlim([0,25])
        ylim([0,50])

        # lgnd = mapreduce(d->["T₁ $d", "T₂ $d"], vcat, descriptions)
        lgnd = ["Cartesian", "Radial"]
        legend(lgnd)
        grid(true)


    figure()

        for j = 1:length(GS.T₁.mean)

            subplot(3,4,j)

            for (i,r) in enumerate(reconstructions)
                plot(vec(standard_deviation[i].T₁[j,:]), vec(bias[i].T₁[j,:]), marker="o")
                plot(vec(standard_deviation[i].T₂[j,:]), vec(bias[i].T₂[j,:]), marker="o")
            end

            xlabel("Standard Deviation (%)")
            ylabel("Bias (%)")
            # title("Mean bias & standard deviation for tube $j")
            # ylim([0,25])

            if j == 1
                lgnd = mapreduce(d->["T₁ $d", "T₂ $d"], vcat, descriptions)
                legend(lgnd)
            end
        end

        return cartesian_index, radial_index
end


# function barplot_tube_recons(reconstructions, iterations, descriptions; figtitle="MR-STAT", black_theme=false)

#     # select iteration for each of the different reconstructions
#     # recons = map((recon,iter)->recon[iter].qmaps, reconstructions, iterations)

#     recons = map((recon,iter)->recon[iter], reconstructions, iterations)

#     # for barplots
#     ntubes = length(recons[1].tubes.T₁.mean)
#     width = 0.2
#     nrecons = length(iterations)
#     totalwidth = (nrecons+1) * width
#     barcenters = LinRange(-totalwidth/2 + width/2, totalwidth/2 - width/2, nrecons + 1)
#     # barcolors = ["#D36060" ,"#669BBC", "#264653", "#FFD275", "#2A9D8F", "#E85D75", "#E1CE7A"]
#     # barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E1CE7A"]

#     temp = 294
#     tube_idx = [8,3,14,15,6,10,7,12,13,11,16,17]
#     # tube_idx = [11,7,17,8,15,6,13,10,12,16,14,3]

#     # choose clims and cmaps for each reconstruction parameter
#     clims = Dict(:T₁ => (0.0, 2.0), :T₂ => (0.0, 0.4), :ρ => (0.0, 2.0))
#     cmaps = Dict(:T₁ => "inferno", :T₂ => "parula",   :ρ => "gray")

#     if black_theme
#         clrbg  = "black"
#         clrtxt = "white"
#         clrerr = "lightgray"
#         barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E1CE7A"]
#     else
#         clrbg = "white"
#         clrtxt = "black"
#         clrerr = "black"
#         barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E1CE7A"]
#     end

#     # barcolors = ["#D36060" ,"#669BBC", "#264653", "#FFD275", "#2A9D8F", "#E85D75", "#E1CE7A"]

#     fig = figure(facecolor=clrbg)

#     parcnt = 1
#     nrows = 1
#     ncols = 2

#     # make plots, one figure for each reconstruction parameter
#     for par in [:T₁, :T₂]

#         stats = [@eval $recon.tubes.$par for recon ∈ recons]

#         for i in 1:length(recons)

#             ax = subplot(nrows,ncols,parcnt);

#                 bar( (1:ntubes) .+ barcenters[i], stats[i].mean, width = width, yerr = stats[i].std, color=barcolors[i], ecolor=clrerr, label=descriptions[i])

#                 ax.xaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
#                 ax.yaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
#                 ax.tick_params(labelbottom=true, labeltop=false, labelleft=true, labelright=false, color="white")

#                 xlabel("# tube", color=clrtxt);
#                 ylabel("$(par) (s)", color=clrtxt);
#                 title("Mean $(par) values per tube", color=clrtxt);

#                 ylim(clims[par])

#                 ax.set_facecolor(clrbg)
#                 ax.spines["left"].set_color(clrtxt)
#                 ax.spines["bottom"].set_color(clrtxt)
#                 ax.tick_params(axis="x", colors=clrtxt)
#                 ax.tick_params(axis="y", colors=clrtxt)
#                 legend()
#         end

#         parcnt +=1
#     end

#     # add ground truth values
#     ax = subplot(nrows,ncols, 1);
#     bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth_T₁[temp][tube_idx], width = width, label="Ground-truth")
#     legend()
#     ax = subplot(nrows,ncols, 2);
#     bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth_T₂[temp][tube_idx], width = width, label="Ground-truth")
#     legend()

#     fig.suptitle(figtitle, fontsize=16, color=clrtxt)
# end


# ##

# function barplot_compare_iterations(recons, descriptions, black_theme=false)

#     # select iteration for each of the different reconstructions
#     maxits = maximum(length(r) for r in recons)

#     # for barplots
#     ntubes = length(recons[1][1].tubes.T₁.mean)
#     width = 0.2
#     nrecons = length(recons)
#     totalwidth = (nrecons+1) * width
#     barcenters = LinRange(-totalwidth/2 + width/2, totalwidth/2 - width/2, nrecons + 1)

#     temp = 294
#     # tube_idx = [8,3,14,15,6,10,7,12,13,11,16,17]
#     tube_idx = [11,13,7,10,17,12,8,16,15,14,6,3]


#     # choose clims and cmaps for each reconstruction parameter
#     clims = Dict(:T₁ => (0.0, 2.0), :T₂ => (0.0, 0.4), :ρ => (0.0, 2.0))
#     cmaps = Dict(:T₁ => "inferno", :T₂ => "parula",   :ρ => "gray")

#     if black_theme
#         clrbg  = "black"
#         clrtxt = "white"
#         clrerr = "lightgray"
#         barcolors = ["#E2798E" ,"#861D32", "#669BBC", "#264653", "#2A9D8F", "#E85D75", "#E08E45"]
#     else
#         clrbg = "white"
#         clrtxt = "black"
#         clrerr = "black"
#         barcolors = ["#E2798E" ,"#669BBC", "#861D32", "#264653", "#2A9D8F", "#E85D75", "#E08E45"]
#     end

#     for par in [:T₁, :T₂]
#         fig = figure(facecolor=clrbg)

#         for tube in 1:ntubes

#         ax = subplot(6,2,tube)
#         parcnt = 1
#         nrows = 1
#         ncols = 2


#         # make plots, one figure for each reconstruction parameter


#             for i in 1:length(recons)

#                 recon = recons[i]

#                 @eval means = round.(Int, 1000 * [$recon[it].tubes.$(par).mean[$(tube)] for it in 1:$(maxits)])
#                 @eval stds  = round.(Int, 1000 * [$recon[it].tubes.$(par).std[$(tube)] for it in 1:$(maxits)])

#                 # ax = subplot(nrows,ncols,parcnt);

#                 bar( (2:maxits) .+ barcenters[i], means[2:end], width = width, yerr = stds[2:end], color=barcolors[i], ecolor=clrerr, label=descriptions[i])
#                 legend()

#                 plot(0.5:1:maxits+0.5, fill(means[end], maxits+1), linestyle="dotted", color=barcolors[i])

#                 ax.xaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
#                 ax.yaxis.set_tick_params(which="both", bottom=true, top=false, color="white")
#                 ax.tick_params(labelbottom=true, labeltop=false, labelleft=true, labelright=false, color="white")

#                 xlabel("Iteration", color=clrtxt);
#                 ylabel("Tube $tube", color=clrtxt);
#                 # title("Mean $(par) values per tube", color=clrtxt);

#                 # ylim(clims[par])

#                 # ax.set_facecolor(clrbg)
#                 # ax.spines["left"].set_color(clrtxt)
#                 # ax.spines["bottom"].set_color(clrtxt)
#                 # ax.tick_params(axis="x", colors=clrtxt)
#                 # ax.tick_params(axis="y", colors=clrtxt)
#                 # legend()
#             end

#             @eval gt = groundtruth.$(par)[$temp]

#             plot(0.5:1:maxits+0.5, fill(gt[tube_idx[tube]], maxits+1), linestyle="dashed", color=barcolors[end], label="Tabulated")
#             # bar( (1:maxits) .+ barcenters[end],gt[tube_idx[tube]] * ones(maxits), width = width, label="Ground-truth", color=barcolors[end-3])
#             legend()
#             suptitle("$par: mean and std per iteration per tube")
#             # parcnt +=1
#         end
#     end

#     # add ground truth values
#     ax = subplot(nrows,ncols, 1);
#     bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth.T₁[temp][tube_idx], width = width, label="Ground-truth")
#     legend()
#     ax = subplot(nrows,ncols, 2);
#     bar( (1:ntubes) .+ barcenters[end], 10^-3 * groundtruth.T₂[temp][tube_idx], width = width, label="Ground-truth")
#     legend()

#     # fig.suptitle(figtitle, fontsize=16, color=clrtxt)
# end

# ## Compute mean and std values per tube for each iteration of each reconstruction

# function tube_statistics(qmaps::StructArray, tubes)

#     #compute mean T₁ and T₂ values and standard deviations
#     recon = (
#         qmaps = qmaps,
#         tubes = (
#             T₁ = (
#                 mean = [mean(qmaps.T₁[tube]) for tube in tubes],
#                 std = [std(qmaps.T₁[tube]) for tube in tubes]
#             ),
#             T₂ = (
#                 mean = [mean(qmaps.T₂[tube]) for tube in tubes],
#                 std = [std(qmaps.T₂[tube]) for tube in tubes]
#             ),
#             ρ = (
#                 mean = [mean(abs.(qmaps.ρ[tube])) for tube in tubes],
#                 std = [std(abs.(qmaps.ρ[tube])) for tube in tubes]
#             )
#         )
#     );

#     return recon
# end