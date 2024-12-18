include("../startup.jl")
include("ScanAnalysis.jl")

sx = 224; sy = 224; nSweeps = 6
nZoom = 1

function figuresReadData!(figurePars)
    folderDetails = figurePars["folderDetails"]      # vector of strings referring to the sequence type    
    fn_common_set = figurePars["fn_common_set"]      # vector of (typically 2 or 3) folder names (refering to recon types); length of 1 is accepted 
    seqDescription= figurePars["seqDescription"]     # vector of description of sequences
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    ntypes        = figurePars["ntypes"]             # vector of numbers 3 (set is not supposed to contain co-reconed B1) or 4 (does contain B1)
    nRecons       = figurePars["nRecons"]            # How many reconstructions are minimally available 
    sliceRange    = figurePars["sliceRange"]         # The range of slices that are available 
                                                        # Note: overloaded with recons dimension, cannot be both >1
    maxIt         = figurePars["maxIt"]               # How many iterations are available in the set containing least of them 
    gels          = figurePars["gels"]               # The set of vials to be used
    comparisonSet = figurePars["comparisonSet"]      # vector of two numbers, referring to indices of folder names; typically [1,2]
    
    # Read data into 6-d array
    mapSet = []
    dreamSet= []
    outputs = []
    dimRange = (length(sliceRange) > 1) ? sliceRange : (1:nRecons)

    testIt = zeros(Int64,length(fn_common_set),length(folderDetails))
    RFdeg  = zeros(length(folderDetails), sy*nSweeps)
    RFphaseDeg = similar(RFdeg)
    for (recType, fn_common) in enumerate(fn_common_set)
        dream= []
        maps = []
        maps = zeros(length(folderDetails), dimRange,4,maxIt,sx,sy);
        maps[:,:,4,:,:,:] .= 1.0  # initialize B1 to 1.0 instead of zero
        if (length(sliceRange) > 1)
            dream= zeros(length(folderDetails),sx,sy, length(sliceRange))
            # maps = zeros(length(folderDetails), dimRange,4,maxIt,sx,sy, length(sliceRange));
            # maps[:,:,4,:,:,:,:] .= 1.0  # initialize B1 to 1.0 instead of zero
        else
            dream= zeros(length(folderDetails),sx,sy)
        end

        for fff in eachindex(folderDetails)
            fn = fn_common*folderDetails[fff]
            folder = fn
            if (length(sliceRange) > 1)
                # Multiple slice, so the 'rec' dimension of the array is actually used for slices;
                # slice information resides one folder deeper in the multi-slice case 
                # If several recons are available, only the last is taken into account
                folder = [sub for sub in filter(x -> isdir(joinpath(fn, x)), readdir(fn, join=true))][end]
            end

            outputs = get_old_recons(folder)
            for rec in dimRange
                for it in eachindex(outputs[rec].qmaps); if it<=maxIt
                    for type = 1:ntypes[recType] # 1:ntypes[fff]
                        @show recType, fff, rec, type, it
                        maps[fff,rec,type,it,:,:] = ExtractImage(outputs, type, rec, it)
                    end
                end; end
            end
            testIt[recType,fff] = min(length(outputs[1].qmaps),maxIt)

            # Subsequent four lines worked with UMCU-export but not if I have reconframe2mrstat format
            # U = matread(fn*"/U.mat")
            # dream[fff,:,:] = abs.(U["B1"])
            # RFdeg[fff,:]      = U["Sequence"]["Flip_angle"][1,:]
            # RFphaseDeg[fff,:] = U["Sequence"]["RF_phase"][1,:]

            if (length(sliceRange) > 1)
                for s in 1:length(sliceRange)
                    dream[fff,:,:,s] = abs.(outputs[s].B₁map)
                end
            else
                dream[fff,:,:] = abs.(outputs[1].B₁map)
            end
            RFdeg[fff,:]      =   abs.(outputs[1].sequence.RF_train)
            RFphaseDeg[fff,:] = angle.(outputs[1].sequence.RF_train).*180/π
            
        end
        # maps = fftshift(maps,6) # unfortunately, depends very much on recon version
        reverse!(maps,dims=5)
        #reverse!(dream,dims=1)  # seen on 2024-06-27 as pretty dubious and inactive in most cases; is dims=2 intended??
        reverse!(dream,dims=2)
        push!(mapSet,maps)
        push!(dreamSet,dream)
    end

    # means and standard deviations
    dynDim = (length(sliceRange) > 1) ? 99 : 2
    meanoverdyn= [mean(mapSet[i],dims=dynDim) for i in eachindex(mapSet)];
    devoverdyn = [std(mapSet[i], dims=dynDim) for i in eachindex(mapSet)]; # case [dummy] type iter x y
    # devoverdyn = meanoverdyn
    # Measured by Oscar using MIR-MSE measurements
    goldstandard_3T = (
            #         1    2    3    4    5    6    7    8    9    10   11   12    13    14    15    16    17    18
        T1 = Float64[216, 323, 316, 500, 488, 480, 646, 629, 814, 801, 965, 1645, 1047, 1120, 1266, 1431, 1394, 1555],
        T2 = Float64[ 44,  61,  96,  45,  78, 128,  78, 109,  96, 127, 108,  303,  179,  132,  167,  161,  145,  140]
    )

    # manual tube placement
    # symbolicPos = [(0,1), (0,2), (1,0), (1,1), (1,2), (1,3), (2,0), (2,1), (2,2), (2,3), (3,1), (3,2), (2.6,2.6), (1.5,1.5)];
    # tubeLabel   = [   1,     2,     3,     4,     5,     6,     7,    8,    16,     18,     14,    15,          0,        0]
    # origin = (35,71)
    # distance = 44.1
    # angulation = -0.23
    # radius = 8 
    # refSlice = 10
    symbolicPos = [(0,1), (0,2), (1,0), (1,1), (1,2), (1,3), (2,0), (2,1), (2,2), (2,3), (3,1), (3,2), (2.6,2.6), (1.5,1.5)];
    tubeLabel   = [  14,     9,    11,     1,     2,     4,    18,   10,    13,     16,      5,     7,          0,        0]
    origin = (67,33.5)
    distance = 43.8
    angulation = 0.24
    radius = 7 
    refSlice = 1 # 10
    usableSliceRange = 8:14 # 1:1 # 10:10 # 8:14
    centerPos = [origin.+(distance*x*cos(angulation)-distance*y*sin(angulation), distance*x*sin(angulation)+distance*y*cos(angulation)) for (x,y) in symbolicPos]
    tube_centers = [CartesianIndex(round(Int64,x),round(Int64,y)) for (x,y) in centerPos]

    refRec = (length(sliceRange) > 1) ? refSlice : 1
    ρ = meanoverdyn[1][1,refRec,3,testIt[1,1],:,:]

    d(x::CartesianIndex, y::CartesianIndex) = √( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    allidx = CartesianIndices(ρ)

    tubes = [allidx[ d.((center,), allidx) .< radius] for center ∈ tube_centers];
    ShowTubes(ρ, tubes)

    # Take mean values over ROIs
    indexRange= (length(sliceRange) > 1) ? sliceRange : (1:1)
    meanMeans = zeros(length(fn_common_set),length(folderDetails),2,length(tubes))
    meanDevs  = zeros(length(fn_common_set),length(folderDetails),2,length(tubes))
    devDevs   = zeros(length(fn_common_set),length(folderDetails),2,length(tubes)) # disused
    for recType in eachindex(fn_common_set)
        for case in 1:length(folderDetails)
            for m in 1:2
                # case zoom [dummy] type iter x y
                meanIm = meanoverdyn[recType][case,indexRange,m,testIt[recType,case],:,:] 
                devIm  =  devoverdyn[recType][case,indexRange,m,testIt[recType,case],:,:] 
                for (i,tube) in enumerate(tubes)
                    @show recType, case, m, i
                    meanMeans[recType,case,m,i] = mean(meanIm[usableSliceRange,tube])
                    meanDevs[recType,case,m,i] = mean(devIm[usableSliceRange,tube])
                end        
            end
        end
    end

    # out: meanoverdyn, meanMeans, meanDevs, goldstandard_3T, tubes, dream, RFdeg, RFphaseDeg
    figurePars["meanoverdyn"]     = meanoverdyn
    figurePars["meanMeans"]       = meanMeans
    figurePars["meanDevs"]        = meanDevs
    figurePars["goldstandard_3T"] = goldstandard_3T
    figurePars["tubes"]           = tubes
    figurePars["dreamSet"]        = dreamSet
    figurePars["RFdeg"]           = RFdeg
    figurePars["RFphaseDeg"]      = RFphaseDeg
    figurePars["testIt"]          = testIt
    figurePars["tubeLabel"]       = tubeLabel
    figurePars["usableSliceRange"]= usableSliceRange
end

function figuresPhantomFig1(figurePars)
    folderDetails = figurePars["folderDetails"]      # vector of strings referring to the sequence type    
    meanoverdyn   = figurePars["meanoverdyn"]     
    dreamSet      = figurePars["dreamSet"]          
    RFdeg         = figurePars["RFdeg"]           
    RFphaseDeg    = figurePars["RFphaseDeg"] 
    testIt        = figurePars["testIt"]
    dispMax       = figurePars["dispMax"]
    sliceRange    = figurePars["sliceRange"]
    usableSliceRange = figurePars["usableSliceRange"] 
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = length(sliceRange) > 1 ? midSlice : 1

    recType = 1

    lines_color_cycle = [p["color"] for p in plt.rcParams["axes.prop_cycle"]]

    wRat = ones(length(folderDetails))
    wRat[end] = 1.25
    ddd=Dict("width_ratios" => wRat)
    fig_seqs,ax_seqs = subplots(1,length(folderDetails),figsize=(8,3))
    fig_seqs.suptitle("RF sequences")
    fig_seqs.subplots_adjust(wspace=0.25,hspace=0.0,bottom=0.2)

    fig_maps,ax_maps = subplots(2,length(folderDetails),figsize=(10,8),gridspec_kw=ddd)
    fig_maps.subplots_adjust(wspace=-0.1,hspace=0.0)
    fig_maps.suptitle("T1 and T2 maps of "*figurePars["recDescription"][recType])

    fig_B1s, ax_B1s  = subplots(1,3,figsize=(8,3),gridspec_kw=Dict("width_ratios" => [1, 1, 1.25]))
    fig_B1s.subplots_adjust(wspace=-0.1,hspace=0.0)
    fig_B1s.suptitle("B1 maps: unmodified, modified and no correction")

    for case in eachindex(folderDetails)
        # plot RF sequence
        ax_seqs[case].set_ylim(-10.0,90.0)
        ax_seqs[case].plot(RFdeg[case,:],label="amplitude")
        anglesdd = zeros(size(RFdeg)[2])
        for i in 1:size(RFdeg)[2]-2
            anglesdd[i] = ((-RFphaseDeg[case,i]+2*RFphaseDeg[case,i+1]-RFphaseDeg[case,i+2]+720.0+270.0)) % 180.0 - 90.0
        end
        ax_seqs[case].plot(anglesdd,label="ϕ''")  
        ax_seqs[case].set_xlabel("pulse number")
        ax_seqs[1].set_ylabel("amplitude [deg] or ϕ'' [deg/TR²]")
        
        # plot T1 and T2 maps reconstructed via first recon type
        for m in 1:2
            meanIm = meanoverdyn[recType][case,recPick,m,testIt[recType,case],:,:] 

            imClip, rgb_vec = relaxationColorMap("T$m", meanIm, 0.0, dispMax[m])  # call to resource, generating a colormap 
            cmap = PyPlot.ColorMap("relaxationColor", rgb_vec, length(rgb_vec), 1.0) 
            pcm = ax_maps[m,case].imshow(imClip, vmax=dispMax[m], cmap=cmap);
            if case==length(folderDetails)
                clb = colorbar(pcm,ax=ax_maps[m,case], shrink=0.8, ticks = [0, dispMax[m]]); 
                clb.ax.tick_params(labelsize=20) 

                dispMaxVal = 1000*dispMax[m]
                clb.ax.set_yticklabels(["0ms", "$dispMaxVal ms"], fontsize=14)
            end
            ax_maps[m,case].set_xticks([]); 
            ax_maps[m,case].set_yticks([]);
            ax_maps[m,case].set_yticklabels([]);  
        end
    end
    ax_seqs[2].legend()

    # plot Dream maps 
    B1map = dreamSet[1][1,:,:,midSlice]
    maskOnes = ones(size(B1map)) .* (B1map .> 0.01)

    vminB1 = 0.75; vmaxB1=1.25
    pcm = ax_B1s[1].imshow(B1map,  vmin=vminB1,vmax=vmaxB1)
    ax_B1s[2].imshow(1.05 .* B1map,vmin=vminB1,vmax=vmaxB1)
    ax_B1s[3].imshow(maskOnes,     vmin=vminB1,vmax=vmaxB1)
    clb = colorbar(pcm,ax=ax_B1s[3], shrink=0.8); 
    clb.ax.tick_params(labelsize=14) 
    for i in 1:3
        ax_B1s[i].set_xticks([]); 
        ax_B1s[i].set_yticks([]);
        ax_B1s[i].set_yticklabels([]);  
    end

    # txt_seqs = ["(a)", "(b)", "(c)"]
    # txt_maps = ["(d)" "(e)" "(f)"; "(g)" "(h)" "(i)"; "(j)" "(k)" "(l)"]
    # for i in 1:3
    #     ax_seqs[i].text(0,0.9,txt_seqs[i], fontsize=18, transform=ax_seqs[i].transAxes)
    #     for j in 1:3
    #         ax_maps[i,j].text(0.05,0.85,txt_maps[i,j], fontsize=18, transform=ax_maps[i,j].transAxes, color="white")
    #     end
    # end
end

function figuresPhantomFig2plus(figurePars, diffusion_corrected=false)
    meanoverdyn = figurePars["meanoverdyn"]  
    meanMeans   = figurePars["meanMeans"]       
    meanDevs    = figurePars["meanDevs"]
    goldstandard_3T=figurePars["goldstandard_3T"]
    tubes       = figurePars["tubes"]
    tubeLabel   = figurePars["tubeLabel"] 
    folderDetails = figurePars["folderDetails"]      # vector of strings referring to the sequence type    
    gels          = figurePars["gels"]               # The set of vials to be used
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    comparisonSet = figurePars["comparisonSet"]      # vector of two numbers, referring to indices of folder names; typically [1,2]
    diffuCorrFile = figurePars["diffuCorrFile"]      # location/name of file containing diffusion correction data
    dctxt = ""
    if diffusion_corrected
        @load diffuCorrFile mapCollection
        dctxt = " (diffusion-corrected)"
    end

    lines_color_cycle = [p["color"] for p in plt.rcParams["axes.prop_cycle"]]

    # Prepare Bias list and relative bias 
    type = 2 # T2
    gtT2 = goldstandard_3T[type]
    xxx = ([gtT2[tubeLabel[gel]] for gel in gels])
    biasT2    = zeros(length(recDescription),length(folderDetails),length(gels))
    biasT2rel = zeros(length(recDescription),length(folderDetails),length(gels))
    for recType in eachindex(recDescription)
        for case in eachindex(folderDetails)
            biasT2[recType,case,:] = 1000.0 .* meanMeans[recType,case,type,gels].-xxx
            if diffusion_corrected
                diffCorr = [mapCollection[case][type,tubeLabel[gel]] for gel in gels]
                biasT2[recType,case,:] = biasT2[recType,case,:] .+ 1000.0 .* diffCorr
            end 
            biasT2rel[recType,case,:] = biasT2[recType,case,:]./xxx
        end
    end

    fig,ax = subplots(1,3,figsize=(17,6))
    markers = ["o","x","+"]
    reconNames = recDescription

    m = 2
    for case in eachindex(folderDetails)
        label = seqDescription[case] # chop(folderDetails[case],head=1,tail=0)
        color = lines_color_cycle[case]
        yyy = 1000.0 .* (meanMeans[comparisonSet[2],case,m,gels].-meanMeans[comparisonSet[1],case,m,gels])
        ax[1].scatter(xxx,yyy,label=label, color=color)
        slope = mean(yyy) / mean(xxx)
        @show label, slope 
        xl = minimum(xxx); xh = maximum(xxx)
        ax[1].plot([xl,xh],slope.*[xl,xh], color=color)        
    end
    ax[1].set_xlabel("Gold standard T2 [ms]")
    ytext = "T2 estimation difference [ms] due to 5% offset in B1"*dctxt
    ax[1].set_ylabel(ytext)
    ax[1].legend()

    for recType in eachindex(recDescription)
        for case in eachindex(folderDetails)
            label = ""
            label = recType==1 ? folderDetails[case] : (case==1 ? recDescription[recType] : "") 
            color = lines_color_cycle[case]
            marker = markers[recType]
            #yyy = 1000.0 .* meanMeans[recType,case,m,gels].-xxx
            yyy = biasT2[recType,case,:]
            ax[2].scatter(xxx,yyy,label=label, color=color, marker=marker)
        end
    end
    ax[2].set_xlabel("Gold standard T2 [ms]")
    ax[2].set_ylabel("T2 estimation error [ms]"*dctxt)
    ax[2].legend()

    # Alessandro's suggested deviation-bar graphs
    spacer = 3 + 2
    for recType in eachindex(recDescription)
        for case in eachindex(folderDetails)
            # pos = 3*(case-1) + recType-1
            pos = spacer*(recType-1) + (case-1)
            ave = 100.0 * mean(biasT2rel[recType,case,:])
            dev = 100.0 * std(biasT2rel[recType,case,:])
            color = lines_color_cycle[case]
            ax[3].bar(pos,ave,color=color)
            ax[3].errorbar(pos,ave,dev,linewidth=4.0,capsize=8.0,color="black")
        end
    end
    ax[3].set_xticks([]); 
    ax[3].plot([-0.8,(3-1)*spacer+4.0],[0.0,0.0],color="black")
    ax[3].yaxis.set_label_position("right")
    ax[3].yaxis.tick_right()
    ax[3].set_ylabel("T2 bias in %", fontsize=18)
    
    for case in eachindex(folderDetails)
        color = lines_color_cycle[case]
        ax[3].text(0.0,-13.0-2*case, seqDescription[case], color=color, fontsize = 14)
    end 
    
    yPositions = [3,-6,1]
    for recType in eachindex(recDescription)
        ax[3].text(spacer*(recType-1)-1.2,yPositions[recType],recDescription[recType], fontsize=12)
    end

    # Prepare statistics for plotting of noise bars 
    goldValues = copy(meanMeans)
    for case in eachindex(folderDetails)
        for r in eachindex(recDescription)
            for m in 1:2
                goldValues[r,case,m,gels] = (0.001 .* [goldstandard_3T[m][tubeLabel[gel]] for gel in gels])
            end
        end
    end
    goldSpread = std(goldValues[:,:,:,gels], dims=4)
    reconSpread = std(meanMeans[:,:,:,gels],dims=4)

    qqq = meanDevs./goldValues
    # Correct for bias in results 
    for case in eachindex(folderDetails)
        for r in eachindex(fn_common_set)
            for m in 1:2
                thisSpread = reconSpread[r,case,m,1]
                qqq[r,case,m,:] = qqq[r,case,m,:] .* goldSpread[r,case,m,1] ./ thisSpread 
            end
        end
    end    

    if ~any(isnan,qqq)
        pos = 1:length(folderDetails)
        # Noise bars 
        figure()
        for m in 1:2
            qqqMean = mean(qqq[:,:,:,gels],dims=4)
            qqqMeanAct = qqqMean[1,:,m,1]
            bar(pos.+0.3*(m-1),qqqMeanAct)
            for case in eachindex(folderDetails)
                @printf("Relative noise level of of T%d of sequence %s is %.3f\n", m, folderDetails[case], qqqMeanAct[case])
            end
        end
        xticks(pos,seqDescription)
    end
end

function figureAllIm(figurePars)
    folderDetails = figurePars["folderDetails"]      # vector of strings referring to the sequence type    
    fn_common_set = figurePars["fn_common_set"]
    meanoverdyn   = figurePars["meanoverdyn"]     
    testIt        = figurePars["testIt"]
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    dispMax       = figurePars["dispMax"]
    sliceRange    = figurePars["sliceRange"]
    usableSliceRange = figurePars["usableSliceRange"] 
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = length(sliceRange) > 1 ? midSlice : 1

    # (all of T1 maps and T2 maps for scan*recon combinations)

    wRat = ones(length(folderDetails))
    wRat[end] = 1.25
    ddd=Dict("width_ratios" => wRat)

    for m in 1:2
        fig_maps,ax_maps = subplots(length(fn_common_set),length(folderDetails),figsize=(7,6),gridspec_kw=ddd)
        fig_maps.subplots_adjust(wspace=-0.1,hspace=0.0,right=0.8, top=1.0, bottom=0.0, left=0.0)
        for case in eachindex(folderDetails)
            for r in eachindex(fn_common_set)
                # row = 3*(m-1)+r
                row = r
                # case zoom [dummy] type iter x y
                meanIm = meanoverdyn[r][case,recPick,m,testIt[r,case],:,:] 

                imClip, rgb_vec = relaxationColorMap("T$m", meanIm, 0.0, dispMax[m])  # call to resource, generating a colormap 
                cmap = PyPlot.ColorMap("relaxationColor", rgb_vec, length(rgb_vec), 1.0) 
                pcm = ax_maps[row,case].imshow(imClip, vmax=dispMax[m], cmap=cmap);
                if case==length(folderDetails)
                    clb = colorbar(pcm,ax=ax_maps[row,case], shrink=0.5, ticks = [0, dispMax[m]]); 
                    clb.ax.tick_params(labelsize=20) 
                    #      clb = colorbar(pcm, ax=thisax, shrink=0.6, ticks=[loLev, (upLev+loLev)*0.5, upLev])

                    dispMaxVal = 1000*dispMax[m]
                    clb.ax.set_yticklabels(["0ms", "$dispMaxVal ms"], fontsize=14)
                end
                ax_maps[row,case].set_xticks([]); 
                ax_maps[row,case].set_yticks([]);
                ax_maps[row,case].set_yticklabels([]);  
            end
            ax_maps[1,case].set_title(seqDescription[case], fontsize=16)
        end
        ax_maps[2,2].text(0.0,0.95,"T$m", fontsize=24, transform=ax_maps[2,2].transAxes, color="white")
    end

    # casenames = seqDescription
    # for case in 1:3
    #     ax_maps[1,case].text(0.12,0.85,casenames[case], fontsize=14, transform=ax_maps[1,case].transAxes, color="white")
    # end

    # reconnames = recDescription
    # reconTextYstart=[0.05, 0.05, 0.05]
    # for i in 1:6
    #     r = (i-1)%3+1
    #     ax_maps[i,1].text(0.01,reconTextYstart[r],reconnames[r], fontsize=14, transform=ax_maps[i,1].transAxes, rotation="vertical",color="white")
    # end
end   


function figureDifferences(figurePars)
    folderDetails = figurePars["folderDetails"]      # vector of strings referring to the sequence type    
    meanoverdyn   = figurePars["meanoverdyn"]     
    testIt        = figurePars["testIt"]
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    dispMax       = figurePars["dispMax"]
    sliceRange    = figurePars["sliceRange"]
    usableSliceRange = figurePars["usableSliceRange"] 
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = length(sliceRange) > 1 ? midSlice : 1

    # For the T2 maps, differences between the first reconstruction mode and the rest  
    if length(recDescription) > 1   # only meaningful if there is a 'rest'
        wRat = ones(length(folderDetails))
        wRat[end] = 1.25
        ddd=Dict("width_ratios" => wRat)

        fig_rows = min(2,length(recDescription)-1)
        fig_maps,ax_maps = subplots(2,length(folderDetails),figsize=(7,8),gridspec_kw=ddd)  
        fig_maps.subplots_adjust(wspace=-0.1,hspace=0.0,right=0.9, top=1.0, bottom=0.0, left=0.0)

        for case in eachindex(folderDetails)
            m = 2
            meanImRef = meanoverdyn[1][case,recPick,m,testIt[1,case],:,:] 
            for r in 2:length(recDescription)
                row = r-1
                # case zoom [dummy] type iter x y
                meanIm = meanoverdyn[r][case,recPick,m,testIt[r,case],:,:] 

                #vmax = 0.05*dispMax[m]; vmin = -vmax;
                vmax = 0.012; vmin = -vmax;

                pcm = ax_maps[row,case].imshow(meanIm .- meanImRef, vmin=vmin, vmax=vmax, cmap="RdBu");
                if case==length(folderDetails)
                    clb = colorbar(pcm,ax=ax_maps[row,case], shrink=0.4, ticks = [vmin, 0, vmax]); 
                    clb.ax.tick_params(labelsize=20) 
                    #      clb = colorbar(pcm, ax=thisax, shrink=0.6, ticks=[loLev, (upLev+loLev)*0.5, upLev])

                    clb.ax.set_yticklabels(["$(1000*vmin) ms", "0 ms", "$(1000*vmax) ms"], fontsize=14)
                end
                ax_maps[row,case].set_xticks([]); 
                ax_maps[row,case].set_yticks([]);
                ax_maps[row,case].set_yticklabels([]);  
            end
            ax_maps[1,case].set_title(seqDescription[case], fontsize=16)    
        end
    end

end
