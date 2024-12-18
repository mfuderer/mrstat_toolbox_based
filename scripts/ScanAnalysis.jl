using Statistics
include("../plot/tubes.jl")

function GetMultipleReconsets(recon_folder_list)
    dirList = [filter(x -> isdir(joinpath(recfold, x)), readdir(recfold, join=true)) for recfold in recon_folder_list]
    return [[load(joinpath(recon, "mrstat_output.jld2"))["mrstat_output"] for recon in thisdir] for thisdir in dirList] 
end

function get_old_recons(recfold)
    #return [load(joinpath(recon, "mrstat_output.jld2"))["mrstat_output"] for recon in readdir(recfold, join=true)]
    return [load(joinpath(recon, "mrstat_output.jld2"))["mrstat_output"] 
                for recon in filter(x -> isdir(joinpath(recfold, x)), readdir(recfold, join=true))]
end

function GetRequiredRecons()
    outputs = get_old_recons(saParams["tag"]);      
    return outputs;  
end

function ShowSpectrum(outputs)
    nCases = saParams["nCases"]
    nRecons = saParams["nRecons"]
    latestIts = [length(outputs[case][i].qmaps) for i in 1:nRecons for case in 1:nCases]
    maxIt = minimum(latestIts) 
    largestIt = maximum(latestIts)
    casenames = saParams["casenames"]
    colortag = saParams["colortag"]
    mapnames = ["T1", "T2", "rho"]
    cset     = saParams["displSet"] 


    # Reshape data into 6-d array
    (sx,sy) = size(outputs[1][1].qmaps[1].T₁)
    maps = zeros(nCases,nRecons,3,largestIt,sx,sy);

    for case in 1:nCases
        for rec in 1:nRecons
            for it in 1:largestIt
                takeIt = min(it, length(outputs[case][rec].qmaps))
                maps[case,rec,1,it,:,:] = outputs[case][rec].qmaps[takeIt].T₁ 
                maps[case,rec,2,it,:,:] = outputs[case][rec].qmaps[takeIt].T₂
                maps[case,rec,3,it,:,:] = abs.(outputs[case][rec].qmaps[takeIt].ρ)
            end
        end
    end
    meanoverdyn= mean(maps,dims=2);

    # spectral analysis
    kxr         = saParams["kxRange"]     # selects a range of kx-values for analysis of ky-spectra 
    centerrange = saParams["centerrange"] # allows to select a central region for analysis 
    centerMaps=maps[:,:,:,:,centerrange,centerrange]         # case dyn map iter x y
    centerMapsTweak = centerMaps

    # possibly remove "suspicious" values (mask and CSF)
    if saParams["tweakForSpectrum"]
        defval = [1.0, 0.05]  # default values for T1 and T2
        centerMapsTweak = similar(centerMaps)
        for c in 1:nCases; for d in 1:nRecons; for m in 1:2; for it = 1:largestIt;
            x      = centerMaps[c,d,m,it,:,:]
            mmmask = (x .> 0.001) .* (x .< (2.5 * defval[m]))
            erode(mmmask)
            for i in eachindex(x); 
                bbb = mmmask[i];
                xxx = x[i];
                yyy = bbb ? xxx : defval[m];
                x[i] = yyy; 
            end;
            centerMapsTweak[c,d,m,it,:,:] = x;
        end; end; end; end;
    end

    halfsize = size(centerMapsTweak,6)÷2; quartersize=halfsize÷2
    fmaps= zeros(ComplexF64,size(centerMapsTweak));
    fmaps = fft(centerMapsTweak,[5,6])
    fmaps = fftshift(fftshift(fmaps, 5),6)
    devfmaps = std(fmaps,dims=2)
    figure(); imshow((abs.(devfmaps[1,1,2,largestIt,:,:])), cmap="gray");
    figure(); imshow((abs.(devfmaps[1,1,2,largestIt,1:quartersize,(halfsize+1):end])), cmap="gray");
    #ppp = mean((abs.(devfmaps[:,1,:,:,1:quartersize,:]).^2),dims=4)
    ppp = mean((abs.(devfmaps[:,1,:,:,kxr,:]).^2),dims=4)
    figure(); plot(sqrt.(ppp[1,2,largestIt,1,(halfsize+1):end]));
    #
    for c in 1:nCases
        mt1 = mean(meanoverdyn[c,1,1,7,centerrange,centerrange]);
        mt2 = mean(meanoverdyn[c,1,2,7,centerrange,centerrange])
        @show casenames[c], mt1, mt2
    end
    for c in 1:nCases
        mt1 = mean(meanoverdyn[c,1,1,largestIt,centerrange,centerrange]);
        mt2 = mean(meanoverdyn[c,1,2,largestIt,centerrange,centerrange])
        @show casenames[c], mt1, mt2
    end

    # Todo: read-in BLAKJac-predicted spectra 
    blj = saParams["compareBLAKJac"]
    if blj
        tag = saParams["tag"]
        fn = "/home/mfuderer/Documents/Julia/Capture/IEEETMI_BLJspectra$tag.jld2"
        vars = FileIO.load(fn)
        savedSpectra = vars["savedSpectra"] #cases maptype ky
        #bljSpectra = reverse(savedSpectra, dims=1) # now assumed to run as [2DW, 2DB, 0DW, 0DB]
        bljSpectra = savedSpectra
    end
    bljOrder = [2,1,4,3]

    fudgeScaleFactor = saParams["scaleFactors"]
    plotEnd = Int(floor((1.0+0.8)*halfsize))                               
    (fig,ax)=(subplots(1,2,figsize=(6,4)))
    for (i,c) in enumerate(cset)   
        for m in 1:2
            x = fudgeScaleFactor[m] .*sqrt.(ppp[c,m,largestIt,1,(halfsize+1):plotEnd])
            ax[m].plot( x , label=casenames[c], color=colortag[c])     
            if blj
                ax[m].plot( sqrt.(abs.(bljSpectra[bljOrder[c],m,:])), color=colortag[c], "--")
            end
        end
    end
    for i in 1:2
        ax[i].set_ylim(0,40)
        ax[i].set(title="T$i noise spectrum")
        ax[i].set(xlabel="spatial frequency [FOV⁻¹]");
        ax[i].legend()
    end
    ax[1].set(ylabel="root of noise power density [a.u.]")
end



function VolunteerScanAnalysis(outputs)
    nCases = saParams["nCases"]
    nRecons = saParams["nRecons"]
    latestIts = [length(outputs[case][i].qmaps) for i in 1:nRecons for case in 1:nCases]
    maxIt = minimum(latestIts) 
    largestIt = maximum(latestIts)
    casenames = saParams["casenames"]
    colortag = saParams["colortag"] 
    cset     = saParams["displSet"] 
    mapnames = ["T1", "T2", "rho"]
    cmaps    = ["MRF_T1","MRF_T2"]
    dispMax  = saParams["displMax"]

    # Reshape data into 6-d array
    (sx,sy) = size(outputs[1][1].qmaps[1].T₁)
    maps = zeros(nCases,nRecons,3,largestIt,sx,sy);

    for case in 1:nCases
        for rec in 1:nRecons
            for it in 1:largestIt
                takeIt = min(it, length(outputs[case][rec].qmaps))
                maps[case,rec,1,it,:,:] = outputs[case][rec].qmaps[takeIt].T₁ 
                maps[case,rec,2,it,:,:] = outputs[case][rec].qmaps[takeIt].T₂
                maps[case,rec,3,it,:,:] = abs.(outputs[case][rec].qmaps[takeIt].ρ)
            end
        end
    end
    maps = reverse(maps, dims=5)
    meanoverdyn= mean(maps,dims=2);
    devoverdyn = std(maps,dims=2); # case [dummy] map iter x y

    # display selected ROI
    roi=saParams["ROI"];

    # trend over iterations
    tttt = copy(meanoverdyn[1,1,1,largestIt,:,:]); tttt[roi...] .+=0.6;
    figure(); imshow(tttt, vmax=2.0);
    for map in 1:2
        figure(); 
        for i in cset
            mmm = [mean(devoverdyn[i,1,map,it,roi...]) for it in 1:largestIt];
            plot(mmm, label="$(casenames[i])", color=colortag[i]); 
        end; 
        xlabel("iterations")
        legend()
        title("deviation of $(mapnames[map]) in ROI");
    end
    for map in 1:2
        figure(); 
        for i in cset
            mmm = [mean(meanoverdyn[i,1,map,it,roi...]) for it in 1:largestIt];
            plot(mmm, label="$(casenames[i])", color=colortag[i]); 
        end; 
        xlabel("iterations")
        legend()
        title("mean of $(mapnames[map]) of ROI");
    end

    # trend over dynamics (added 2022-06-28)
    for map in 1:2
        figure(); 
        for i in cset
            mmm = [mean(maps[i,dyn,map,largestIt,roi...]) for dyn in 1:nRecons];
            plot(mmm, label="$(casenames[i])", color=colortag[i]); 
            @show casenames[i], mean(mmm), std(mmm)
        end; 
        xlabel("dynamic")
        legend()
        title("mean of $(mapnames[map]) of ROI");
    end    



    # Display of resulting images (mean)
    for m in 1:2
        (fig,ax)=(subplots(2,2,figsize=(9,9)))
        subplots_adjust(wspace=0.01,hspace=0.01)
        for i in eachindex(cset)
            pcm = ax[i].imshow(meanoverdyn[cset[i],1,m,largestIt,:,:],vmax=dispMax[m], cmap=cmaps[m]);
            ax[i].set_xticks([]); 
            ax[i].set_yticks([]);
            ax[i].set_yticklabels([]);
            ax[i].set_title(casenames[cset[i]]) # casenames[cset[i]]
            @show i, casenames[cset[i]];
            if i==length(cset); colorbar(pcm,ax=ax[:], shrink=0.8); end; # orientation="horizontal", 
        end
        fig.suptitle(mapnames[m])
    end

    # Save images for "manual" assembly
    tag = saParams["tag"];
    roi_x = saParams["roi_x"]
    fn = "/home/mfuderer/Documents/Julia/Capture/IEEETMIimages$tag.jld2"
    save_im = zeros(4,2,224,length(roi_x))
    for m in 1:2
        for (i,c) in enumerate(cset)   
            save_im[i,m,:,:] = meanoverdyn[c,1,m,largestIt,:,roi_x];
        end
    end
    FileIO.save(fn,"save_im",save_im)    

    # barplot of deviations
    devRoi = zeros(2,nCases)
    devAll = zeros(2,nCases)
    devSave = zeros(2,2,nCases)
    for m in 1:2
        for (i,c) in enumerate(cset)   
        #for i in 1:nCases
            devAll[m,c] = mean(devoverdyn[c,1,m,largestIt,:,:]);
            devRoi[m,c] = mean(devoverdyn[c,1,m,largestIt,roi...])
            devSave[1,m,i]=devAll[m,c]
            devSave[2,m,i]=devRoi[m,c]
        end; 
    end

    @show devAll
    @show devRoi
    #FileIO.save(fn,"devSave",devSave)    

    nnn = casenames[cset] |> vec
    ttt = [Float64(i) for i in 0:(length(cset)-1)]
    (fig,ax)=(subplots(1,2,figsize=(10,4)))
    for m in 1:2
        ax[m].set_xticks(ttt)
        ax[m].set_xticklabels(nnn, rotation=25)
        ax[m].plot(1000.0 .*devAll[m,cset],label="mean noise over image","x")
        ax[m].plot(1000.0 .*devRoi[m,cset],label="mean noise over WM ROI","x")
        ax[m].set_ylabel("standard deviation [ms]")
        ax[m].legend()
        ax[m].set_title(mapnames[m])
    end
end

function ExtractImage(outputs, type::Int, element::Int, iteration::Int)
    fullIm = similar(outputs[element].qmaps[iteration].T₁)
    if     type==1;   fullIm = outputs[element].qmaps[iteration].T₁;
    elseif type==2;   fullIm = outputs[element].qmaps[iteration].T₂;
    elseif type==3
        ρˣ = outputs[element].qmaps[iteration].ρˣ
        ρʸ = outputs[element].qmaps[iteration].ρʸ
        fullIm = abs.(ρˣ.+im.*ρʸ);
    elseif type==4;   fullIm = outputs[element].qmaps[iteration].B₁;
    else @assert(false);
    end
    return fullIm
end

function ZoomImage(fullIm, zoomfac::Float64)
    fs = size(fullIm,1); 
    firstPix = 1 + round(Int64, fs*0.5*(zoomfac-1.0)/zoomfac)
    lastPix  = round(Int64, fs*(1.0-0.5*(zoomfac-1.0)/zoomfac))
    cropIm = fullIm[firstPix:lastPix, firstPix:lastPix]

    x = range(1,size(fullIm,1),size(cropIm,1))
    y = range(1,size(fullIm,2),size(cropIm,2))
    itp = LinearInterpolation((x, y), cropIm)
    zoomedIm = [itp(y,x) for y in 1:size(fullIm,2), x in 1:size(fullIm,1)]
    return zoomedIm
end

function ExtractAndZoom(outputs, type::Int, element::Int, iteration::Int, zoomfac::Float64)
    fullIm = ExtractImage(outputs, type, element, iteration)
    zoomedIm = ZoomImage(fullIm, zoomfac)
    return zoomedIm
end

function ShowTubes(ρ, tubes)
    figure()
    subplot(121)
        # overlay of ROIs on proton density map
        ρ_max = maximum(ρ);
        for i in eachindex(tubes)
            ρ[tubes[i]] .= 2*ρ_max;
        end
        imshow(ρ)
        title("Mask per tube")

    subplot(122)
        # numbered
        x = zeros(size(ρ)...)
        for (i,r) in enumerate(tubes)
            x[r] .= i
        end
        imshow(x)
        colorbar()
    end

function EurospinScanAnalysis(outputs)
    nCases = saParams["nCases"]
    nRecons = saParams["nRecons"]
    latestIts = [length(outputs[case][i].qmaps) for i in 1:nRecons for case in 1:nCases]
    maxIt = minimum(latestIts) 
    largestIt = maximum(latestIts)
    casenames = saParams["casenames"]
    colortag = saParams["colortag"] 
    cset     = saParams["displSet"] 
    mapnames = ["T1", "T2", "rho"]
    cmaps    = ["MRF_T1","MRF_T2"]
    dispMax  = saParams["displMax"]

    # Reshape data into 6-d array
    (sx,sy) = size(outputs[1][1].qmaps[1].T₁)
    maps = zeros(nCases,nRecons,3,largestIt,sx,sy);

    for case in 1:nCases
        for rec in 1:nRecons
            for it in 1:largestIt
                takeIt = min(it, length(outputs[case][rec].qmaps))
                maps[case,rec,1,it,:,:] = outputs[case][rec].qmaps[takeIt].T₁ 
                maps[case,rec,2,it,:,:] = outputs[case][rec].qmaps[takeIt].T₂
                maps[case,rec,3,it,:,:] = abs.(outputs[case][rec].qmaps[takeIt].ρ)
            end
        end
    end
    maps = reverse(maps, dims=5)
    meanoverdyn= mean(maps,dims=2);
    devoverdyn = std(maps,dims=2); # case [dummy] map iter x y

    # detect tubes and number them
    ρ = maps[1,1,3,6,:,:]; 
    mask = ρ .> 0.1 * maximum(ρ);
    tubesize=65;
    tubes = _mask_per_tube(mask, tubesize);
    show_mask_results(ρ, tubes)
    gt1 = goldstandard_3T.T1
    gt2 = goldstandard_3T.T2
    tubenumbers = saParams["tubenumbers"]
    kiwi = 19

    # Per vial, determine std, mean and expected goldstandard
    iterchoices=[5, largestIt]; nit = length(iterchoices)
    snrs = zeros(nCases,2,nit,length(tubes));
    biases = zeros(nCases,2,nit,length(tubes));
    variances = zeros(nCases,2,nit,length(tubes));
    goldstandard = zeros(nCases,2,nit,length(tubes));
    means = zeros(nCases,2,nit,length(tubes));
    for (tnr, tube) in enumerate(tubes)
        if tubenumbers[tnr] < kiwi 
            tt = [gt1[tubenumbers[tnr]], gt2[tubenumbers[tnr]]];
            for map in 1:2
                for c in 1:nCases
                    for (i,it) in enumerate(iterchoices)
                        ddd = mean(devoverdyn[c,1,map,it,tube]);
                        mmm = mean(meanoverdyn[c,1,map,it,tube]);
                        ttt = tt[map]/1000;
                        snrs[c,map,i,tnr] = mmm/ddd;
                        biases[c,map,i,tnr] = abs(ttt-mmm)/ttt;
                        variances[c,map,i,tnr] = ddd^2;
                        goldstandard[c,map,i,tnr] = ttt;                
                        means[c,map,i,tnr] = mmm;                
                    end
                end; 
            end
        else
            # kiwi-related entries are purposely left on 0
            ;
        end 
    end

    # show resulting values
    meansnr = mean(snrs[:,:,:,:],dims=4);
    meanbias= mean(biases[:,:,:,:],dims=4);
    meanvar = mean(variances[:,:,:,:],dims=4)
    for m in 1:2
        for (i,it) in enumerate(iterchoices)
            @show it, mapnames[m], meansnr[cset,m,i,1];
            @show it, mapnames[m], meanbias[cset,m,i,1];
            @show it, mapnames[m], sqrt.(meanvar[cset,m,i,1])
        end
    end

    # plot correspondence of values to glod-standard
    (fig,ax)=(subplots(1,2,figsize=(11,6)))
    subplots_adjust(wspace=0.1,hspace=0.1)
    takeDiff = saParams["BlandAltmanStyle"]
    for map in 1:2; 
        for ccc in cset
            casename= casenames[ccc]
            y     = takeDiff ? means[ccc,map,1,:].-goldstandard[ccc,map,1,:] : means[ccc,map,1,:];
            ax[map].scatter(1000 .*goldstandard[ccc,map,1,:], 1000 .*y , label="$(casename)", color=colortag[ccc]);   
            x=goldstandard[ccc, map,1,:]; y=means[ccc,map,1,:]; mx=mean(x); my=mean(y); mn=mean((x.-mx).*(y.-my)); md=mean((x.-mx).^2); mi=my-mx*mn/md 
            leftY = takeDiff ? mi+dispMax[map]*(mn/md-1) : mi+dispMax[map]*mn/md ;
            ax[map].plot([0, 1000*dispMax[map]],[1000*mi, 1000*leftY], "--", color=colortag[ccc]);
        end
        leftRefY = takeDiff ? 0 : dispMax[map];
        ax[map].plot([0, 1000*dispMax[map]],[0, 1000*leftRefY]);
        ax[map].set_xlabel("gold standard [ms]", fontsize=18)
        ax[map].legend(fontsize=18);
    end;

    #ax[1].set_xticklabels(ax[1].get_xticklabels(),fontsize=18.0)
    #ax[1].set_yticklabels(ax[1].get_yticklabels(),fontsize=18.0)

    #ax[2].set_xticks([0.0, 100.0, 200.0, 300.0])
    #ax[2].set_xticklabels(["0", "100", "200", "300"],fontsize=18.0)
    ax[2].yaxis.set_label_position("right")
    ax[2].yaxis.tick_right()
    #ax[2].set_yticks([0.0, 100.0, 200.0, 300.0])
    #ax[2].set_yticklabels(["0", "100", "200", "300"],fontsize=18.0)
    label    = takeDiff ? "recon difference to gold standard [ms]" : "reconstructed value [ms]";
    ax[1].set_ylabel(label, fontsize=18); 
    ax[2].set_ylabel(label, fontsize=18); 

    # Graphs over iterations
    if saParams["iterGraphs"]
        for (tnr, tube) in enumerate(tubes)
            if tubenumbers[tnr] < kiwi 
                tt1 = gt1[tubenumbers[tnr]]
                tt2 = gt2[tubenumbers[tnr]]
                for map in 1:2
                    figure(); 
                    for i in cset
                        mmm = [mean(meanoverdyn[i,1,map,it,tube]) for it in 1:largestIt];
                        plot(mmm, label="$(casenames[i])"); 
                    end; 
                    xlabel("iterations")
                    legend()
                    title("mean of $(mapnames[map]) in ROI of tube ($tt1,$tt2)");
                end
            end
        end
    end

    # Display of resulting images (mean)
    for m in 1:2
        (fig,ax)=(subplots(2,2,figsize=(9,9)))
        subplots_adjust(wspace=0.01,hspace=0.01)
        for i in eachindex(cset)
            pcm = ax[i].imshow(meanoverdyn[cset[i],1,m,largestIt,:,:],vmax=dispMax[m], cmap=cmaps[m]);
            ax[i].set_xticks([]); 
            ax[i].set_yticks([]);
            ax[i].set_yticklabels([]);
            ax[i].set_title(casenames[cset[i]]) # casenames[cset[i]]
            @show i, casenames[cset[i]];
            if i==length(cset); colorbar(pcm,ax=ax[:], shrink=0.8); end; # orientation="horizontal", 
        end
        fig.suptitle(mapnames[m])
    end

    # barplots
    blakjacDevs = [3.2, 2.8, 5.2, 3.2]
    factor = [0.0125, 0.0022]
    width=0.35
    ct = [colortag[c] for c in cset]
    stdDev = sqrt.(meanvar[:,:,2,1])  # case map iterchoice dummy
    (fig,ax)=(subplots(1,2,figsize=(8,4)))
    fig.subplots_adjust(bottom=0.2)
    for m in 1:2
        sd = [stdDev[c,m] for c in cset]
        ax[m].set_xticks(1:length(cset))
        ax[m].set_xticklabels(casenames[cset],rotation=25, fontsize=14)
        ax[m].bar((1:4).+width/2,            1000.0 .*blakjacDevs.*factor[m], width=0.7*width, hatch="///"  ,color=ct)
        ax[m].bar((1:length(cset)).-width/2, 1000.0 .*sd,                     width,                         color=ct)
        ax[m].set_title(mapnames[m]);
    end
    ax[2].yaxis.set_label_position("right")
    ax[2].yaxis.tick_right()
    ax[1].set_ylabel("std.dev[ms]", fontsize=14); 
    ax[2].set_ylabel("std.dev[ms]", fontsize=14); 
end


function ReadInto7DArray(folderDetails, nZoom, nRecons, nTypes, maxIt, sx, sy, fn_common)
    maps = zeros(length(folderDetails), nZoom, nRecons, nTypes, maxIt, sx, sy);
    for fff in eachindex(folderDetails)
        fn = fn_common*folderDetails[fff]      
        outputs = get_old_recons(fn)
        for rec in 1:nRecons
            for it in 1:maxIt
                for type = 1:3
                    thisMap = ExtractAndZoom(outputs, type, rec, it, 1.0)
                    #@show fff, zoomtype, rec, type, it
                    maps[fff,zoomtype,rec,type,it,:,:] = thisMap
                end
            end
        end
    end
    maps=reverse(maps,dims=6);

    return maps;
end

function DefineTissueRegionsKnee(maps, maxIt, roi_marrow, roi_muscle)
    # get something like average-over-all T1, T2 and rho-derived mask
    mmm = mean(mean(mean(maps,dims=3),dims=2),dims=1);
    allMeansT1 = mmm[1,1,1,1,maxIt,:,:];
    allMeansT2 = mmm[1,1,1,2,maxIt,:,:];
    rhoMask    =1.0 *(abs.(mmm[1,1,1,3,maxIt,:,:])) .> 0.2;
    allidx=CartesianIndices(allMeansT1)

    # remove the skin components from mask
    width = 3.0
    sz = size(allMeansT1)[1]; centre = sz÷2+1;
    d2(x::CartesianIndex, y::CartesianIndex) = ( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    ggg   = exp.(-d2.((CartesianIndex(centre,centre),), allidx) / (2.0 *width^2)) ./(2*pi*width^2); #gaussian kernel image
    rhoSmooth = conv_psf(rhoMask,ggg);
    interior  = (rhoSmooth .> 0.7)

    # segment the principal tissues
    tiss    = [vec(allidx[roi_marrow...]) for i in 1:9] 
    # two lines commented out due to unintelligible error
    #tiss[1] =  vec(allidx[allMeansT1.<0.5 .&& allMeansT1.>0.3 .&& allMeansT2 .> 0.1 .&& allMeansT2 .<0.3 .&& interior ]) # marrow and fat
    #tiss[2] =  vec(allidx[allMeansT1.>1.2 .&& allMeansT1.<2.0 .&& allMeansT2 .> 0.01 .&& allMeansT2 .<0.045 .&& interior])   # muscle
    tiss[3] =  vec(allidx[roi_marrow...])
    tiss[4] =  vec(allidx[roi_muscle...])
    roi_musc_ant=[167:172,116:121]
    roi_musc_post=[167:172,164:169]
    tiss[5] =  vec(allidx[roi_musc_ant...])
    tiss[6] =  vec(allidx[roi_musc_post...])
    roi_marr_inf=[143:148,104:109]
    roi_marr_sup=[30:35,  100:105]
    tiss[7] =  vec(allidx[roi_marr_inf...])
    tiss[8] =  vec(allidx[roi_marr_sup...])
    roi_marr_small = [160:167,95:101] 
    tiss[9] =  vec(allidx[roi_marr_small...])

    roiNames = ["marrow", "muscle", "ROI in marrow", "ROI in muscle", "muscle anterior", "muscle posterior", "marrow inferior", "marrow superior",
                "roi_marr_small"]

    x = zeros(size(allMeansT1)...)
    for (i,r) in enumerate(tiss)
        x[r] .= i
    end
    figure(); imshow(x);

    return tiss, roiNames
end

function DefineTissueRegionsBrain(maps, maxIt, roi_wm)
    # get something like average-over-all T1, T2 and rho-derived mask
    mmm = mean(mean(mean(maps,dims=3),dims=2),dims=1);
    allMeansT1 = mmm[1,1,1,1,maxIt,:,:];
    allMeansT2 = mmm[1,1,1,2,maxIt,:,:];
    rhoMask    =1.0 *(abs.(mmm[1,1,1,3,maxIt,:,:])) .> 0.2;
    allidx=CartesianIndices(allMeansT1)

    # remove the skin components from mask
    width = 10.0
    sz = size(allMeansT1)[1]; centre = sz÷2+1;
    d2(x::CartesianIndex, y::CartesianIndex) = ( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    ggg   = exp.(-d2.((CartesianIndex(centre,centre),), allidx) / (2.0 *width^2)) ./(2*pi*width^2); #gaussian kernel image
    rhoSmooth = conv_psf(rhoMask,ggg);
    interior  = (rhoSmooth .> 0.7)

    # segment the principal three tissues
    tiss    = [allidx[allMeansT1.<1.15 .&& allMeansT1.>0.8 .&& allMeansT2 .> 0.03 .&& allMeansT2 .<0.1 .&& interior ] for i in 1:4]
    tiss[2] =  allidx[allMeansT1.>1.35 .&& allMeansT1.<1.5 .&& allMeansT2 .> 0.03 .&& allMeansT2 .<0.1 .&& interior]
    tiss[3] =  allidx[allMeansT2.>0.1.&& interior]
    tiss[4] =  vec(allidx[roi_wm...])

    roiNames = ["white matter", "gray matter", "CSF", "wmROI"]

    x = zeros(size(allMeansT1)...)
    for (i,r) in enumerate(tiss)
        x[r] .= i
    end
    figure(); imshow(x);

    return tiss, roiNames
end

function DefineTissueRegionsPhantom(meanoverdyn, maxIt)
    symbolicPos = [(0,1), (0,2), (1,0), (1,1), (1,2), (1,3), (2,0), (2,1), (2,2), (2,3), (3,1), (3,2), (2.6,2.6), (1.5,1.5)];
    origin = (49,39)
    distance = 44
    angulation = 0.14
    radius = 8
    centerPos = [origin.+(distance*x*cos(angulation)-distance*y*sin(angulation), distance*x*sin(angulation)+distance*y*cos(angulation)) for (x,y) in symbolicPos]
    tube_centers = [CartesianIndex(round(Int64,x),round(Int64,y)) for (x,y) in centerPos]
    
    ρ = meanoverdyn[1,1,1,3,maxIt,:,:]
    
    d(x::CartesianIndex, y::CartesianIndex) = √( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    allidx = CartesianIndices(ρ)
    
    tubes = [allidx[ d.((center,), allidx) .< radius] for center ∈ tube_centers];
    ShowTubes(ρ, tubes)
    return tubes
end

function ShowMeansDeviations(maps, description, tiss, nTypesShown, gels)
    for gel in gels
        tissue = tiss[gel]
        acc = 0.0
        for m in 1:nTypesShown
            for case in eachindex(description)
                md   = [mean(mean(maps[case,1,d,m,end,tissue]) ) for d in 1:size(maps)[3]]
                mstd  = mean(devoverdyn[case,1,1,m,end,tissue])
                re = 100.0 * std(md)/mean(md)
                acc += re*re
                varFromMean = [(maps[case,1,d,m,end,:,:].-meanoverdyn[case,1,1,m,end,:,:]) for d in 1:size(maps)[3]]
                relVarFromMean = [100.0 .*((maps[case,1,d,m,end,:,:]./(max.(meanoverdyn[case,1,1,m,end,:,:],0.01)).-1.0)) for d in 1:size(maps)[3]]
                if (case==1); 
                    figure(); imshow(relVarFromMean[5], vmin=-100, vmax=100, cmap="RdBu" ); colorbar(); 
                end
                stdRoi = [std(varFromMean[d][tissue]) for d in 1:size(maps)[3]]
                sss = description[case]
                @printf("Sequence %s in tissue %d for T%d: std_of_var=%.4f, mean(stddev)=%.4f, stddev(mean)=%.4f, mean=%.4f, rel.err=%.3f  \n",
                     sss, gel, m, mean(stdRoi), mstd[1], std(md), mean(md), re)
            end 
        end
        @printf "RMS of relative errors (in percent): %.1f \n " sqrt(acc/nTypesShown/length(description))
    end
end

function ShowStabilityOverDynamics(maps, description, tiss, nTypesShown, gels)
    for m in 1:nTypesShown
        figure()
        for case in 1:length(folderDetails)
            mdt  = [mean([mean(maps[case,1,d,m,end,tissue]) for tissue in tiss][gels]) for d in 1:size(maps)[3]]
            # @show mean(mdt), description[case]
            # t1vals[case] = mean(mdt2_1)
            plot(mdt, label=chop(description[case], head=1, tail=0))
        end
        legend()
    end
end

function ShowStabilityOverIterations(meanoverdyn, description, tiss, nTypesShown, roiNames, mapnames)
    for m in 1:nTypesShown
        gels = eachindex(tiss)
        nGels = length(gels)
        rows = Int64(round(sqrt(0.8*nGels)))
        cols = Int64(ceil(nGels/rows))
        (fig,ax) = subplots(rows,cols,figsize=(12,8))
        fig.subplots_adjust(wspace=0.15,hspace=0.25)
        for gel in gels
            row = (gel-1)÷cols +1
            col = (gel-1)%cols +1
            tissue = tiss[gel]
            ax[row,col].set_title("$(roiNames[gel]), $(mapnames[m])")
            for case in eachindex(description)
                meant2iter = [mean(meanoverdyn[case,1,1,m,it,tissue]) for it in 1:size(meanoverdyn)[5]]
                ax[row,col].plot(meant2iter, label=chop(description[case], head=1, tail=0));
            end
            if (gel==nGels) ax[row,col].legend(); end;
        end
    end
end

function ShowStabilityOverIterationsSingle(meanoverdyn, description, tissue, nTypesShown, roiName, mapnames)
    for m in 1:nTypesShown

        figure()
        title("$roiName, $(mapnames[m])")
        for case in eachindex(description)
            meant2iter = [mean(meanoverdyn[case,1,1,m,it,tissue]) for it in 1:size(meanoverdyn)[5]]
            plot(meant2iter, label=chop(description[case], head=1, tail=0));
        end
        legend()
    end
end

function DisplayImageSet(meanoverdyn, description, nTypesShown, nSubtypes)
    cmaps    = ["lipari","navia","gray"]
    mapnames = ["T1", "T2", "rho"]
    for m=1:nTypesShown
        (fig,ax) = subplots(2,nSubtypes,figsize=(16,7.3))
        #fig.subplots_adjust(wspace=0.01,hspace=0.01)
        fig.subplots_adjust(wspace=0.01,hspace=0.01,left=0.01,right=0.9,bottom=0.04,top=0.99)
        pcm = ax[1,nSubtypes]
        for case in eachindex(description)
            row = (case-1)÷nSubtypes +1;
            col = (case-1)%nSubtypes +1    
            thisMean = meanoverdyn[case,1,1,m,end,:,:] 
            pcm = ax[row,col].imshow(thisMean, vmax=dispMax[m], cmap=cmaps[m])
            ax[row,col].set_xticks([]); 
            ax[row,col].set_yticks([]);
            ax[row,col].set_yticklabels([]);
            ax[row,col].text(0.05, 0.8, description[case], fontsize=14, color="red", transform=ax[row,col].transAxes)
        end
    
        p0 = ax[1,nSubtypes].get_position().get_points()
        p1 = ax[2,nSubtypes].get_position().get_points()
        # p2 = ax[2].get_position().get_points().flatten()
        # left bottom width height
        @show p0, p1
        #ax_cbar = fig.add_axes([0.93, p1[1,2], 0.02, 0.7])
        ax_cbar = fig.add_axes([0.91, p1[1,2]+0.05, 0.02, 0.7])
        clb = plt.colorbar(pcm, cax=ax_cbar)
        clb.ax.tick_params(labelsize=20)
        fig.text(0.93,0.85,fontsize=20,mapnames[m]*"[s]",color="black")  # added 2023-01-17
    end
end

function DisplaySequences(fn_common, folder_deatils, nSubtypes)
    (fig,ax) = subplots(2,nSubtypes,figsize=(12,6))
    fig.subplots_adjust(wspace=0.01,hspace=0.01)
    for case in eachindex(folderDetails)
        row = (case-1)÷nSubtypes +1;
        col = (case-1)%nSubtypes +1    
        fn = fn_common*folderDetails[case]      
        outputs = get_old_recons(fn)
        RFdeg = outputs[1].sequence.RF_train

        ax[row,col].plot(abs.(RFdeg))
        anglesdd = zeros(length(RFdeg))
        for i in 1:length(RFdeg)-2
            anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
        end
        ax[row,col].plot(anglesdd)

        for i in 1:6; ax[row,col].plot([224*i-112, 224*i-112],[0,10],color="green"); end
        ax[row,col].text(0.05, 0.8, description[case], fontsize=14, color="red", transform=ax[row,col].transAxes)
    end
end



