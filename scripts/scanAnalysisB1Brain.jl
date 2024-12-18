## Created 2024-05-06 for Julia-calls to FSL functions BET and FAST 
#  Style attempts to follow figuresB1phantom.jl
using ImageMorphology
using ImageFiltering

function BrainSegment(case, r, figurePars, label=2, threshold=0.7)
    meanoverdyn   = figurePars["meanoverdyn"]   
    sliceRange    = figurePars["sliceRange"]   
    testIt        = figurePars["testIt"]
    t1map = meanoverdyn[r][case,sliceRange,1,testIt[r,case],:,:] 
    t2map = meanoverdyn[r][case,sliceRange,2,testIt[r,case],:,:] 
    rhomap = meanoverdyn[r][case,sliceRange,3,testIt[r,case],:,:]
    t1wemul = t1map .* rhomap ./ (t1map.^2 .+ 0.01^2) 
    t1wemul = permutedims(t1wemul,(3,2,1))
    t1wemul = reverse(t1wemul,dims=2)
    #t1wemul = reshape(t1wemul,size(t1wemul)[1],size(t1wemul)[2],size(t1wemul)[3])
    
    # Get header structure from similar task and fill in appropriate dimensions
    fnsample = "/smb/user/mfuderer/BLD_RT_RESEARCH_DATA/PROJECT/MRSTAT/experiments/Jordi/Studies/Demonstrator/EPT/All_seperated/MRSTAT01/MRSTAT01_Conv_FLAIRw.nii.gz"
    ni = niread(fnsample)
    ni_ext = ni.extensions
    ni_h   = ni.header
    t1wemul_dims = size(t1wemul)
    ni_h.dim = (t1wemul_dims..., ones(Int64, 8-length(t1wemul_dims))...)
    
    fsl_fn    = "/home/mfuderer/tmp/nisti_tmp.nii"
    fslbet_fn = "/home/mfuderer/tmp/nistibet_tmp.nii"
    
    # Cleanup tmp directory before running FSL
    tmp_dir = "/home/mfuderer/tmp"
    for file in readdir(tmp_dir)
        rm(joinpath(tmp_dir, file))
    end

    # write t1wemul to file and run FSL BET and FAST
    ni_out = NIVolume(ni_h,ni_ext,t1wemul)
    niwrite(fsl_fn, ni_out)
    shellCall = `/usr/local/fsl/share/fsl/bin/bet $fsl_fn $fslbet_fn`
    run(shellCall)
    shellCall = `/usr/local/fsl/share/fsl/bin/fast $fslbet_fn.gz`
    run(shellCall)
    
    # labels = niread("/home/mfuderer/tmp/nistibet_tmp_seg.nii.gz")
    # figure(); imshow(labels[:,:,10], cmap="gray", vmin=0, vmax=5)
    wmfraction = niread("/home/mfuderer/tmp/nistibet_tmp_pve_$label.nii.gz")
    wmfraction = reverse(wmfraction,dims=2)
    wmfraction = permutedims(wmfraction, (3,2,1))
    allidx = CartesianIndices(wmfraction)
    mask = (wmfraction[allidx] .> threshold)

    mask = erode(mask)
    
    wmroi  = allidx[mask]

    figurePars["wmroi"] = wmroi
    return mean(t2map[wmroi]), std(t2map[wmroi]);
end

# Smoothing of a volume, weighted with a mask
function SmoothOverMask(case, r, figurePars, type)
    meanoverdyn   = figurePars["meanoverdyn"]   
    sliceRange    = figurePars["sliceRange"]   
    testIt        = figurePars["testIt"]
    map           = meanoverdyn[r][case,sliceRange,type,testIt[r,case],:,:] 
    wmroi         = figurePars["wmroi"]

    justRoi   = zeros(size(map))
    mapRoi    = zeros(size(map)) # ; mapRoi .= mean(map[wmroi])
    smoothRoi = copy(mapRoi)
    for p in wmroi
        justRoi[p] = 1.0
        mapRoi[p]  = map[p]  
    end 

    kernel = Kernel.gaussian((1,1,1))
    convolvedMapRoi  = imfilter(mapRoi, kernel)
    convolvedJustRoi = imfilter(justRoi, kernel)

    # Normalize the convolution for points within ROI (the rest stays at the mean value)
    for p in wmroi
        smoothRoi[p] = convolvedMapRoi[p] / convolvedJustRoi[p]
    end 
    return smoothRoi
end
