# Collection of various (ky,kz)-pattern definitions

function k_ky_inner()
    nky = recon_options["nky"]
    nkz = recon_options["nkz"]
    nTR = recon_options["nTR"]
    ky = 1.0 .*(repeat(1:nky, inner=1, outer= Int64(ceil(nTR/nky)) ))      
    kz = 1.0 .*(repeat(1:nkz, inner=nky, outer= Int64(ceil(nTR/nky/nkz)))) 
    return ky[1:nTR],kz[1:nTR]
end

function k_yzzy()
    nky = recon_options["nky"]
    nkz = recon_options["nkz"]
    nTR = recon_options["nTR"]
    ky_yz = 1.0 .*(repeat(1:nky, inner=1, outer= Int64(ceil(nTR/nky)) ))      
    kz_yz = 1.0 .*(repeat(1:nkz, inner=nky, outer= Int64(ceil(nTR/nky/nkz)))) 
    ky_zy = 1.0 .*(repeat(1:nky, inner=nkz, outer= Int64(ceil(nTR/nky/nkz)) ))      
    kz_zy = 1.0 .*(repeat(1:nkz, inner=1, outer= Int64(ceil(nTR/nkz)))) 
    switch = 3*nky*nkz
    ky[1:switch]     = ky_yz[1:switch]
    ky[switch+1:nTR] = ky_zy[switch+1:nTR]
    kz[1:switch]     = kz_yz[1:switch] 
    kz[switch+1:nTR] = kz_zy[switch+1:nTR]
    return ky, kz
end

function kz_as_trajectories()
    nky = recon_options["nky"]
    nkz = recon_options["nkz"]
    nTR = recon_options["nTR"]
    ky = 1.0 .*(repeat(1:nky, inner=1, outer= Int64(ceil(nTR/nky)) ))      
   
    trajectorySet = Vector{Vector{TrajectoryElement}}(undef,length(ky))
    for i in 1:length(ky)
        tr = [TrajectoryElement(ky[i]-nky÷2-1, Float64(kz-nkz÷2-1)) for kz in 1:nkz]
        trajectorySet[i] = Vector{TrajectoryElement}(tr)

    end
    return trajectorySet
end

function spirals_just_as_a_try(turns::Float64)
    nky = recon_options["nky"]
    nkz = recon_options["nkz"]
    nTR = recon_options["nTR"]
    oversample = 2.0
    r = nky/2.0
    nSamples = Int(ceil(pi*turns* r * oversample));
    q = sqrt(1.0/nSamples)

    angle = 0
    ga = pi*(3.0-sqrt(5.0));

    trajectorySet = Vector{Vector{TrajectoryElement}}(undef,nTR)
    tr            = Vector{TrajectoryElement}(undef,nSamples)
    for TR in 1:nTR
        angle += ga
        for i in 1:nSamples
            u = q*(sqrt(i)-1)  # approximation of an approximation ...
            θ = 2*pi*turns * u + angle;
            tr[i] = TrajectoryElement(u*r*cos(θ), u*r*sin(θ))
        end
        trajectorySet[TR] = copy(tr);
    end
    return trajectorySet
end

function single_voxel()
    nTR = recon_options["nTR"]

    trajectorySet = Vector{Vector{TrajectoryElement}}(undef,nTR)
    tr            = Vector{TrajectoryElement}(undef,1)
    for TR in 1:nTR
        tr=[TrajectoryElement(0,0)]
        trajectorySet[TR] = copy(tr);
    end
    return trajectorySet
end