@info "Running mrstat_toolbox_based startup.jl"

# Load external packages

    using Revise
    using JLD2, MAT, Serialization
    using StaticArrays, StructArrays
    using Random
    using FFTW, FourierTools
    using Interpolations
    using Printf


# Include sources on Lipari and Navia color maps 
include("/home/mfuderer/colorResources/RelaxationColor.jl") 



# # "Auto-select available hardware resource (single CPU, multi-CPU or GPU)"

#     resource = has_cuda_gpu() ? CUDALibs() : CPU1()

# # Plotting functions  
#     using Colors
#     using PyPlot

#     includet("plot/pyplot.jl");
#     @info "Set plotting function: PyPlot backend"
#     qplot = pyplot;



