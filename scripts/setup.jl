## Define packages
using LinearAlgebra # package for linear algebra operations (dot, )
using DelimitedFiles # For input of .txt files containing RF-patterns
using FFTW
using FileIO
using ComputationalResources, CUDA

# Before adding PyPlot, do
#     ENV["PYTHON"] = "" in the Julia REPL.
#     Add PyPlot in Pkg mode
#     Pkg.build("PyCall")
#     In vscode, deactivate settings > Extensions > Julia > Use Plot Plane
#     Restart Julia
using PyPlot # same plotting backend as commonly used in python
using PyCall

cpu  = ComputationalResources.CPU1()
using BlochSimulators
using Colors
using Statistics
using BLAKJac
using MRSTATToolbox
using Optim
include("recon_options.jl")

include("../numerical_phantom/RF_Shapes.jl")
include("../numerical_phantom/k_shapes.jl")
# src/MRSTAT/src/Utils/RF_Shapes.jl
if (!@isdefined(gelphantoms)) include("gelphantoms.jl") end
if (!@isdefined(tissues)) include("tissues.jl") end
@pyimport matplotlib.animation as anim
