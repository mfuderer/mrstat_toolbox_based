## Prepare for numerical phantom MR-STAT reconstruction

# Load external packages
using Revise
using JLD2, MAT, Serialization
using StaticArrays, StructArrays
using Random
using ComputationalResources, CUDA

# Load main MRSTAT package
using MRSTAT

# Select hardware resource
resource = has_cuda_gpu() ? CUDALibs() : CPU1()

# Code
include("../numerical_phantom/load_data_phantom.jl");
include("../numerical_phantom/RF_Shapes.jl");

# Numerical phantom setup
recon_options = load_default_recon_options();

recon_options["description"] = "";
recon_options["recon_folder"] = "tmp";

recon_options["numphantom"] = true;
recon_options["numphantom_type"] = "brainweb"; # "line", "brainweb", "checkerboard", "shepp_logan", "tubes"
recon_options["numphantom_size"] = [32,32];
recon_options["numphantom_sequence"] = "Spoiled"; # "Spoiled", "Balanced"
recon_options["numphantom_trajectory"] = "Cartesian"; # "Cartesian", "Radial"
recon_options["numphantom_rf_shape"] = "insync_varpeaks";
recon_options["numphantom_noise"] = false;

recon_options["trf_max_iter"] = 20;

# Set plotting function

# use pyplot to plot parameter maps at each iteration

# includet("../plot/pyplot.jl");
# pygui(true); # needed because of bug in vscode extension, will be fixed soonᵀᴹ
# plotfn = pyplot;

# on cluster without graphical user interface, plot in terminal directly

# includet("../plot/terminalplot.jl")
# @info "Set plotting function: ImageInTerminal backend"
# plotfn = terminalplot;

# or don't plot at all

plotfn(x; figtitle="") = println("no plotting")

## Run MR-STAT reconstruction

sequence, coordinates, coilmaps, trajectory, mask, phantom = load_data_phantom(recon_options);

output = mrstat_recon(sequence, trajectory, phantom, mask, coordinates, coilmaps, recon_options, resource, plotfn);

julia_to_matlab(output);