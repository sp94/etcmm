using JLD
using LinearAlgebra
using LorentzDrudeMetals
using ProgressMeter
using Random
BLAS.set_num_threads(7)
include("FiniteDifferenceHomogenisation.jl")

# decrease the resolution of wavelength scans, etc, by factor of sf
# e.g. for a quick test plot, use sf = 10
sf = 1


######################################
############# Figure 2c ##############
######################################

println("Begin Figure 2c; square cylinders on a square lattice")
begin
    res = 0.2 # grid resolution in nanometers
    lam0s = union(2:1*sf:10, 10:10*sf:240) # microns
    metal = LorentzDrudeMetals.Ti
    dname = joinpath("./data/figure2c/")
    @showprogress for lam0 in shuffle(lam0s)
        k0 = 2pi/lam0/1000 # [1/nm]
        fname = joinpath(dname, "lam0=$(lam0).jld")
        if isfile(fname)
            continue
        end
        # Create cubic cylinders by first creating cubes
        ep, mu, dx, dy, dz = makecubes([6.5,1], [metal[lam0],1], res)
        # And then taking a cross-section through the middle
        mid = size(ep,2)÷2
        ep = ep[:,mid:mid,:,:]
        mu = mu[:,mid:mid,:,:]
        neffs = getneff(ep, mu, dx, dy, dz, k0)
        save(fname, "neffs", neffs, "lam0", lam0, "dx", dx, "dy", dy, "dz", dz)
    end
end


######################################
############# Figure 2b ##############
######################################

println("Begin Figure 2b; spheres on a hexagonal lattice, d=20nm, L=22nm (G=2nm)")
println("Note that a resolution of dx=0.75nm was used in the published work")
begin
    #res = 0.75 # grid resolution in nanometers
    res = 1.5; println("Here we have reduced the resolution to dx=1.5nm so the code can be run in about an hour, but the result will not be well converged")
    lam0s = union(2:0.25*sf:10, 10:1*sf:20, 20:10*sf:240) # microns
    metal = LorentzDrudeMetals.Ti
    dname = joinpath("./data/figure2b/")
    @showprogress for lam0 in shuffle(lam0s)
        k0 = 2pi/lam0/1000 # [1/nm]
        fname = joinpath(dname, "lam0=$(lam0).jld")
        if isfile(fname)
            continue
        end
        ep, mu, dx, dy, dz = makesphereshex([10,1], [metal[lam0],1], res)
        neffs = getneff(ep, mu, dx, dy, dz, k0)
        save(fname, "neffs", neffs, "lam0", lam0, "dx", dx, "dy", dy, "dz", dz)
    end
end
