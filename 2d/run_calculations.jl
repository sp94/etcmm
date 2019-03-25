using JLD
using PyPlot
using PyCall
using ProgressMeter
using LinearAlgebra
using LorentzDrudeMetals
BLAS.set_num_threads(7)
include("FiniteElementHomogenisation.jl")
include("GmshTemplates.jl")
meshio = pyimport("meshio")

# decrease the resolution of wavelength scans, etc, by factor of sf
# e.g. for a quick test plot, use sf = 10
sf = 1


######################################
####### Define mixing formulas #######
######################################

function maxwellgarnett(epm,ep1,f)
    @assert epm == 1
    return sqrt.((ep1+1+f*ep1-f)/(ep1+1-f*ep1+f)+0im)
end

function maxwellgarnett_hex(epm,ep1,r,R)
    @assert epm == 1
    f = pi/2/sqrt(3) * (r/R)^2
    return maxwellgarnett(epm,ep1,f)
end

function maxwellgarnett_sqr(epm,ep1,r,R)
    @assert epm == 1
    f = pi/4 * (r/R)^2
    return maxwellgarnett(epm,ep1,f)
end

function rayleigh_sqr(epm,ep1,r,R)
    @assert epm == 1
    f = pi/4 * (r/R)^2
    return sqrt.(1 + 2*f/((ep1+1)/(ep1-1)-f-(ep1-1)/(ep1+1)*(0.3058*f^2+0.0134*f^8)))
end


######################################
############# Figure 2a ##############
######################################

# println("Begin Figure 2a; hexagonal lattice, d=38nm, L=40nm (G=2nm); both polarisations")
# jldopen("./data/figure2a.jld", "w") do f
#     lam0s = union(2:0.5*sf:10, 10:2*sf:50, 50:5*sf:100, 100:10*sf:240) # microns
#     neffs_TM = ComplexF64[]
#     neffs_TE = ComplexF64[]
#     ns_bulk = ComplexF64[]
#     @showprogress for lam0 in lam0s
#         k0 = 2pi/lam0/1000 # [1/nm]
#         ep = LorentzDrudeMetals.Ti[lam0]
#         n = sqrt(ep)
#         push!(ns_bulk, n)
#         delta_s = lam0 / (-4pi*imag(n)) * 1000 # skin depth [nm]
#         msh = make_msh_hexlat_cyl(38, 38, 40, 40*sqrt(3), delta_s, 20)
#         m = Mesh(msh, [1.0,ep], [1.0,1.0])
#         neff = homogenize(m, k0, :TM)[1][1]
#         push!(neffs_TM, neff)
#         neff = homogenize(m, k0, :TE)[1][1]
#         push!(neffs_TE, neff)
#     end
#     write(f, "lam0s", lam0s)
#     write(f, "neffs_TM", neffs_TM)
#     write(f, "neffs_TE", neffs_TE)
#     write(f, "ns_bulk", ns_bulk)
# end


######################################
############# Figure 2c ##############
######################################

#println("Begin Figure 2c; ")


######################################
############# Figure 4a ##############
######################################

# println("Begin Figure 4a; square lattice, d=38nm, L=40nm (G=2nm); different metals")
# jldopen("./data/figure4a.jld", "w") do f
#     lam0s = union(2:0.5*sf:10, 10:2*sf:50, 50:5*sf:100, 100:10*sf:240)
#     write(f, "lam0s", lam0s)
#     for (label,metal) in [("Al",LorentzDrudeMetals.Al),
#                           ("Ag",LorentzDrudeMetals.Ag),
#                           ("Au",LorentzDrudeMetals.Au),
#                           ("Ti",LorentzDrudeMetals.Ti)]
#         neffs = []
#         @showprogress for lam0 in lam0s
#             k0 = 2pi/lam0/1000 # [1/nm]
#             ep = metal[lam0]
#             n = sqrt(ep)
#             delta_s = lam0 / (-4pi*imag(n)) * 1000 # skin depth [nm]
#             @assert isapprox(delta_s, -1/2/imag(n)/k0)
#             msh = make_msh_sqrlat_cyl(38, 38, 40, 40, delta_s, 20)
#             m = Mesh(msh, [1.0,ep], [1.0,1.0])
#             neff = homogenize(m, k0, :TE)[1][1] ##################### TODO CHANGE THIS TO TM ONCE WE FIX LABELS
#             push!(neffs, neff)
#         end
#         write(f, "neffs_$(label)", neffs)
#     end
# end


######################################
############# Figure 4b ##############
######################################

println("Begin Figure 4b; square lattice, fixed wavelength; different filling fractions")
jldopen("./data/figure4b.jld", "w") do f
    lam0 = 200 # microns
    k0 = 2pi/lam0/1000 # [1/nm]
    write(f, "lam0", lam0)
    ep = LorentzDrudeMetals.Al[lam0]
    n = sqrt(ep)
    delta_s = lam0 / (-4pi*imag(n)) * 1000 # skin depth [nm]
    d_by_delta_s = 10 .^ range(-1.8, stop=1.8, length=15÷sf)
    write(f, "d_by_delta_s", d_by_delta_s)
    # MG and PC limits
    for (i,G_by_d) in enumerate([1/2,1/4,1/8])
        d = d_by_delta_s[end]*delta_s
        G = G_by_d*d
        @show d, d+G
        msh = make_msh_sqrlat_cylhole(d, d, d+G, d+G)
        m = Mesh(msh, [1], [1], mirror=true)
        neff_PC = homogenize(m, k0, :TE)[1][1] ##################### TODO CHANGE THIS TO TM ONCE WE FIX LABELS
        neff_MG = maxwellgarnett_sqr(1,ep,d,d+G)
        write(f, "neff_PC_G_by_d=$(G_by_d)", neff_PC)
        write(f, "neff_MG_G_by_d=$(G_by_d)", neff_MG)
    end
    # Numerical homogenisation
    for (i,G_by_d) in enumerate([1/2,1/4,1/8])
        neffs = ComplexF64[]
        for d in d_by_delta_s.*delta_s
            G = G_by_d*d
            msh = make_msh_sqrlat_cyl(d, d, d+G, d+G, delta_s, 20)
            @show d, d+G
            m = Mesh(msh, [1,ep], [1,1])
            @time neff = homogenize(m, k0, :TE)[1][1] ##################### TODO CHANGE THIS TO TM ONCE WE FIX LABELS
            push!(neffs, neff)
        end
        write(f, "neffs_G_by_d=$(G_by_d)", neffs)
    end
end


######################################
############# Figure 5a ##############
######################################

println("Begin Figure 5a; square lattice, fixed gap size; changing aspect ratio")
jldopen("./data/figure5a.jld", "w") do f
    lam0 = 200 # microns
    k0 = 2pi/lam0/1000 # [1/nm]
    ep = LorentzDrudeMetals.Au[lam0]
    n = sqrt(ep)
    delta_s = lam0 / (-4pi*imag(n)) * 1000 # skin depth [nm]
    d = 30 # nm
    aspects = 2 .^ range(-3, stop=3, length=21÷sf)
    write(f, "aspects", aspects)
    for G_by_d in [1/2, 1/4, 1/8]
        neffs = ComplexF64[]
        neffs_MG = ComplexF64[]
        @showprogress for aspect in aspects
            lx = d*sqrt(aspect)
            ly = d/sqrt(aspect)
            Lx = lx + G_by_d*d
            Ly = ly + G_by_d*d
            msh = make_msh_sqrlat_cyl(lx, ly, Lx, Ly, delta_s, 20)
            m = Mesh(msh, [1,ep], [1,1])
            neff = homogenize(m, k0, :TE)[1][1]
            push!(neffs, neff)
            frac = pi/4*lx*ly/Lx/Ly
            # MG for elliptical cylinders (see supplementary material)
            neff_MG = sqrt((ep+aspect+aspect*frac*(ep-1))/(ep+aspect-frac*(ep-1)))
            push!(neffs_MG, neff_MG)
        end
        write(f, "neffs_G_by_d=$(G_by_d)", neffs)
        write(f, "neffs_MG_G_by_d=$(G_by_d)", neffs_MG)
    end
end


######################################
############# Figure 5b ##############
######################################

println("Begin Figure 5b; hexagonal lattice, fixed lattice constant; changing diameter")
jldopen("./data/figure5b.jld", "w") do f
    lam0 = 2 # microns
    k0 = 2pi/lam0/1000 # [1/nm]
    ep = LorentzDrudeMetals.Au[lam0]
    n = sqrt(ep)
    delta_s = lam0 / (-4pi*imag(n)) * 1000 # skin depth [nm]
    ds = 0:sf:48
    neffs = ComplexF64[]
    neffs_MG = ComplexF64[]
    @showprogress for d in ds
        if d == 0
            push!(neffs, 1)
            push!(neffs_MG, 1)
            continue
        end
        msh = make_msh_hexlat_cyl(d, d, 50, 50*sqrt(3), delta_s, 20)
        m = Mesh(msh, [1.0,ep], [1.0,1.0])
        neff = homogenize(m, k0, :TE)[1][1]
        push!(neffs, neff)
        neff_MG = maxwellgarnett_hex(1,ep,d,50)
        push!(neffs_MG, neff_MG)
    end
    write(f, "ds", ds)
    write(f, "neffs", neffs)
    write(f, "neffs_MG", neffs_MG)
    write(f, "delta_s", delta_s)
end
