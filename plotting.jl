using PyPlot
using PyCall
using JLD
using LorentzDrudeMetals
patches = pyimport("matplotlib.patches")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")
image = pyimport("matplotlib.image")
gs = pyimport("matplotlib.gridspec")
# https://github.com/JuliaPy/PyPlot.jl/issues/149
slice(x...) = pycall(pybuiltin("slice"), PyObject, x...)

######################################
########## Configure plots ###########
######################################

rc("font", size=8)
rc("figure", dpi=300)
# Width of one column figure
ONECOL = 88 / 25.4 # 88mm in inches
# Width of two column figure
TWOCOL = 180 / 25.4 # 180mm in inches
rc("xtick", direction="in")
rc("ytick", direction="in")


######################################
############## Figure 2 ##############
######################################

# CONFIGURE GRIDSPEC

fig = figure(figsize=(ONECOL,0.75ONECOL))
gs0 = gs.GridSpec(2, 2, figure=fig, width_ratios=[7,4], height_ratios=[4,3])

a = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(slice(0,2),0)), hspace=0.05)
a1 = subplot(get(a,(0,0)))
title("a", loc="left", ha="right", weight="bold")
a2 = subplot(get(a,(1,0)))

b = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(0,1)), hspace=0.1)
b1 = subplot(get(b,(0,0)))
title("b", loc="left", ha="right", weight="bold")
b2 = subplot(get(b,(1,0)))

c = subplot(get(gs0,(1,1)))
title("c", loc="left", ha="right", weight="bold")

# Figure 2a
data = load("./fem/data/figure2a.jld")
lam0s = data["lam0s"]
neffs_TE = data["neffs_TE"]
neffs_TM = data["neffs_TM"]
ns_bulk = data["ns_bulk"]
subplot(a1)
plot(lam0s, +real(ns_bulk.^2)/1e3, color="grey", ls="--")
plot(lam0s, +real(neffs_TE.^2)/1e3, color="#be2b2b")
plot(lam0s, +real(neffs_TM.^2)/1e3, color="#2bbebe")
text(240,-0.05,L"\hat{\epsilon}_\mathrm{eff}^\mathrm{TE}\approx +\! 10", ha="right", va="top", color="#be2b2b", fontsize=8)
text(240,-0.80,L"\hat{\epsilon}_\mathrm{eff}^\mathrm{TM}\ll0", ha="right", va="bottom", color="#2bbebe", fontsize=8)
ylim(-1.2,0.2)
xticks([])
yticks([0,-1],[0,L"-10^3"])
xticks([2,100,240], ["","","",""])
subplot(a2)
semilogy(lam0s, -imag(ns_bulk.^2), color="grey", ls="--")
semilogy(lam0s, -imag(neffs_TE.^2), color="#be2b2b")
semilogy(lam0s, -imag(neffs_TM.^2), color="#2bbebe")
xlabel("Wavelength [μm]")
xticks([2,100,240])
ylim(0.2e-2, 5e5)
yticks([1e-2,1e5])

# Figure 2b
lam0s = Float64[]
neffs = ComplexF64[]
dname = "./fdfd/data/figure2b/"
for fname in readdir(dname)
    data = load(joinpath(dname,fname))
    lam0 = data["lam0"]
    neff_lo, neff_hi = data["neffs"][1:2]
    push!(lam0s, lam0)
    push!(neffs, (neff_lo+neff_hi)/2)
end
idx = sortperm(lam0s)
lam0s = lam0s[idx]
neffs = neffs[idx]
subplot(b1)
plot(lam0s, +real(neffs.^2), "-", color="#be2b2b")
ylabel(L"\mathrm{Re}[\hat{\epsilon}_\mathrm{eff}]")
xticks([2,100,240], ["","","",""])
subplot(b2)
semilogy(lam0s, -imag(neffs.^2), color="#be2b2b")
xlabel("Wavelength [μm]")
ylabel(L"\mathrm{Im}[\hat{\epsilon}_\mathrm{eff}]")
xticks([2,100,240])
ylim(10^-3.5, 10^0.5)
yticks([1e-3, 1e0])

# Figure 2c
subplot(c)
# # Germanium
# d_ge = lam0s_ge*1e-6./ks_ge/4/pi
# semilogy(lam0s_ge, d_ge, "C2--")
xlabel("Wavelength [μm]")
ylabel("Penetration\nlength [m]")
xlim(-10,250)
ylim(1e-4,10^-0.75)
xticks([2,100,240])
# text(240,10^(-2.5),"Ge", ha="right", va="bottom", color="C2", fontsize=8)
# Titanium cylinder
dname = "./fdfd/data/figure2c/"
lam0s = []
neffs = []
for fname in readdir(dname)
    data = load(joinpath(dname,fname))
    lam0 = data["lam0"]
    neff = data["neffs"][1]
    push!(lam0s, lam0)
    push!(neffs, neff)
end
idx = sortperm(lam0s)
lam0s = lam0s[idx]
neffs = neffs[idx]
ks = -imag(neffs)
ds = lam0s*1e-6./(4pi*ks)
semilogy(lam0s, ds, "-", color="#be2b2b")
# text(240,1e-1, "Ti", ha="right", va="top", color="#be2b2b", fontsize=8)
yticks([1e-4,1e-1])

fig.align_ylabels([a1,a2])
fig.align_ylabels([b1,b2,c])

tight_layout()

savefig("./figure2.pdf"); println("Figure 2 saved")


######################################
############## Figure 4 ##############
######################################

fig = figure(figsize=(ONECOL,1.75ONECOL))
gs0 = gs.GridSpec(2, 1)

a = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(0,slice(0,2))), hspace=0)
a1 = subplot(get(a,(0,0)))
title("a", loc="left", ha="right", weight="bold")
title("Cylinders, "*L"d=38\,\mathrm{nm}, G=2\,\mathrm{nm}", fontsize="small")
ylabel(L"\mathrm{Re}[n_\mathrm{eff}]")
xticks([])
ylim(2.68,2.9)
yticks([2.7,2.8,2.9])
a2 = subplot(get(a,(1,0)))
xlabel("Wavelength [μm]")
ylabel(L"\mathrm{Im}[n_\mathrm{eff}]")
xticks([2,100,240])
ylim(0,0.11)
yticks([0,0.05,0.1])
data = load("./fem/data/figure4a.jld")
lam0s = data["lam0s"]
for label in ["Al", "Ag", "Au", "Ti"]
    neffs = data["neffs_$(label)"]
    a1.plot(lam0s, +real(neffs), label=label)
    a2.plot(lam0s, -imag(neffs), label=label)
end
subplot(a1)
legend(loc="upper right", ncol=4, handlelength=1.0)


b = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(1,slice(0,2))), hspace=0.05)
b1 = subplot(get(b,(0,0)))
b2 = subplot(get(b,(1,0)))
data = load("./fem/data/figure4b.jld")
d_by_delta_s = data["d_by_delta_s"]
for (i,G_by_d) in enumerate([1/2, 1/4, 1/8])
    neffs = data["neffs_G_by_d=$(G_by_d)"]
    neff_PC = data["neff_PC_G_by_d=$(G_by_d)"]
    neff_MG = data["neff_MG_G_by_d=$(G_by_d)"]
    subplot(b1)
    semilogx(d_by_delta_s, +real(neffs), color="C$(i+2)")
    semilogx([10^-2.5, 10^-1.8], real([neff_MG, neff_MG]), ls="--", color="C$(i+2)")
    semilogx([10^1.8, 10^2.5], real([neff_PC, neff_PC]), ls="--", color="C$(i+2)")
    subplot(b2)
    loglog(d_by_delta_s, -imag(neffs), color="C$(i+2)")
    loglog([10^-2.5, 10^-1.8], -imag([neff_MG, neff_MG]), ls="--", color="C$(i+2)")
end
subplot(b1)
title("b", loc="left", ha="right", weight="bold")
xlim(10^-2.5, 10^2.5)
ylim(1, 2.5)
text(10^-2.4, 2.1, "MG")
text(10^+2.35, 1.4, "PC", ha="right")
xticks([])
yticks(1:0.5:2.5)
ylabel(L"\mathrm{Re}[n_\mathrm{eff}]")
subplot(b2)
xlim(10^-2.5, 10^2.5)
ylim(10^-6.0, 10^0)
yticks([10^-5.0, 10^-3.0, 10^-1.0])
text(10^-2.4, 10^-4.75, "MG")
xlabel(L"d/\delta_s")
ylabel(L"\mathrm{Im}[n_\mathrm{eff}]")
plot([], [], label=L"G=d/8", color="C5")
plot([], [], label=L"G=d/4", color="C4")
plot([], [], label=L"G=d/2", color="C3")
legend(loc="lower right")

tight_layout()
fig.align_ylabels([a1,a2])
fig.align_ylabels([b1,b2])

# Add skin depth inset
x0 = a2.get_position().x0
y0 = a2.get_position().y0
x1 = a1.get_position().x1
y1 = a1.get_position().y1
w = x1-x0
h = y1-y0
inset = fig.add_axes([x0+w/2.25,y0+h/3,w/2.2,h/3])
lam0s = 2:248 # microns
inset.set_xlabel("Wavelength [μm]", fontsize=7)
inset.set_title("Bulk metal skin depth [nm]", fontsize=7)
inset.set_xticks([2,100,240])
inset.set_xticklabels([2,100,240], fontsize=7)
inset.set_ylim(0,200)
inset.set_yticks([0,200])
inset.set_yticklabels([0,200], fontsize=7)
for (label,metal) in [("Al",LorentzDrudeMetals.Al), ("Ag",LorentzDrudeMetals.Ag), ("Au",LorentzDrudeMetals.Au), ("Ti",LorentzDrudeMetals.Ti)]
    skin_ds = []
    for lam0 in lam0s
        ep = metal[lam0]
        n = conj(sqrt(ep))
        k0 = 2pi/lam0/1000 # 1/nm
        skin_d = 1/2/imag(n)/k0 # nm
        push!(skin_ds, skin_d)
    end
    inset.plot(lam0s, skin_ds, label=label)
end

savefig("figure4.pdf"); println("Figure 4 saved")


######################################
############## Figure 5 ##############
######################################

fig = figure(figsize=(ONECOL,0.5ONECOL))
gs0 = gs.GridSpec(1, 2, wspace=0.9)

a = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(0,0)), hspace=0.1)
a1 = subplot(get(a,(0,0)))
a2 = subplot(get(a,(1,0)))
data = load("./fem/data/figure5a.jld")
aspects = data["aspects"]
for (i,G_by_d) in enumerate([1/2, 1/4, 1/8])
    neffs = data["neffs_G_by_d=$(G_by_d)"]
    neffs_MG = data["neffs_MG_G_by_d=$(G_by_d)"]
    a1.semilogx(aspects, +real(neffs_MG), "--", color="C$(i+2)", alpha=0.5)
    a1.semilogx(aspects, +real(neffs), color="C$(i+2)")
    a2.loglog(aspects, -imag(neffs_MG), "--", color="C$(i+2)", alpha=0.5)
    a2.loglog(aspects, -imag(neffs), color="C$(i+2)")
end
subplot(a1)
title("a", loc="left", ha="right", weight="bold")
ylabel(L"\mathrm{Re}[n_\mathrm{eff}]")
xticks([1/4,1,4], ["","",""])
gca().minorticks_off()
subplot(a2)
xlabel("Aspect ratio "*L"(L_\perp/L_\parallel)")
ylabel(L"\mathrm{Im}[n_\mathrm{eff}]")
xticks([1/4,1,4], [L"1/4",L"1",L"4"])
ylim(1e-6,1e-2)
gca().minorticks_off()
plot([], [], label=L"G=d/2", color="C3")
plot([], [], label=L"G=d/4", color="C4")
plot([], [], label=L"G=d/8", color="C5")
legend(ncol=3, fontsize="small", handlelength=0.5, handletextpad=0.5, columnspacing=1)
bbox = subplot(a1).get_position()
ax1 = fig.add_axes((bbox.x0, bbox.y0, bbox.width, 2*bbox.height))
ax1.set_aspect("equal")
ax1.axis("off")
xlim(-3.3,3.3)
ylim(-1.2,1.2)
scale = 1
rect = patches.Rectangle((-0.5*scale,-0.5*scale), scale, scale, fc="white", ec="black", linewidth=0.75)
ax1.add_artist(rect)
for (ox,oy) in [(-0.25,-0.25),(0.25,-0.25),(-0.25,0.25),(0.25,0.25)]
    circ = patches.Ellipse((ox*scale,oy*scale),0.4*scale,0.4*scale, facecolor="grey")
    ax1.add_artist(circ)
end
rect = patches.Rectangle((2-0.9*scale,-0.3*scale), 1.8*scale, 0.6*scale, fc="white", ec="black", linewidth=0.75)
ax1.add_artist(rect)
for (ox,oy) in [(-0.45,-0.15),(0.45,-0.15),(-0.45,0.15),(0.45,0.15)]
    circ = patches.Ellipse((ox*scale+2,oy*scale),0.4*scale*2,0.4*scale/2, facecolor="grey")
    #circ = patches.Ellipse((ox+2,oy),0.4*2,0.4/2, facecolor="grey")
    ax1.add_artist(circ)
end
rect = patches.Rectangle((-2-0.3*scale,-0.9*scale), 0.6*scale, 1.8*scale, fc="white", ec="black", linewidth=0.75)
ax1.add_artist(rect)
for (ox,oy) in [(-0.15,-0.45),(0.15,-0.45),(-0.15,0.45),(0.15,0.45)]
    circ = patches.Ellipse((ox*scale-2,oy*scale),0.4*scale/2,0.4*scale*2, facecolor="grey")
    ax1.add_artist(circ)
end
ax1.text(0, 1, L"E", ha="center")
ax1.annotate("", xy=(1.5, 0.75), xytext=(-1.5, 0.75), arrowprops=Dict("arrowstyle"=>"<->", "lw"=>0.8))


b = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=get(gs0,(0,1)), hspace=0.1)
b1 = subplot(get(b,(0,0)))
b2 = subplot(get(b,(1,0)))
data = load("./fem/data/figure5b.jld")
ds = data["ds"]
neffs = data["neffs"]
neffs_MG = data["neffs_MG"]
delta_s = data["delta_s"]
subplot(b1)
title("b", loc="left", ha="right", weight="bold")
axvline(delta_s, 0, 0.4, color="black", lw=0.75, ls="--")
text(10, 2, L"\delta_s", va="bottom")
plot(ds, +real(neffs_MG), "--", color="C0", alpha=0.7)
plot(ds, +real(neffs), color="C0")
ylabel(L"\mathrm{Re}(n_\mathrm{eff})")
xticks([0,25,48], ["","",""]);
subplot(b2)
axvline(delta_s, color="black", lw=0.75, ls="--")
semilogy(ds, -imag(neffs_MG), "--", color="C0", alpha=0.7, label="Maxwell Garnett")
semilogy(ds, -imag(neffs), color="C0", label="Numerics")
xlabel("Diameter [nm]")
ylabel(L"\mathrm{Im}(n_\mathrm{eff})")
subplots_adjust(hspace=0.1)
xticks([0,25,48], ["0","25","48"]);
ylim(1e-6,1e-1)
bbox = subplot(b1).get_position()
ax1 = fig.add_axes((bbox.x0, bbox.y0, bbox.width, 2*bbox.height))
ax1.set_aspect("equal")
ax1.axis.("off")
xlim(0,48)
ylim(-10,10)
for d in [10, 25, 40]
    hex = patches.RegularPolygon((d,0), 6, 1.6*3, fc="white", ec="black", lw=0.75, orientation=pi/2)
    ax1.add_artist(hex)
    for (ix,iy) in [(0,0),(1,0),(-1,0),(-0.5,0.5),(0.5,0.5),(-0.5,-0.5),(0.5,-0.5)]
        ox = ix*3
        oy = iy*3*sqrt(3)
        circ = patches.Circle((d+ox,oy), 1.5*d/50, fc="grey")
        ax1.add_artist(circ)
    end
end

fig.align_ylabels([a1,a2])
fig.align_ylabels([b1,b2])
fig.align_xlabels([a2,b2])
savefig("figure5.pdf", bbox_inches="tight"); println("Figure 5 saved")


#show()
