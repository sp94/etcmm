using LightGraphs
using SparseArrays
using SuiteSparse

mutable struct Node
    x::Real
    y::Real
    xy::Array{Real,1}
    region::Symbol # region (:boundary, :below, or :above)
    mi::Int # mesh index
    ri::Int # region index
    Node(x::Real,y::Real) = new(x,y,[x,y])
end


mutable struct Element
    no1::Node
    no2::Node
    no3::Node
    nos::Array{Node}
    ep::ComplexF64
    mu::ComplexF64
    function Element(no1::Node,no2::Node,no3::Node,ep::Number,mu::Number)
        s1 = no3.xy - no2.xy
        s2 = no1.xy - no3.xy
        s3 = no2.xy - no1.xy
        detJ = s2[1]*s3[2] - s2[2]*s3[1] # = 2 * element area
        if detJ > 0
            return new(no1,no2,no3,[no1,no2,no3],ep,mu)
        else
            return new(no3,no2,no1,[no3,no2,no1],ep,mu)
        end
    end
end

function getAe(el::Element, k0::Real, polarization::Symbol)
    # Construct the local element 'stiffness' matrix
    if polarization == :TE
        alpha = 1/el.ep
        beta = -el.mu*k0^2
    elseif polarization == :TM
        alpha = 1/el.mu
        beta = -el.ep*k0^2
    else
        error("Unknown polarization, $(polarization). Must be :TE or :TM.")
    end
    
    s1 = el.no3.xy - el.no2.xy
    s2 = el.no1.xy - el.no3.xy
    s3 = el.no2.xy - el.no1.xy
    detJ = s2[1]*s3[2] - s2[2]*s3[1] # = 2 * element area
    @assert detJ > 0
    
    # Precomputed integrals of the local basis functions
    # int_phi2[i,j] = int(phi_i*phi_j)
    int_phi2 = [2 1 1; 1 2 1; 1 1 2] / 24
    
    # Gradients of the local basis functions
    grad_phi1e = [-s1[2]; s1[1]] / detJ
    grad_phi2e = [-s2[2]; s2[1]] / detJ
    grad_phi3e = [-s3[2]; s3[1]] / detJ
    grad_phi = [grad_phi1e grad_phi2e grad_phi3e]

    # Construct element matrix
    Ae = Array{ComplexF64}(undef, 3, 3)
    for row in 1:3, col in 1:3
        Ae[row,col] = alpha*detJ*dot(grad_phi[:,row],grad_phi[:,col])/2
    end
    Ae += beta*detJ*int_phi2
    return Ae
end


mutable struct Mesh
    nos::Array{Node}
    els::Array{Element}
    xlim::Array{Real}
    ylim::Array{Real}
    function Mesh(nos::Array{Node}, els::Array{Element})
        # Get mesh dimensions
        xlim = [+Inf, -Inf]
        ylim = [+Inf, -Inf]
        for no in nos
            xlim[1] = min(no.x, xlim[1])
            xlim[2] = max(no.x, xlim[2])
            ylim[1] = min(no.y, ylim[1])
            ylim[2] = max(no.y, ylim[2])
        end
        # Check periodicity
        lnos, rnos = Node[], Node[]
        bnos, tnos = Node[], Node[]
        for no in nos
            if no.x == xlim[1]
                push!(lnos, no)
            elseif no.x == xlim[2]
                push!(rnos, no)
            end
            if no.y == ylim[1]
                push!(bnos, no)
            elseif no.y == ylim[2]
                push!(tnos, no)
            end
        end
        sort!(lnos, by=node->node.y)
        sort!(rnos, by=node->node.y)
        sort!(bnos, by=node->node.x)
        sort!(tnos, by=node->node.x)
        @assert isapprox([no.y for no in lnos], [no.y for no in rnos])
        @assert isapprox([no.x for no in bnos], [no.x for no in tnos])
        rnos2lnos = Dict(zip(rnos,lnos))
        tnos2bnos = Dict(zip(tnos,bnos))

        ri_boundary = 1
        ri_notboundary = 1
        for node in nos
            # Periodic wrapping
            wrapped_node = get(rnos2lnos, node, node)
            wrapped_node = get(tnos2bnos, wrapped_node, wrapped_node)
            if wrapped_node.ri == 0
                if wrapped_node.region == :boundary
                    wrapped_node.ri = ri_boundary
                    ri_boundary += 1
                else
                    wrapped_node.ri = ri_notboundary
                    ri_notboundary += 1
                end
            end
            node.ri = wrapped_node.ri
        end

        new(nos, els, xlim, ylim)
    end
end

function Mesh(meshname::String, eps::Array{<:Number,1}, mus::Array{<:Number,1}; mirror::Bool=false)
    
    msh = meshio.read("$(meshname)")
    p = msh.points
    cells = msh.cells
    cell_data = msh.cell_data

    t = cells["triangle"] .+ 1 # compensate for zero based indexing in meshio
    m = cell_data["triangle"]["gmsh:physical"]
    boundary = unique(cells["line"]) .+ 1

    # Store nodes
    nos = Node[]
    for ino in 1:size(p,1)
        x, y = p[ino,:]
        node = Node(x,y)
        node.mi = ino
        node.region = :undefined
        node.ri = 0
        push!(nos, node)
    end
    # Store elements
    els = Element[]
    for iel in 1:size(t,1)
        no1, no2, no3 = nos[t[iel,:]]
        ep = eps[m[iel]]
        mu = mus[m[iel]]
        push!(els, Element(no1, no2, no3, ep, mu))
    end
    # Apply mirror operator to mesh
    if mirror
        # Mirror in x = 0 and y = 0
        mirrors = [Node[] for no in nos]
        next_mi = length(nos)+1
        for node in nos[:]
            if node.x > 0
                new_node = Node(-node.x,node.y)
                new_node.mi = next_mi
                new_node.region = :undefined
                new_node.ri = 0
                next_mi += 1
                push!(nos, new_node)
                push!(mirrors[node.mi], new_node)
                if node.mi in boundary
                    push!(boundary, new_node.mi)
                end
            else
                push!(mirrors[node.mi], node)
            end
            if node.y > 0
                new_node = Node(node.x,-node.y)
                new_node.mi = next_mi
                new_node.region = :undefined
                new_node.ri = 0
                next_mi += 1
                push!(nos, new_node)
                push!(mirrors[node.mi], new_node)
            else
                push!(mirrors[node.mi], node)
            end
            if node.x > 0 && node.y > 0
                new_node = Node(-node.x,-node.y)
                new_node.mi = next_mi
                new_node.region = :undefined
                new_node.ri = 0
                next_mi += 1
                push!(nos, new_node)
                push!(mirrors[node.mi], new_node)
            elseif node.x > 0
                push!(mirrors[node.mi], mirrors[node.mi][1])
            elseif node.y > 0
                push!(mirrors[node.mi], mirrors[node.mi][2])
            else
                push!(mirrors[node.mi], node)
            end
        end
        for el in els[:]
            node1, node2, node3 = el.nos
            for (new_node1,new_node2,new_node3) in zip(mirrors[node1.mi],mirrors[node2.mi],mirrors[node3.mi])
                # Skip the case where no nodes have been mirrored
                if [new_node1,new_node2,new_node3] != [node1,node2,node3]
                    # If at least one node has been mirrored, create a new element
                    new_element = Element(new_node1, new_node2, new_node3, el.ep, el.mu)
                    push!(els, new_element)
                end
            end
        end
    end
    # Detect regions
    for mi in boundary
        nos[mi].region = :boundary
    end
    # Construct graph
    g = Graph(length(nos))
    for iel in 1:length(els)
        el = els[iel]
        for (no1,no2) in zip(el.nos,circshift(el.nos,1))
            add_edge!(g, no1.mi, no2.mi)
        end
    end
    # Separate the regions in our graph
    for edge in collect(edges(g))
        node1 = nos[src(edge)]
        node2 = nos[dst(edge)]
        # Delete edges that connect boundary nodes to non-boundary nodes
        if xor(node1.region == :boundary, node2.region == :boundary)
            rem_edge!(g, edge)
        end
    end
    regions = connected_components(g)
    #@show length(regions)
    for region in regions
        if nos[region[1]].region == :boundary
            continue
        end
        ys_boundary = [nos[mi].y for mi in boundary]
        ys_region = [nos[mi].y for mi in region]
        if minimum(ys_region) < minimum(ys_boundary)
            region_symbol = :below
        else
            region_symbol = :above
        end
        for mi in region
            nos[mi].region = region_symbol
        end
    end
    #@assert length(regions) == 3
    Mesh(nos, els)
end

Base.show(io::IO, m::Mesh) = print(io, "Mesh with $(length(m.nos)) nodes, $(length(m.els)) elements")

function mplot(m::Mesh; tile::Int=1)
    figure()
    subplot(1,2,1)
    gca().set_aspect("equal", "datalim")
    subplot(1,2,2)
    gca().set_aspect("equal", "datalim")
    for ix in 1:tile, iy in 1:tile
        Lx = m.xlim[2] - m.xlim[1]
        Ly = m.ylim[2] - m.ylim[1]
        xs = [no.x for no in m.nos] .+ (ix-1)*Lx
        ys = [no.y for no in m.nos] .+ (iy-1)*Ly
        eps = [el.ep for el in m.els]
        mus = [el.mu for el in m.els]
        tris = [el.nos[ino].mi for el in m.els, ino in 1:3]
        subplot(1,2,1)
        tripcolor(xs, ys, tris-1, real(eps), edgecolors="k", alpha=0.5)
        subplot(1,2,2)
        tripcolor(xs, ys, tris-1, real(mus), edgecolors="k", alpha=0.5)
        xs = []
        ys = []
        for no in m.nos
            if no.region == :boundary
                push!(xs, no.x)
                push!(ys, no.y)
            end
        end
        plot(xs, ys, "r.")
    end
    title("$(length(m.nos)) nodes")
end


function uniformtris(Nx, Ny, dx, dy; distortion=0)
    
    next_mi = 1
    nodes = Dict{Tuple{Int,Int},Node}()
    nos = Node[]
    for ix in 1:Nx+1, iy in 1:Ny+1
        node = Node((ix-1)*dx,(iy-1)*dy)
        if 1 < ix < Nx+1 && 1 < iy < Ny+1
            node.x += distortion*rand()
            node.y += distortion*rand()
        end
        node.mi = next_mi
        node.ri = 0
        next_mi += 1
        nodes[ix,iy] = node
        push!(nos, node)
    end

    for ix in 1:Nx+1, iy in 1:Ny+1
        node = nodes[ix,iy]
        if iy == 1+div(Ny,2)
            node.region = :boundary
        end
        if iy < 1+div(Ny,2)
            node.region = :below
        end
        if iy > 1+div(Ny,2)
            node.region = :above
        end
    end

    els = Element[]
    for ix in 1:Nx, iy in 1:Ny
        no1 = nodes[ix,iy]
        no2 = nodes[ix+1,iy]
        no3 = nodes[ix,iy+1]
        no1, no2, no3 = shuffle([no1, no2, no3])
        push!(els, Element(no1,no2,no3,1,1))
        no1 = nodes[ix+1,iy]
        no2 = nodes[ix+1,iy+1]
        no3 = nodes[ix,iy+1]
        no1, no2, no3 = shuffle([no1, no2, no3])
        push!(els, Element(no1,no2,no3,1,1))
    end

    return Mesh(nos, els)
end


function getabc(m::Mesh, k0::Real, polarization::Symbol)
    
    boundary_nodes = filter(node->node.region==:boundary, m.nos)
    notboundary_nodes = filter(node->node.region!=:boundary, m.nos)
    nk = length(unique([node.ri for node in boundary_nodes]))
    na = length(unique([node.ri for node in notboundary_nodes]))
    
    # It would be better to initialise these from sparse(row,col,val)
    b1 = spzeros(ComplexF64, na, na)
    a2 = spzeros(ComplexF64, nk, na)
    b2 = spzeros(ComplexF64, nk, nk)
    c2 = spzeros(ComplexF64, nk, na)
    for el in m.els
        Ae = getAe(el, k0, polarization)
        for row in 1:3, col in 1:3
            node1 = el.nos[row]
            node2 = el.nos[col]
            if node1.region == :boundary
                if node2.region == :boundary
                    b2[node1.ri,node2.ri] += Ae[row,col]
                elseif node2.region == :above
                    a2[node1.ri,node2.ri] += Ae[row,col]
                elseif node2.region == :below
                    c2[node1.ri,node2.ri] += Ae[row,col]
                end
            elseif node1.region != :boundary && node2.region != :boundary
                b1[node1.ri,node2.ri] += Ae[row,col]
            end
        end
    end
    a1 = transpose(c2)
    c1 = transpose(a2)
    return a1, b1, c1, a2, b2, c2
end

function getblochabc(a1, b1, c1, a2, b2, c2)
    b1fact = lu(b1)
    function mysolve(A,B)
        out = zeros(ComplexF64, size(A,1), size(B,2))
        for col in 1:size(B,2)
            out[:,col] = A \ Array(B[:,col])
        end
        return out
    end
    b1_inv_a1 = mysolve(b1fact, a1)
    b1_inv_c1 = mysolve(b1fact, c1)
    a = -a2*b1_inv_a1
    b = b2 - a2*b1_inv_c1 - c2*b1_inv_a1
    c = -c2*b1_inv_c1
    return a, b, c
end

function homogenize(m::Mesh, k0::Real, polarization::Symbol)
    a1, b1, c1, a2, b2, c2 = getabc(m, k0, polarization)
    a, b, c = getblochabc(a1, b1, c1, a2, b2, c2)
    #@show isapprox(a,c'), isapprox(b,b')
    #@show isapprox(a,transpose(c)), isapprox(b,transpose(b))
    vals = eigvals(-b\(a+c))
    neffs = acos.(1.0./vals) / k0 / (m.ylim[2]-m.ylim[1])
    # Sort by speed of decay (or growth)
    sort!(neffs, by=x->abs(imag(x)))
    # Make real part positive
    map!(neff->real(neff)<0 ? -neff : neff, neffs, neffs)
    # But swap sign if imag(neff) > 0 ~ 1e-15 with tolerance for numerical error
    map!(neff->imag(neff)>1e-15 ? -neff : neff, neffs, neffs)
    return neffs, a, b, c
end

function visualize(m::Mesh, k0::Real, polarization::Symbol)

    @assert polarization == :TE # (calculating E_x, E_y from H_z)
    # Might need a sign change for calculating H_x, H_y from E_z, also switch ep to mu

    # Get Fz
    a1, b1, c1, a2, b2, c2 = getabc(m, k0, polarization)
    a, b, c = getblochabc(a1, b1, c1, a2, b2, c2)
    vals, vecs = eig(-b\(a+c))
    KLs = acos.(1.0./vals)
    # Sort by speed of decay (or growth)
    idx = sortperm(KLs, by=x->abs(imag(x)))
    KL = KLs[idx][1]
    F2 = vecs[:,idx][:,1]
    F2 = F2 / maximum(abs.(F2))
    # Solve boundary value problem
    RHS = -a1*F2*exp(1im*KL)-c1*F2
    F1 = b1\RHS
    Fz = zeros(ComplexF64, length(m.nos))
    for no in m.nos
        if no.region == :boundary
            Fz[no.mi] = F2[no.ri]
        else
            if no.region == :above
                Fz[no.mi] = F1[no.ri]
            else
                Fz[no.mi] = F1[no.ri]*exp(-1im*KL)
            end
        end
    end

    # Get Fx, Fy
    Fx = zeros(ComplexF64, length(m.els))
    Fy = zeros(ComplexF64, length(m.els))
    for (iel,el) in enumerate(m.els)
        s1 = el.no3.xy - el.no2.xy
        s2 = el.no1.xy - el.no3.xy
        s3 = el.no2.xy - el.no1.xy
        detJ = s2[1]*s3[2] - s2[2]*s3[1] # = 2 * element area
        @assert detJ > 0
        # here: switched omega -> k0, assuming the normalisation ep0 = mu0 = 1
        Fx[iel] = (Fz[el.no1.mi]*s1[1] + Fz[el.no2.mi]*s2[1] + Fz[el.no3.mi]*s3[1]) / (detJ * 1im*k0*el.ep)
        Fy[iel] = (Fz[el.no1.mi]*s1[2] + Fz[el.no2.mi]*s2[2] + Fz[el.no3.mi]*s3[2]) / (detJ * 1im*k0*el.ep)
    end

    xs = [no.x for no in m.nos]
    ys = [no.y for no in m.nos]
    tris = [el.nos[ino].mi for el in m.els, ino in 1:3]
    return xs, ys, tris, Fx, Fy, Fz
end
