using ProgressMeter
using SparseArrays
using LinearAlgebra

function getyeemats(Nx,Ny,Nz,k0dx,k0dy,k0dz)
    N = Nx*Ny*Nz
    Dx = ones(Nx,Ny,Nz)
    Dy = ones(Nx,Ny,Nz)
    Dz = ones(Nx,Ny,Nz)
    p1 = CartesianIndex(1,1,1)
    p2 = CartesianIndex(Nx,Ny,Nz)
    for p in p1:p2
        if p[1] == Nx
            Dx[p] = 0
        end
        if p[2] == Ny
            Dy[p] = 0
        end
        if p[3] == Nz
            Dz[p] = 0
        end
    end
    DxE = spdiagm(0=>-ones(N), 1=>Dx[1:N-1], 1-Nx=>(1.0.-Dx)[Nx:N]) / k0dx
    DyE = spdiagm(0=>-ones(N), Nx=>Dy[1:N-Nx], Nx-Nx*Ny=>(1.0.-Dy)[1+Nx*Ny-Nx:N]) / k0dy
    DzE = spdiagm(0=>-ones(N), Nx*Ny=>Dz[1:N-Nx*Ny]) / k0dz # no periodic BC for z direction
    #DzE = spdiagm([-ones(N),dz[1:N-Nx*Ny],(1-dz)[1+Nx*Ny*Nz-Nx*Ny:N]],[0,Nx*Ny,Nx*Ny-Nx*Ny*Nz]) / k0dz
    DxH, DyH, DzH = -DxE', -DyE', -DzE'
    return DxE, DyE, DzE, DxH, DyH, DzH, spzeros(N,N)
end

function reorder(A)
    N = size(A, 1)
    order = map(i->1+div(N,3)*mod(i-1,3)+div(i-1,3), 1:N)
    if ndims(A) == 1
        return A[order]
    elseif ndims(A) == 2
        return A[order,order]
    end
end

function getA(ep,mu,yeemats)
    Nx = size(ep, 1)
    Ny = size(ep, 2)
    Nz = size(ep, 3)
    DxE, DyE, DzE, DxH, DyH, DzH, Z = yeemats
    if Ny == 1
        epxx = ep[:,:,:,1]
        muyy = mu[:,:,:,2]
        epzz = ep[:,:,:,3]
        return DzE*spdiagm(0=>1.0./epxx[:])*DzH + DxE*spdiagm(0=>1.0./epzz[:])*DxH + spdiagm(0=>muyy[:])
    else
        CE = [Z -DzE DyE; DzE Z -DxE; -DyE DxE Z]
        CH = [Z -DzH DyH; DzH Z -DxH; -DyH DxH Z]
        #A = CH*spdiagm(1.0/mu[:])*CE - spdiagm(ep[:])
        A = CE*spdiagm(0=>1.0./ep[:])*CH - spdiagm(0=>mu[:])
        return reorder(A)
    end
end

function getsliceabc(ep, mu, iz, yeemats)
    crop(arr,iz) = circshift(arr,(0,0,2-iz,0))[:,:,1:3,:]
    Nx = size(ep, 1)
    Ny = size(ep, 2)
    Nz = size(ep, 3)
    if Ny == 1
        G = Nx*Ny
    else
        G = 3*Nx*Ny
    end
    A = getA(crop(ep,iz),crop(mu,iz),yeemats)
    ###### would be better to directly fill a, b, c
    ###### from A without reordering a sparse matrix
    a = A[(1:G).+1G, (1:G).+0G]
    b = A[(1:G).+1G, (1:G).+1G]
    c = A[(1:G).+1G, (1:G).+2G]
    return Array(a), Array(b), Array(c)
end

function getblochabc(ep, mu, yeemats)
    Nx = size(ep, 1)
    Ny = size(ep, 2)
    Nz = size(ep, 3)
    b0 = getsliceabc(ep, mu, 1, yeemats)[2]
    a1, b1, c1 = getsliceabc(ep, mu, 1, yeemats)
    a2, b2, c2 = getsliceabc(ep, mu, 2, yeemats)
    for iz in 3:Nz+1
        # Calculate next slice
        a3, b3, c3 = getsliceabc(ep, mu, iz, yeemats)
        # Absorb middle slice
        b2lu = factorize(b2)
        b2lua2 = b2lu\a2
        b2luc2 = b2lu\c2
        b1 = b1 - c1*b2lua2
        c1 =    - c1*b2luc2
        b3 = b3 - a3*b2luc2
        a3 =    - a3*b2lua2
        # Discard and relabel
        a2, b2, c2 = a3, b3, c3
    end
    a = a2
    c = c1
    b = b0 + (b1-b0) + (b2-b0)
    return a, b, c
end

function getneff(ep, mu, dx, dy, dz, k0; withabc=false)
    Nx = size(ep, 1)
    Ny = size(ep, 2)
    Nz = size(ep, 3)
    yeemats = getyeemats(Nx, Ny, 3, k0*dx, k0*dy, k0*dz)
    a, b, c = getblochabc(ep, mu, yeemats)
    # First estimate
    myeigs, myvecs = eigen(-b\(a+c))
    neffs = acos.(1.0./myeigs.+0im) / (k0*Nz*dz)
    sort!(neffs, by=x->abs(imag(x)))
    # Make real part positive
    neffs = map(x->real(x)>0 ? x : -x, neffs)
    # Ensure imag part is negative
    neffs = map(x->imag(x)<1e-10 ? x : -x, neffs)
    # Least decaying modes first
    sort!(neffs, by=x->abs(imag(x)))
    return withabc ? (neffs,a,b,c) : neffs
end

function iterateneff(neff, a, b, c, k0, Lz)
    neffs = [neff]
    K = neff*k0
    for i in 1:100
        M = a*exp(1im*K*Lz) + b + c*exp(-1im*K*Lz)
        Λ = sort(eigvals(M), by=abs)[1]
        #@show K
        #@show Λ
        tmp = eigvals(a*exp(1im*K*Lz)-c*exp(-1im*K*Lz))
        Δs = -Λ ./ tmp / 1im / Lz
        #@show Δs[1]
        #@show sort(eigvals(M(K+Δs[1])), by=abs)[1]
        K = K + sort(Δs, by=abs)[1]
        push!(neffs, K/k0) ##### would be better to do the calc directly with neff?
    end
    return neffs
end


function makecyls(rs, eps, dr)
    Lx, Lz = 2*sum(rs), 2*sum(rs)
    dx, dz = dr, dr

    Nx, Nz = ceil(Int,Lx/dx), ceil(Int,Lz/dz)
    dx, dz = Lx/Nx, Lz/Nz
    
    Ny = 1
    dy = dx
    Ly = Ny*dy

    ep = ones(ComplexF64, Nx, Ny, Nz, 3)
    mu = ones(ComplexF64, Nx, Ny, Nz, 3)
    ep2 = ones(ComplexF64, 2Nx, 2Ny, 2Nz)
    x2axis = range(0, stop=(Nx-1)*dx, length=2Nx)
    y2axis = range(0, stop=(Ny-1)*dy, length=2Ny)
    z2axis = range(0, stop=(Nz-1)*dz, length=2Nz)
    x2axis = x2axis
    y2axis = y2axis
    z2axis = z2axis
    origins = [
        (Lx/2,0,Lz/2),
    ]
    ep2 .= eps[end]
    for (x0,y0,z0) in origins
        for (iz,z) in enumerate(z2axis), (iy,y) in enumerate(y2axis), (ix,x) in enumerate(x2axis)
            for (r,epr) in zip(cumsum(rs),eps)
                if (x-x0)^2+(z-z0)^2 < r^2
                    ep2[ix,iy,iz] = epr
                    break
                end
            end
        end
    end
    ep[:,:,:,1] = ep2[2:2:end,1:2:end,1:2:end]
    ep[:,:,:,2] = ep2[1:2:end,2:2:end,1:2:end]
    ep[:,:,:,3] = ep2[1:2:end,1:2:end,2:2:end];
    return ep, mu, dx, dy, dz
end

function makecylshex(rs, eps, dr)
    Lx, Lz = 2*sum(rs), 2*sum(rs)*sqrt(3)
    dx, dz = dr, dr

    Nx, Nz = ceil(Int,Lx/dx), ceil(Int,Lz/dz)
    dx, dz = Lx/Nx, Lz/Nz
    
    Ny = 1
    dy = dx
    Ly = Ny*dy

    ep = ones(ComplexF64, Nx, Ny, Nz, 3)
    mu = ones(ComplexF64, Nx, Ny, Nz, 3)
    ep2 = ones(ComplexF64, 2Nx, 2Ny, 2Nz)
    x2axis = range(0, stop=(Nx-1)*dx, length=2Nx)
    y2axis = range(0, stop=(Ny-1)*dy, length=2Ny)
    z2axis = range(0, stop=(Nz-1)*dz, length=2Nz)
    x2axis = x2axis
    y2axis = y2axis
    z2axis = z2axis
    origins = [
        (Lx/2,0,Lz/2),
        (0,0,0),
        (Lx,0,0),
        (0,0,Lz),
        (Lx,0,Lz),
    ]
    ep2 .= eps[end]
    for (x0,y0,z0) in origins
        for (iz,z) in enumerate(z2axis), (iy,y) in enumerate(y2axis), (ix,x) in enumerate(x2axis)
            for (r,epr) in zip(cumsum(rs),eps)
                if (x-x0)^2+(z-z0)^2 < r^2
                    ep2[ix,iy,iz] = epr
                    break
                end
            end
        end
    end
    ep[:,:,:,1] = ep2[2:2:end,1:2:end,1:2:end]
    ep[:,:,:,2] = ep2[1:2:end,2:2:end,1:2:end]
    ep[:,:,:,3] = ep2[1:2:end,1:2:end,2:2:end];
    return ep, mu, dx, dy, dz
end

function makespheres(rs, eps, dr)
    throw("Not implemented")
end

function makesphereshex(rs, eps, dr)
    Lx, Ly, Lz = 2*sum(rs), 2*sum(rs)*sqrt(3), 2*sum(rs)*sqrt(3)
    dx, dy, dz = dr, dr, dr

    Nx, Ny, Nz = ceil(Int,Lx/dx), ceil(Int,Ly/dy), ceil(Int,Lz/dz)
    dx, dy, dz = Lx/Nx, Ly/Ny, Lz/Nz

    ep = ones(ComplexF64, Nx, Ny, Nz, 3)
    mu = ones(ComplexF64, Nx, Ny, Nz, 3)
    ep2 = ones(ComplexF64, 2Nx, 2Ny, 2Nz)
    x2axis = range(0, stop=(Nx-1)*dx, length=2Nx)
    y2axis = range(0, stop=(Ny-1)*dy, length=2Ny)
    z2axis = range(0, stop=(Nz-1)*dz, length=2Nz)
    x2axis = x2axis
    y2axis = y2axis
    z2axis = z2axis
    origins = [
    (0,0,0),
    (Lx,0,0),
    (0,Ly,0),
    (Lx,Ly,0),
    (Lx/2,Ly/2,0),
    (Lx/2,Ly/6,Lz/2),
    (Lx/2,7Ly/6,Lz/2),
    (0,2Ly/3,Lz/2),
    (Lx,2Ly/3,Lz/2),
    (0,0,Lz),
    (Lx,0,Lz),
    (0,Ly,Lz),
    (Lx,Ly,Lz),
    (Lx/2,Ly/2,Lz)
    ]
    ep2 .= eps[end]
    for (x0,y0,z0) in origins
        for (iz,z) in enumerate(z2axis), (iy,y) in enumerate(y2axis), (ix,x) in enumerate(x2axis)
            for (r,epr) in zip(cumsum(rs),eps)
                if (x-x0)^2+(y-y0)^2+(z-z0)^2 < r^2
                    ep2[ix,iy,iz] = epr
                    break
                end
            end
        end
    end
    ep[:,:,:,1] = ep2[2:2:end,1:2:end,1:2:end]
    ep[:,:,:,2] = ep2[1:2:end,2:2:end,1:2:end]
    ep[:,:,:,3] = ep2[1:2:end,1:2:end,2:2:end];
    return ep, mu, dx, dy, dz
end


function makecubes(rs, eps, dr)
    Lx, Ly, Lz = 2*sum(rs), 2*sum(rs), 2*sum(rs)
    dx, dy, dz = dr, dr, dr

    Nx, Ny, Nz = ceil(Int,Lx/dx), ceil(Int,Ly/dy), ceil(Int,Lz/dz)
    dx, dy, dz = Lx/Nx, Ly/Ny, Lz/Nz

    ep = ones(ComplexF64, Nx, Ny, Nz, 3)
    mu = ones(ComplexF64, Nx, Ny, Nz, 3)
    ep2 = ones(ComplexF64, 2Nx, 2Ny, 2Nz)
    x2axis = range(0, stop=(Nx-1)*dx, length=2Nx)
    y2axis = range(0, stop=(Ny-1)*dy, length=2Ny)
    z2axis = range(0, stop=(Nz-1)*dz, length=2Nz)
    x2axis = x2axis
    y2axis = y2axis
    z2axis = z2axis
    ep2 .= eps[end]
    for (x0,y0,z0) in [(Lx/2,Ly/2,Lz/2)]
        for (iz,z) in enumerate(z2axis), (iy,y) in enumerate(y2axis), (ix,x) in enumerate(x2axis)
            for (r,epr) in zip(cumsum(rs),eps)
                if abs(x-x0)<r && abs(y-y0)<r && abs(z-z0)<r
                    ep2[ix,iy,iz] = epr
                    break
                end
            end
        end
    end
    ep[:,:,:,1] = ep2[2:2:end,1:2:end,1:2:end]
    ep[:,:,:,2] = ep2[1:2:end,2:2:end,1:2:end]
    ep[:,:,:,3] = ep2[1:2:end,1:2:end,2:2:end];
    return ep, mu, dx, dy, dz
end
