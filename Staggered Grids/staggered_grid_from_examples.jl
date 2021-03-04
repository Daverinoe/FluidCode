using Plots
# using GLMakie

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# Domain size and physical variables
Lx = 1
Ly = 1
gx = 0.0
gy = -600.0
r1 = 1.0
r2 = 2
rro = r1

# Dynamic viscosity
m0 = 0.005

# Tangential velocities
uNorth = 0.0
uSouth = 2.0
vWest = 0.0
vEast = 0.0

# Numerical variables
nx = 100
ny = 100
xc = floor(Int, nx/2)
yc = floor(Int, ny/2)
# rad = floor(Int, sqrt(nx^2 + ny^2)/10)
rad = 12
dt = 0.1e-3
nstep = 5000
maxiter = 200
maxError = 2e-3
beta = 1.2

u = zeros(nx+1, ny+2)
v = zeros(nx+2, ny+1)
p = zeros(nx+2, ny+2)

ut = zeros(nx+1, ny+2)
vt = zeros(nx+2, ny+1)

tmp1 = zeros(nx+2, ny+2)
tmp2 = zeros(nx+2, ny+2)

# Velocities at center of grid for plotting
uu = zeros(nx+1, ny+1)
vv = zeros(nx+1, ny+1)

# Define the grid
dx = Lx / nx
dy = Ly / ny

x = zeros(nx+2)
y = zeros(nx+2)

for i=1:nx+2
    x[i] = dx*(i-1.5)
end
for j=1:ny+2
    y[j] = dy*(j-1.5)
end

xh = zeros(nx+1)
yh = zeros(ny+1)
for i = 1:nx+1
    xh[i] = dx * (i-1)
end
for j = 1:ny+1
    yh[j] = dy * (j-1)
end

r = zeros(nx+2, ny+2) .+ r1
# i_base = floor(Int, nx/2)
# for j = 1:ny+2
#     for i = 1:nx+2
#         if i > i_base + 0.1*ny*cos(1*2*pi*j/ny)
#             r[i, j] = r2 # Set density of drop
#         end
#     end
# end

for j = 1:ny+2
    for i = 1:nx+2
        if (i - xc)^2 + (j - yc)^2 <= rad^2
            r[i, j] = r2 # Set density of drop
        end
    end
end

time = 0.0

@time anim = @animate for steps in range(1, stop = nstep)

    # Set the tangential velocities at the boundaries.
    u[:, 1] = 2 * uSouth .- u[:, 2]
    u[:, end] = 2 * uNorth .- u[:, end-1]
    v[1, :] = 2 * vWest .- v[1, :]
    v[end, :] = 2 * vEast .- v[end-1, :]

    # u[:, 1] = u[:, 2]
    # u[:, end] = u[:, end-1]
    # v[1, :] = v[2, :]
    # v[end, :] = v[end-1, :]

    # u[:, 1] = u[:, end-1]
    # u[:, end] = u[:, 2]
    # v[1, :] = v[end-1, :]
    # v[end, :] = v[2, :]

    # r[:, 2] = r[:, end-2]
    # r[:, end-1] = r[:, 3]

    # Temporary velocities
    # U-velocity
    for i in range(2, stop = nx)
        for j in range(2, stop = ny+1)
            reynolds = m0/(0.5*(r[i+1,j]+r[i,j]))
            # reynolds = u[i,j]*Lx/m0
            ut[i,j]=u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])^2-(u[i,j]+u[i-1,j])^2)/dx+((u[i,j+1]+u[i,j])*(v[i+1,j]+v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/dy)+reynolds*((u[i+1,j]-2*u[i,j]+u[i-1,j])/dx^2+(u[i,j+1]-2*u[i,j]+u[i,j-1])/dy^2 )+gx)
        end
    end

    for i in range(2, stop = nx+1)
        for j in range(2, stop = ny)
            reynolds = m0/(0.5*(r[i+1,j]+r[i,j]))
            # reynolds = v[i,j]*Ly/m0
            vt[i,j]=v[i,j]+dt*(-0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j]+v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/dx+((v[i,j+1]+v[i,j])^2-(v[i,j]+v[i,j-1])^2)/dy)+reynolds*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx^2+(v[i,j+1]-2*v[i,j]+v[i,j-1])/dy^2 )+gy)
        end
    end

    # Compute source term and coefficient for p[i, j]
    rt = copy(r)
    lrg = 1000.0
    rt[:, 1] .= lrg
    rt[:, end] .= lrg
    rt[1, :] .= lrg
    rt[end, :] .= lrg

    for i in range(2, stop = nx+1)
        for j in range(2, stop = ny+1)
            tmp1[i,j]= (0.5/dt)*( (ut[i,j]-ut[i-1,j])/dx+(vt[i,j]-vt[i,j-1])/dy )
            tmp2[i,j]=1.0/( (1 ./ dx)*( 1 ./ (dx*(rt[i+1,j]+rt[i,j]))+1 ./ (dx*(rt[i-1,j]+rt[i,j])))+(1 ./ dy)*(1 ./ (dy*(rt[i,j+1]+rt[i,j]))+1 ./ (dy*(rt[i,j-1]+rt[i,j]))))
        end
    end

    iter = 0
    while true
        pn = copy(p)
        iter = iter+1
        for i in range(2, stop = nx+1)
            for j in range(2, stop = ny+1)
                p[i,j]=(1.0-beta)*p[i,j]+beta* tmp2[i,j]*((1 ./ dx)*( p[i+1,j]/(dx*(rt[i+1,j]+rt[i,j]))+p[i-1,j]/(dx*(rt[i-1,j]+rt[i,j])))+(1 ./ dy)*( p[i,j+1]/(dy*(rt[i,j+1]+rt[i,j]))+p[i,j-1]/(dy*(rt[i,j-1]+rt[i,j])))-tmp1[i,j])
            end
        end
        if maximum(abs.(pn.-p))<maxError
            break
        end
        if iter>maxiter
            break
        end
    end

    #CORRECT THE u-velocity
    for i in range(2, stop = nx)
        for j in range(2, stop = ny+1)
            u[i,j]=ut[i,j]-dt*(2.0/dx)*(p[i+1,j]-p[i,j])/(r[i+1,j]+r[i,j])
        end
    end

    #CORRECT THE v-velocity
    for i in range(2, stop = nx+1)
        for j in range(2, stop = ny)
            v[i,j] = vt[i,j]-dt*(2.0/dy)*(p[i,j+1]-p[i,j])/(r[i,j+1]+r[i,j])
        end
    end

# ADVECT DENSITY using centered difference plus diffusion
    ro = copy(r)
    for i in range(2, stop = nx+1)
        for j in range(2, stop = ny+1)
            reynolds = m0
            # reynolds = (u[i,j])*(v[i,j])*r[i,j]*Lx*Ly/m0
            r[i,j]=ro[i,j]-(0.5*dt/dx)*(u[i,j]*(ro[i+1,j]+ro[i,j])-u[i-1,j]*(ro[i-1,j]+ro[i,j]))-(0.5* dt/dy)*(v[i,j]*(ro[i,j+1]+ro[i,j])-v[i,j-1]*(ro[i,j-1]+ro[i,j]))+(reynolds*dt/dx/dx)*(ro[i+1,j]-2.0*ro[i,j]+ro[i-1,j])+(reynolds*dt/dy/dy)*(ro[i,j+1]-2.0*ro[i,j]+ro[i,j-1])
        end
    end
    time = time+dt

    uu = 0.5 * (u[:, 2:ny+2] + u[:, 1:ny+1])
    vv = 0.5 * (v[2:nx+2, :] + v[1:nx+1, :])

    p1 = heatmap(
        x[2:end-1],
        y[2:end-1],
        r'[2:end-1,2:end-1],
        c = :blues,
        clim = (r1, r2),
        # aspect_ratio = :equal,
        # title = "Density",
        # legend = false
    )
    # quiver!(
    # xhm,
    # yhm,
    # quiver = (reshape(u, (1, xSize^2)), reshape(v, (1, ySize^2))), # Be careful here to take into account how arrays work vs how heatmaps are plotted.
    # aspect_ratio = :equal,
    # # title = "Velocity field",
    # arrow = arrow(0.01, 0.01)
    # )
    plot(p1, size = (1000, 1000))

    println(steps/nstep*100)
end

gif(anim, "./variable_density.gif", fps = 50)

# fig = Figure(resolution = (1200, 900))
# ax, hm = heatmap(fig[1,1][1,1], x, y, r) #log.(sqrt.(uu.^2 + vv.^2))
# arrows!(xh, yh, uu, vv, lengthscale = 0.1, arrowsize = 0.01)
# Colorbar(fig[1,1][1,2], hm, width = 25, colormap = :viridis)
# fig
