<comment>
problem = jeans problem (standing wave; stable or unstable)
author  = E.C. Ostriker
journal =
config  = --with-problem=jeans --with-gravity=fft --enable-fft

<job>
problem_id      = Jeans     # problem ID: basename of output filenames
maxout          = 4         # Output blocks number from 1 -> maxout
num_domains     = 1         # number of Domains in Mesh

<output1>
out_fmt = hst               # History data dump
dt      = 0.0001              # time increment between outputs

<output2>
out_fmt = vtk               # Binary data dump
dt      = 0.0001
              # time increment between outputs

<output3>
out_fmt = tab
out     = d
id      = d
dt      = 0.0001             # time increment between outputs

<output4>
out_fmt = ppm
out     = d
id      = d
dt      = 0.0001             # time increment between outputs
x3      = 0.0


<time>
cour_no         = 0.5 #for 3D # The Courant, Friedrichs, & Lewy (CFL) Number
nlim            = 10000    # cycle limit
tlim            = 10.0       # time limit (units lambda/ciso)

<domain1>
level           = 0         # refinement level this Domain (root=0)
Nx1             = 128        # Number of zones in X-direction
x1min           = -4*PI       # minimum value of X
x1max           = 4*PI      # maximum value of X
bc_ix1          = 4         # boundary condition flag for inner-I (X1)
bc_ox1          = 4         # boundary condition flag for outer-I (X1)

Nx2             = 128       # Number of zones in X2-direction
x2min           = -4*PI       # minimum value of X2
x2max           = 4*PI       # maximum value of X2
bc_ix2          = 4         # boundary condition flag for inner-J (X2)
bc_ox2          = 4         # boundary condition flag for outer-J (X2)

Nx3             = 128        # Number of zones in X3-direction
x3min           = -4*PI       # minimum value of X3
x3max           = 4*PI       # maximum value of X3
bc_ix3          = 4         # boundary condition flag for inner-K (X3)
bc_ox3          = 4         # boundary condition flag for outer-K (X3)

AutoWithNroc    = 4

<problem>
iso_csound      = 1.0	    # isothermal sound speed (used if ISOTHERMAL)
gamma = 1.6666666666666667  # gamma = C_p/C_v
amp             = 0.0    # Initial wave amplitude
ampb            = 0.0       #P Wave amplitude of B perturbation
beta		= 4.0	    # Pth/Pmag (if MHD is defined)
njeans		= 4.0	    # lambda/lambda_jeans
kdir		= 2	    # direction of wavenumber (1, 2, or 3)
