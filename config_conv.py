import sys
# problem type      0           1           2           3               4               5
problem_list = ["gaussian", "ep-shock", "v-shock", "v-gaussian", "ep-shocktube", "ep-v-gaussian"]
problem = problem_list[0]

# scheme
scheme_list =  ["kt-minmod-tvdrk2", "kt-superbee-tvdrk2", "kt-vanleer-tvdrk2", "kt-mc-tvdrk2"]
scheme = scheme_list[0]

# mesh
# ==========================
x0 = True # include x = 0 in grid?
cCFL = 0.5 # courant factor <=0.5

hx = 0.05/2 # step sizes (0.2,0.1,0.05,0.025,0.0125,0.0125/2,0.0125/4,0.0125/8)
ht = 0.5*hx # for now make this a fixed (SMALL) number independent of hx so that time steps line up in convergence test
if ht > cCFL*hx:
    print("CFL Condition not met. Exiting program.")
    sys.exit()

N = 75  # length of grid
iterations = N/(ht) 
# tmod = 200 # int(1/ht)
tmod = int(1/ht)
# ==========================

# viscosity
# let's test {1,3,6}
Neta = 1 # multiple of 1/4pi for eta/s

# conformal
ep_coeff = 5.26379 # alpha in ep = alpha * T^4

# hydro frames
# ==========================
choose_cplus_a1 = False
frame = 2
    # PAPER FRAMES
    # =============
    # Teaney/Princeton frames
    # c_+ = 0.85, a1 = 25/2 = 12.5, a2 = 25/3 ~ 8.3 -> A
    # c_+ = 0.99, a1 = 25/4 = 6.25 , a2 = 25/7 ~ 3.6 -> B
    # OURS
    # c_+ = 0.99, a1 = 5, a2 = 3.97 -> C?

if choose_cplus_a1 == True: 

    if frame == 1: # B - Princeton
        cplus = 0.9999
        a1 = 25/4 
        # a2 =~ 25/7
    
    if frame == 2: # A - Princeton
        cplus = 0.85
        a1 = 25/2
        # a2 = 25/3

    if frame == 3: # pick something
        cplus = 0.74
        a1 = 25
        # a2 =~ 25

else:

    if frame == 1: # B - Princeton
        # cplus = 1.0
        a1 = 25/4 
        a2 = 25/7
    
    if frame == 2: # A - Princeton
        # cplus = 0.85
        a1 = 25/2
        a2 = 25/3

    if frame == 3: # pick something
        # cplus = 0.74
        a1 = 25
        a2 = 25

    # if frame == 4:
    #     a1 = 5
    #     a2 = 5
    # if frame == 5:
    #     a1 = 10
    #     a2 = 5
    # if frame == 6:
    #     a1 = 20
    #     a2 = 5

    # if frame == 5: # small c_+
    #     a1 = 100
    #     a2 = 100


 
# ==========================

# Kurganov-Tadmor
theta_kt = 1.0

# metadata 
# ==========================
# alarms on or off?
alarms = False
# automatically pass checks and prints and generate file names for sucomp cluster/speed
sucomp = True
# animation bounds?
anim_bounds = True
# zoom in on hydro frame animations?
anim_zoom = False
# ==========================

# Initial conditions
# ==========================
if problem == problem_list[0]:
    epL = 0.48 # epLarger
    epR = 0.12 #epsmalleR
    vL = 0.0
    vR = vL
    Delta = 5.0 # width, WILL BE SQUARED

# smooth shocktube epsilon
if problem == problem_list[1]:
    epL = 0.48
    epR = 0.12
    vL = 0.0
    vR = vL
    Delta = 1.0 # 0.2 

# smooth shocktube velocity
if problem == problem_list[2]:
    # constant background eps
    epL = 0.5
    epR = epL
    # L/R distribution in v
    vL = 0.90
    vR = 0.45
    # distribution gradient width
    Delta = 1.0

# negative gaussian in velocity
if problem == problem_list[3]:
    epL = 0.5
    epR = epL
    vL = 0.6 # max
    vR = 0.3 # height below vL
    Delta = 5.0

# shocktube in epsilon
if problem == problem_list[4]:
    epL = 0.48
    epR = 0.12
    vL = 0.0
    vR = vL
    Delta = 0.0

# gaussian in epsilon and velocity
if problem == problem_list[5]:
    epL = 0.48 # epLarger
    epR = 0.12 #epsmalleR
    vL = 0.3
    vR = 0.0
    Delta = 5.0 # width, WILL BE SQUARED

# ==========================
