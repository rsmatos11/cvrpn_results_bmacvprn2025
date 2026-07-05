import ast
import math
import os
import sys
from docplex.mp.model import Model
from mip_start_cplex import setup_initial_solution_cvrpn

"""
Script for solving the CVRPN with the CPLEX solver and the DOcplex Python API.

The script supports two formulations:
- Model (2)--(20)
- Model (2)--(26) with additional constraints
"""

#===================== MODEL CONFIGURATION ====================#
# False: Model (2)--(20); True : Model (2)--(26)
use_additional_constraints = False
#==============================================================#

#================ INITIAL SOLUTION CONFIGURATION ==============#
# If True, an initial solution file (.txt) must be provided
# as the second command-line argument.
use_initial_solution = False

# If True, the route file (.txt) must be provided as the
# second command-line argument. The routing variables are
# fixed and only the visit points are optimized.
use_fixed_initial_routes = False
#==============================================================#

# Prevent inconsistent configuration
if use_initial_solution and use_fixed_initial_routes:
    print("Error: choose either use_initial_solution or use_fixed_initial_routes, not both.")
    sys.exit(1)

# ----------------------------------------------------------
# Check command-line arguments
#
# Usage:
#   python cvrpn_models_cplex.py instance.dat
#   python cvrpn_models_cplex.py instance.dat initial_solution.txt
# ----------------------------------------------------------

if len(sys.argv) < 2:
    print("Usage:")
    print("  python cvrpn_models_cplex.py instance_file.dat")
    print("  python cvrpn_models_cplex.py instance_file.dat initial_solution.txt")
    sys.exit(1)

# Store the instance filename
filename = sys.argv[1]
print(f"Reading data from file: {filename}")

# Base instance name without file extension
name_inst = os.path.splitext(os.path.basename(filename))[0]

# Optional initial solution file (.txt)
initial_solution_path = sys.argv[2] if len(sys.argv) >= 3 else None

#================= INSTANCE DATA LOADING =================#
def load_dat(filename,variables_needed):
    with open(filename,'r') as file:
        content = file.read().replace('{', '[').replace('}', ']').replace(';', '')
    variables = {}
    for line in content.splitlines():
        line = line.strip()
        if '=' not in line or not line:
            continue
        key, value = line.split('=',1)
        key = key.strip()
        if key in variables_needed:
            try:
                variables[key] = ast.literal_eval(value.strip())
            except (ValueError, SyntaxError) as e:
                print(f"Error while loading variable {key}: {str(e)}")
    return variables

# Variables required for the CVRPN model
var_cvrpn = ['m','maxcap','dem','depot','scalefactor','n','h','U','xp','yp']

# Load data from the .dat instance file
data = load_dat(filename,var_cvrpn)

# If the instance was not scaled, assume scalefactor = 1.0
data.setdefault('scalefactor',1.0)

# Check if any required variable is missing
missing_vars = [var for var in var_cvrpn if var not in data]
if missing_vars:
    raise ValueError(f"Missing variables in instance file: {', '.join(missing_vars)}")

#================= ASSIGNMENT OF INSTANCE DATA =================#
# Number of clients (polygons)
m = data['m']
# Vehicle capacity
maxcap = float(data['maxcap'])
# Demand of each client
dem = [float(d) for d in data['dem']]
# Original depot coordinates
depot_orig = tuple(float(v) for v in data['depot'])
# Scale factor
scalefactor = float(data['scalefactor'])
# Total number of polygon vertices
n = data['n']
# Number of convex parts per client
h = data['h']
# Convex decomposition of each polygon
U = [[sorted(list(part)) for part in poly] for poly in data['U']]
# X and Y coordinates of all polygon vertices
xp = data['xp']
yp = data['yp']

#================= COORDINATE SYSTEM SHIFT =================#
# Shift the coordinate system so that the depot is located at the origin
dptx,dpty = depot_orig
# Relative coordinates of vertices with respect to the depot
xp_rel = [x - dptx for x in xp]
yp_rel = [y - dpty for y in yp]
# Depot coordinates after the shift
depot_rel = (0.0,0.0)

# ====== Function to compute coefficients a, b, and c ======
def calc_coefficients(m,h,xp_rel,yp_rel,U):
    a,b,c = {},{},{}
    for i in range(m):
        for k in range(h[i]):
            vertices = U[i][k]
            nv = len(vertices)
            for idx, j in enumerate(vertices):
                j_next = vertices[(idx + 1) % nv]
                aj = yp_rel[j_next] - yp_rel[j]
                bj = xp_rel[j] - xp_rel[j_next]
                cj = xp_rel[j] * yp_rel[j_next] - yp_rel[j] * xp_rel[j_next]
                norm = abs(aj) + abs(bj) + abs(cj)
                if norm > 0:
                    aj /= norm
                    bj /= norm
                    cj /= norm
                a[j],b[j],c[j] = aj,bj,cj
    return a,b,c

# Compute the coefficients using the shifted coordinates
a,b,c = calc_coefficients(m,h,xp_rel,yp_rel,U)

# ===== Structure of arcs for the decision variables =====
# Create the list of feasible arcs (i,j) such that i != j
clients = list(range(m + 1))  # 0 = depot, 1..m = polygons
arcs = [(i,j) for i in clients for j in clients if i != j]

# Set of convex parts obtained from the decomposition of non-convex polygons
convex_parts = [(i,l) for i in range(1,m + 1) for l in range(1,h[i - 1] + 1) if h[i - 1] > 1]

#==========================================================#
#               CVRPN CPLEX MODEL - MIQCP                  #
#==========================================================#

# Define the CPLEX optimization model for the CVRPN problem
mdl = Model(name="CVRPN_cplex")

# Time limit parameter (1 hour)
mdl.context.cplex_parameters.timelimit = 3600

# Number of threads used by CPLEX (1 = single-thread; comment out for default automatic setting).
mdl.parameters.threads = 1

# Big-M values
# 1. Compute the bounding box limits
xmin,xmax = min(xp_rel + [0.0]),max(xp_rel + [0.0])
ymin,ymax = min(yp_rel + [0.0]),max(yp_rel + [0.0])

# 2. Compute M1 as the Euclidean distance across the bounding box diagonal (upper bound on distance)
M1 = math.sqrt((xmax - xmin)**2 + (ymax - ymin)**2)

# 3. Compute M2 by evaluating the linear functions at the 4 extreme vertices of the bounding box
corners = [(xmin,ymin),(xmin,ymax),(xmax,ymin),(xmax,ymax)]
M2 = max(a[j]*cx + b[j]*cy + c[j] for j in a for cx,cy in corners)

# Depot coordinates (at the origin)
dptx,dpty = depot_rel

#============ MODEL DECISION VARIABLES ============#
x = mdl.binary_var_dict(arcs,name="x")                                                      # 1 if arc (i,j) is used, 0 otherwise
z = mdl.continuous_var_dict(arcs,lb=0.0,name="z")                                           # Load flow on arc (i,j)
px = mdl.continuous_var_dict(range(1, m+1), lb=-float('inf'), ub=float('inf'), name="px")   # X-coordinate of the visit point for polygon i
py = mdl.continuous_var_dict(range(1, m+1), lb=-float('inf'), ub=float('inf'), name="py")   # Y-coordinate of the visit point for polygon i
gama = mdl.continuous_var_dict(arcs, lb=0.0,name="gama")                                    # Euclidean distance associated with arc (i,j)
y = mdl.continuous_var_dict(arcs,lb=0.0,name="y")                                           # Auxiliary variable for the quadratic norm (\xi)
w = mdl.binary_var_dict(convex_parts,name="w")                                              # Activation of convex parts
# Variables for the additional constraints [Model (2)--(26)]
if use_additional_constraints:
    qx = mdl.continuous_var_dict(arcs,lb=0.0,name="qx")                                     # Auxiliary variable \xi_x
    qy = mdl.continuous_var_dict(arcs,lb=0.0,name="qy")                                     # Auxiliary variable \xi_y
#===================================================#

#=================== OBJECTIVE FUNCTION ===================#

mdl.minimize( mdl.sum(gama[i,j] for (i,j) in arcs) )

#======================== CONSTRAINTS =====================#

# Ensure that the number of arcs leaving the depot is equal to the number of arcs entering the depot
mdl.add_constraint( mdl.sum(x[0,j] for j in range(1,m + 1)) == mdl.sum(x[i,0] for i in range(1,m + 1)) )

# In-degree and out-degree constraints for each client
for i in range(1,m + 1):
    mdl.add_constraint( mdl.sum(x[j,i] for j in clients if j != i) == 1 )
    mdl.add_constraint( mdl.sum(x[i,j] for j in clients if j != i) == 1 )

# Load flow balance constraints at each client
for k in range(1,m + 1):
    mdl.add_constraint(
        mdl.sum(z[i,k] for i in clients if i != k) -
        mdl.sum(z[k,j] for j in clients if j != k) == dem[k - 1]
    )

# Capacity limit on each arc
for (i,j) in arcs:
    mdl.add_constraint( z[i,j] <= maxcap * x[i,j] )

# Total load flow leaving the depot must be equal to the total demand
mdl.add_constraint( mdl.sum(z[0,j] for j in range(1,m + 1)) == sum(dem) )

# No load flow enters the depot
mdl.add_constraint( mdl.sum(z[i,0] for i in range(1,m + 1)) == 0 )

# Euclidean distance constraints and linking constraints with y and gama
for (i,j) in arcs:
    pxi = dptx if i == 0 else px[i]
    pyi = dpty if i == 0 else py[i]
    pxj = dptx if j == 0 else px[j]
    pyj = dpty if j == 0 else py[j]

    # Quadratic constraint
    mdl.add_constraint( y[i,j] * y[i,j] >= (pxi - pxj) * (pxi - pxj) + (pyi - pyj) * (pyi - pyj) )

    mdl.add_constraint( gama[i,j] <= M1 * x[i,j] )
    mdl.add_constraint( gama[i,j] >= y[i,j] - M1 * (1 - x[i,j]) )

    # Additional constraints for Model (2)--(26)
    if use_additional_constraints:
        mdl.add_constraint( pxi - pxj <=  qx[i,j] )
        mdl.add_constraint( pxi - pxj >= -qx[i,j] )

        mdl.add_constraint( pyi - pyj <=  qy[i,j] )
        mdl.add_constraint( pyi - pyj >= -qy[i,j] )

        mdl.add_constraint( y[i,j] >= qx[i,j] )
        mdl.add_constraint( y[i,j] >= qy[i,j] )
        mdl.add_constraint( y[i,j] <= qx[i,j] + qy[i,j] )
    

# Geometric constraints for convex polygons
for i in range(1,m + 1):
    if h[i - 1] == 1:
        for k in U[i - 1][0]:
            mdl.add_constraint( a[k] * px[i] + b[k] * py[i] <= c[k] )

# Geometric constraints for non-convex polygons
for i in range(1,m + 1):
    if h[i - 1] > 1:
        mdl.add_constraint( mdl.sum(w[i,l] for l in range(1,h[i - 1] + 1)) == 1 )  # sum of w[i,l] must be equal to 1
        for l in range(1,h[i - 1] + 1):
            for k in U[i - 1][l - 1]:
                mdl.add_constraint( a[k] * px[i] + b[k] * py[i] <= c[k] + M2 * (1 - w[i,l]) )

#===========================================#
#               END OF MODEL                #
#===========================================#

print()
print('=====================================')

# Print general information about the model
print("Problem: CVRPN")

# Print the model type reported by CPLEX
print("Model type:",mdl.problem_type)

if use_additional_constraints:
    print("Formulation: Model (2)--(26)")
else:
    print("Formulation: Model (2)--(20)")
    
print('scalefactor=',scalefactor)

# Indicate whether a warm-start initial solution is used
if use_initial_solution:
    print('USING initial solution (warm start)')

# Indicate whether fixed initial routes are used
if use_fixed_initial_routes:
    print('USING initial solution with FIXED routes: only optimizing visit points')

print('=====================================')
print()

#============ INITIAL SOLUTION ============#
# Load and apply a warm-start solution if requested
if use_initial_solution:

    if initial_solution_path is None:
        print("Error: use_initial_solution = True, but no initial solution file was provided.")
        print("Usage:")
        print("  python cvrpn_models_cplex.py instance_file.dat initial_solution.txt")
        sys.exit(1)

    if os.path.exists(initial_solution_path):
        heuristic,sol_ini,nroutes_ini = setup_initial_solution_cvrpn(
            filename=initial_solution_path,
            mdl=mdl,
            x=x,
            w=w,
            h=h,
            U=U,
            xp_rel=xp_rel,
            yp_rel=yp_rel,
            x_dep=depot_orig[0],
            y_dep=depot_orig[1],
            fix_routes=False
        )
    else:
        print(f"Error: initial solution file '{initial_solution_path}' not found.")
        sys.exit(1)

#========= FIXED INITIAL ROUTES =========#
# Load and apply fixed routes.
# Routing variables are fixed, and the solver only optimizes the visit points inside each polygon.
elif use_fixed_initial_routes:

    if initial_solution_path is None:
        print("Error: use_fixed_initial_routes = True, but no route file was provided.")
        print("Usage:")
        print("  python cvrpn_models_cplex.py instance_file.dat initial_solution.txt")
        sys.exit(1)

    if os.path.exists(initial_solution_path):
        heuristic,sol_ini,nroutes_ini = setup_initial_solution_cvrpn(
            filename=initial_solution_path,
            mdl=mdl,
            x=x,
            w=w,
            h=h,
            U=U,
            xp_rel=xp_rel,
            yp_rel=yp_rel,
            x_dep=depot_orig[0],
            y_dep=depot_orig[1],
            fix_routes=True
        )
    else:
        print(f"Error: route file '{initial_solution_path}' not found.")
        sys.exit(1)

#======= SOLVE THE MODEL =======#
try:
    solution = mdl.solve(log_output=True)
except Exception as erro:
    print("Error during optimization:",erro)
    solution = mdl.solution

# Retrieve solver details
details = mdl.solve_details

print()
print('===========================================================')
print()

# Solution information
if solution:
    print(f"Solution found (scaled): {solution.objective_value:.12f}")
else:
    print("No solution found.")

# Always print solver details
print("\n-------- SOLUTION STATUS --------")
print("Status:", details.status)
print("Total time (s):", details.time)
print("Best bound (scaled):", details.best_bound if details else "-")
print("Gap:", details.mip_relative_gap if details else "-")

print()
print('===========================================================')
print()

#==========================================================#
#                     POST-PROCESSING                      #
#==========================================================#

nroutes = ""

print("\n********** POST-PROCESSING **********")

if solution and solution.objective_value is not None:
    print("\n1. Solution data (original scale):\n")

    print(f"Objective value: {solution.objective_value * scalefactor:.12f}")
    print(f"Best bound: {details.best_bound * scalefactor:.12f}")
    print(f"Gap_percent: {details.mip_relative_gap * 100:.12f}")
    print(f"Elapsed_time: {details.time:.12f}")
    print(f"Number of iterations: {details.nb_iterations}")
    print(f"Total number of B&B nodes: {details.nb_nodes_processed}")

    nroutes = sum(1 for j in range(1,m+1) if x[0,j].solution_value > 0.5)
    print(f"Number of routes: {nroutes}")

    print("\n2. Instance data:\n")
    print(f"Total number of variables: {mdl.number_of_variables}")
    print(f"Binary variables: {mdl.number_of_binary_variables}")
    print(f"Integer variables: {mdl.number_of_integer_variables}")
    print(f"Continuous variables: {mdl.number_of_continuous_variables}")
    print(f"Linear constraints: {mdl.number_of_linear_constraints}")
    print(f"Quadratic constraints: {mdl.number_of_quadratic_constraints}")

    print("\n3. Coordinates of the visit points in each polygon:\n")
    for i in range(1,m+1):
        print(f"Polygon {i}: px = {px[i].solution_value * scalefactor:.15f}, py = {py[i].solution_value * scalefactor:.15f}")

    print("\n4. Route(s) found:\n")
    routes = []
    visited = set()

    for j in range(1,m+1):
        if (0,j) in x and x[0,j].solution_value > 0.5 and j not in visited:
            route = [0,j]
            current = j
            visited.add(j)

            while current != 0:
                for k in clients:
                    if current != k and (current,k) in x and x[current,k].solution_value > 0.5:
                        route.append(k)
                        current = k
                        if k != 0:
                            visited.add(k)
                        break
                else:
                    route.append(0)
                    current = 0
                    break

            length = sum(gama[i,j].solution_value for i, j in zip(route[:-1],route[1:])) * scalefactor
            demand = sum(dem[j-1] for j in route[1:-1])
            routes.append((route,length,demand))

    # Print the routes
    for idx,(route,length,demand) in enumerate(routes):
        print(f"Route {idx+1}: {' -> '.join(map(str, route))}")
        print(f"Length:{length:.12f}")
        print(f"Total demand:{demand}")

    print("\n*************** END OF POST-PROCESSING ***************\n")

else:
    print("The model did not find an optimal or feasible solution, or the time limit was reached.")
