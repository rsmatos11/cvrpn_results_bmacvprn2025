import math
import os

"""
MIP Start and Route Fixing Utilities for CVRPN Models.

This module provides automation functions to inject initial heuristic solutions 
(e.g., Clarke-Wright or Nearest Neighbor) into CPLEX models via DOcplex.

Available Approaches:
1. Warmstart (MIP Start): Injects initial values for decision variables (x,w) 
   to guide the solver, allowing CPLEX to modify the routes if better ones are found.
2. Fixed Routes: Strictly locks the routing variables (x), forcing CPLEX to 
   focus exclusively on optimizing continuous target coordinates (px,py).
"""

#=======================================================================#
#                 INITIAL SOLUTION SETUP FOR CVRPN                      #
#=======================================================================#

def setup_initial_solution_cvrpn(filename,mdl,x,w,h,U,xp_rel,yp_rel,x_dep,y_dep,fix_routes=False):
    """
    Reads a heuristic solution file (CW or NN) and applies it to the CVRPN model.
    
    If fix_routes=True:  Strictly fixes/locks the routing variables (x) in the model.
    If fix_routes=False: Provides a MIP start (warmstart) for x and w variables,
                         allowing CPLEX to modify them during optimization.
    """

    if not os.path.exists(filename):
        print(f"[ERROR] File {filename} not found!")
        return None,None,0

    # 1. READ THE SOLUTION FILE
    with open(filename,"r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # 2. EXTRACT HEURISTIC NAME AND OBJECTIVE VALUE
    heuristic_name = ""
    objective_value = None
    
    for line in lines:
        if "CW SOLUTION" in line:
            heuristic_name = "CW"
        elif "NN SOLUTION" in line:
            heuristic_name = "NN"

        if "Objective:" in line:
            parts = line.split(":")
            if len(parts) >= 2:
                try:
                    objective_value = float(parts[1].strip())
                except ValueError:
                    pass

        if heuristic_name and objective_value is not None:
            break

    # 3. EXTRACT ROUTES AND GEOMETRIC COORDINATES (PATHCOO)
    routes = []
    all_coords = []
    i = 0
    while i < len(lines):
        line = lines[i]
        
        if line.startswith("route ="):
            clusters = list(map(int, line.split("=")[1].split()))
            route = [0] + clusters + [0]  # Add depot (0) at the start and end
            routes.append(route)
            
            while i < len(lines) and not lines[i].startswith("pathcoo ="):
                i += 1
                
            if i < len(lines) and lines[i].startswith("pathcoo ="):
                path_vals = list(map(float, lines[i].split("=")[1].split()))
                coords = [(path_vals[j], path_vals[j+1]) for j in range(0, len(path_vals), 2)]
                all_coords.append(coords)
        i += 1

    nroutes_ini = len(routes)
    
    # ==============================================================
    # INFORMATIVE LOGGING (ONLY DISPLAYED IN ROUTE-FIXING MODE)
    # ==============================================================
    if fix_routes:
        print('==================================================')
        print("Detailed Fixed Routes Applied to Model:")
        for r, route in enumerate(routes, 1):
            print(f"  Route {r}: {' -> '.join(map(str, route))}")
        print(f"  Heuristic: {heuristic_name}")
        if objective_value is not None:
            print(f"  Initial objective value (scaled): {objective_value:.12f}")
        print(f"  Initial number of routes: {nroutes_ini}")
        print('==================================================\n')

    # AUXILIARY INTERNAL FUNCTION (Ray-Casting Algorithm for Warmstart)
    def point_in_polygon(x_p,y_p,polygon,tol=1e-6):
        n = len(polygon)
        if n < 3: return False
        inside = False
        for idx in range(n):
            x1,y1 = polygon[idx]
            x2,y2 = polygon[(idx + 1) % n]
            if abs(x_p - x1) < tol and abs(y_p - y1) < tol:
                return True
            area = abs((x2 - x1) * (y_p - y1) - (x_p - x1) * (y2 - y1))
            if area < tol:
                if (min(x1,x2) - tol <= x_p <= max(x1,x2) + tol and 
                    min(y1,y2) - tol <= y_p <= max(y1,y2) + tol):
                    return True
        p1x,p1y = polygon[0]
        for idx in range(1, n + 1):
            p2x,p2y = polygon[idx % n]
            if y_p > min(p1y,p2y) and y_p <= max(p1y,p2y) and x_p <= max(p1x,p2x):
                if p1y != p2y:
                    xinters = (y_p - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or x_p <= xinters:
                    inside = not inside
            p1x,p1y = p2x,p2y
        return inside

    # 4. APPLY THE SELECTED STRATEGY TO THE CPLEX MODEL
    
    # STRATEGY A: DEFINITIVELY LOCK/FIX ROUTES (set_lb and set_ub)
    if fix_routes:
        for route in routes:
            for k in range(len(route) - 1):
                a_node,b_node = route[k],route[k + 1]
                if (a_node,b_node) in x:
                    x[a_node,b_node].set_lb(1)
                    x[a_node,b_node].set_ub(1)
                    
    # STRATEGY B: INJECT WARMSTART VALUES (MIP Start)
    else:
        mip_start = mdl.new_solution()
        x_activated = set()
        w_activated = set()

        for route_idx, route in enumerate(routes):
            coords = all_coords[route_idx] if route_idx < len(all_coords) else []
            
            # Activate decision variables 'x' from the heuristic solution
            for k in range(len(route) - 1):
                a_node,b_node = route[k],route[k + 1]
                if (a_node,b_node) in x:
                    mip_start.add_var_value(x[(a_node,b_node)],1)
                    x_activated.add((a_node,b_node))
            
            # Map cluster partition binary variables 'w' based on spatial coordinates
            for point_idx, cluster_i in enumerate(route[1:-1]):
                if point_idx < len(coords):
                    # Translate coordinates from global system to relative system (Depot as Origin)
                    x_val = coords[point_idx][0] - x_dep
                    y_val = coords[point_idx][1] - y_dep
                    
                    if cluster_i > 0 and h[cluster_i - 1] > 1:
                        point_found = False
                        convex_parts = U[cluster_i - 1]
                        
                        # Identify which convex part contains the continuous location point
                        for l_idx, polygon_vertices in enumerate(convex_parts,1):
                            polygon_points = [(xp_rel[v],yp_rel[v]) for v in polygon_vertices]
                            
                            if point_in_polygon(x_val,y_val,polygon_points):
                                if (cluster_i,l_idx) in w:
                                    mip_start.add_var_value(w[(cluster_i,l_idx)],1)
                                    w_activated.add((cluster_i,l_idx))
                                    point_found = True
                                    break
                        
                        # Fallback case: assign to the first convex part if ray-casting misses due to precision
                        if not point_found and convex_parts and (cluster_i,1) in w:
                            mip_start.add_var_value(w[(cluster_i,1)],1)
                            w_activated.add((cluster_i,1))

        # Explicitly set all remaining binary variables to 0 to ensure MIP Start consistency
        for (i,j) in x:
            if (i,j) not in x_activated:
                mip_start.add_var_value(x[(i,j)],0)
        for (i,l) in w:
            if (i,l) not in w_activated:
                mip_start.add_var_value(w[(i,l)],0)

        # Inject the constructed solution into the mathematical model
        if mip_start.number_of_var_values > 0:
            mdl.add_mip_start(mip_start)
        else:
            print("[WARNING] No variables activated in starting solution.")

    return heuristic_name,objective_value,nroutes_ini
