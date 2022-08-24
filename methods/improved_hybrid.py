#This code is implements a hybrid differencing scheme to the cells

#it applies CD to the boundary cells
#while it uses a hybrid differencing (CD if Pe < 2, diffusion dominates, UD if Pe >2, convection dominates) 
#in the interior

#this method will give rise to a diagonally dominant matrix iff diffusion is stronger than convection
#quantitatovely, if Pe < 2 in the greatest part of the domain

#particularly problematic are the boundaries where CD is always applied, this is the main cause of divergence of 
#Gauss-Seidel when diffusion is weak (see problem with water)

#when convection is much stronger than diffusion it is a good idea to use a full UD scheme, applied also to the boundaries

#ATTENTION TO PROBLEM GEOMETRY!!!!
#this code is writttern referring to the circular domain problem
#it is a general purpose code in the limits explained above, if the variables
#location_cell, the position of the cell in the geometry according to the sol_mat data structure
#and
# E_cell_U, W_cell_U, S_cell_V, N_cell_V are redefined for the new geometry (STAGGERED VELOCITY ARRANGEMENT)

#we say that this hybrid is imporved because we have used more flexible logical statements that can account for
#smooth reasonably whirling flows

#SPECIFIC TO CIRCULAR DOMAIN GEOMETRY - referred as mug here!
#define a vector of the number of cell difference between adjacent rows, needed for the bottom half
row_diff_vect = np.zeros((1, int(N_d/2)));
for i in range(0, int(N_d/2)):
    diff = cells_per_row[i + 1] - cells_per_row[i];
    row_diff_vect[0, i] = diff;
#half of it    
half_row_diff_vect = row_diff_vect/2;
#we are going to use these two vectors to compute the locations of the V velocities wrt the scalar CV

number_of_cell = 0;

for a_y in range(0, N_d):
    for a_x in range(0, N_d):

        #THIS PART MUST BE REDEFINED ACCORDING TO THE PROBLEM GEOMETRY
        #####################################################################################################
        num_0_left = (N_d - cells_per_row[a_y + 1])/2;
        location_cell = int(sum(cells_per_row[1 : (a_y + 1)]) + a_x - num_0_left)
        
        #the cell number is accounted for in this more efficient way (count the cells inside the 1s domain)
        register = int(location_cell*5);
        S = int(register);
        W = register + 1;
        P = register + 2;
        E = register + 3;
        N = register + 4;
        ###################################################################################################### 
            
        #define the velocity vectors that will surround the cells inside the itaration algorithm
        #first get rid of the BC accounting for the individual case
        #then we have the general algorithm for the remaining cells
        
        #ATTENTION THAT THE DOMAIN IS UPSIDE-DOWN!
        #first of all we account for the boundary cells using CD
        if mug[a_y + 1, a_x + 1] == 0:
            continue
        #Boundaries
        elif mug[a_y + 1, a_x + 1] == 1:
            
            number_of_cell = number_of_cell + 1;
            
            Pe_u = 0;
            Pe_v = 0;
            
            #THIS PART MUST BE REDEFINED ACCORDING TO THE PROBLEM GEOMETRY
            ####################################################################################################
            #define the velocity cells surrounding the scalar CV
            W_cell_U = int(sol_mat[W, 1] - a_y); ###### correction for W and E to be carried out to all lines ######
            E_cell_U = int(sol_mat[E, 1] - a_y - 1);
            #for the S and N velocity cells we need a if-statement
        
            if a_y < N_d/2:
                N_cell_V = int(sol_mat[P, 1]);
                S_cell_V = int(sol_mat[S, 1]);
            elif a_y == N_d/2: #exception
                t = int(a_y - N_d/2)
                N_cell_V = int(sol_mat[P, 1] - sum(row_diff_vect[0, 0 : t]) - half_row_diff_vect[0, t]);
                S_cell_V = int(sol_mat[S, 1]);
            elif a_y > N_d/2:
                t = int(a_y - N_d/2)
                N_cell_V = int(sol_mat[P, 1] - sum(row_diff_vect[0, 0 : t]) - half_row_diff_vect[0, t]);
                S_cell_V = int(sol_mat[S, 1] - sum(row_diff_vect[0, 0 : t - 1]) - half_row_diff_vect[0, t - 1])
            ######################################################################################################
            
            if (mug[a_y + 1, a_x] == 0 and mug[a_y, a_x + 1] == 0): #S, W corner
                sol_mat[S, 2] = 0;
                sol_mat[W, 2] = 0;
                sol_mat[P, 2] = (F_u[0, E_cell_U]/2 + F_v[0, N_cell_V]/2 + 4*D);
                sol_mat[E, 2] = F_u[0, E_cell_U]/2 - D;
                sol_mat[N, 2] = F_v[0, N_cell_V]/2 - D;
            elif (mug[a_y + 1, a_x + 2] == 0 and mug[a_y, a_x + 1] == 0): #S, E corner
                sol_mat[S, 2] = 0;
                sol_mat[W, 2] = - F_u[0, W_cell_U]/2 - D;
                sol_mat[P, 2] = (F_v[0, N_cell_V]/2 - F_u[0, W_cell_U]/2 + 4*D);
                sol_mat[E, 2] = 0;
                sol_mat[N, 2] = F_v[0, N_cell_V]/2 - D;
            elif (mug[a_y + 1, a_x] == 0 and mug[a_y + 2, a_x + 1] == 0): #N, W corner    
                sol_mat[S, 2] = - F_v[0, S_cell_V]/2 - D;
                sol_mat[W, 2] = 0;
                sol_mat[P, 2] = (F_u[0, E_cell_U]/2 - F_v[0, S_cell_V]/2 + 4*D)
                sol_mat[E, 2] = F_u[0, E_cell_U]/2 - D;
                sol_mat[N, 2] = 0;
            elif (mug[a_y + 1, a_x + 2] == 0 and mug[a_y + 2, a_x + 1] == 0): #N, E corner
                sol_mat[S, 2] = - F_v[0, S_cell_V]/2 - D;
                sol_mat[W, 2] = - F_u[0, W_cell_U]/2 - D;
                sol_mat[P, 2] = (- F_u[0, W_cell_U]/2 - F_v[0, S_cell_V]/2 + 4*D);
                sol_mat[E, 2] = 0;
                sol_mat[N, 2] = 0;
            elif mug[a_y + 1, a_x] == 0: #W wall
                sol_mat[S, 2] = - F_v[0, S_cell_V]/2 - D;
                sol_mat[W, 2] = 0;
                sol_mat[P, 2] = (F_u[0, E_cell_U]/2 + F_v[0, N_cell_V]/2 - F_v[0, S_cell_V]/2) + 4*D;
                sol_mat[E, 2] = F_u[0, E_cell_U]/2 - D;
                sol_mat[N, 2] = F_v[0, N_cell_V]/2 - D;
            elif mug[a_y + 1, a_x + 2] == 0: #E wall
                sol_mat[S, 2] = - F_v[0, S_cell_V]/2 - D;
                sol_mat[W, 2] = - F_u[0, W_cell_U]/2 - D;
                sol_mat[P, 2] = (- F_u[0, W_cell_U]/2 - F_v[0, N_cell_V]/2 - F_v[0, S_cell_V]/2) + 4*D;
                sol_mat[E, 2] = 0;
                sol_mat[N, 2] = F_v[0, N_cell_V]/2 - D;
            elif mug[a_y, a_x + 1] == 0: #S wall
                sol_mat[S, 2] = 0;
                sol_mat[W, 2] = - F_u[0, W_cell_U]/2 - D;
                sol_mat[P, 2] = (F_u[0, E_cell_U]/2 - F_u[0, W_cell_U]/2 + F_v[0, N_cell_V]/2) + 4*D;
                sol_mat[E, 2] = F_u[0, E_cell_U]/2 - D;
                sol_mat[N, 2] = F_v[0, N_cell_V]/2 - D;
            elif mug[a_y + 2, a_x + 1] == 0: #N wall
                sol_mat[S, 2] = - F_v[0, S_cell_V]/2 - D;
                sol_mat[W, 2] = - F_u[0, W_cell_U]/2 - D;
                sol_mat[P, 2] = (F_u[0, E_cell_U]/2 - F_u[0, W_cell_U]/2 - - F_v[0, S_cell_V]/2 ) + 4*D;
                sol_mat[E, 2] = F_u[0, E_cell_U]/2 - D;
                sol_mat[N, 2] = 0;
            else: 
            #interior
                Pe_u = abs((F_u[0, E_cell_U]+F_u[0, W_cell_U])/(2*D));
                Pe_v = abs((F_v[0, N_cell_V]+F_v[0, S_cell_V])/(2*D));
                alternative_U = (F_u[0, E_cell_U]+F_u[0, W_cell_U])/2;
                alternative_V = (F_v[0, N_cell_V]+F_v[0, S_cell_V])/2;
                if (Pe_u > 2 and Pe_v > 2): #full upwind
                    if ((F_u[0, E_cell_U] > 0 and F_v[0, N_cell_V] > 0) or (alternative_U > 0 and alternative_V > 0)) :
                        a_S = -F_v[0, S_cell_V] - D;
                        a_W = -F_u[0, W_cell_U] - D;
                        a_E = - D;
                        a_N = - D;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                    elif ((F_u[0, E_cell_U] > 0 and F_v[0, N_cell_V] < 0) or (alternative_U > 0 and alternative_V < 0)):
                        a_S = - D;
                        a_W = -F_u[0, W_cell_U] - D;
                        a_E = - D;
                        a_N = F_v[0, N_cell_V] - D ;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                    elif ((F_u[0, E_cell_U] < 0 and F_v[0, N_cell_V] > 0) or (alternative_U < 0 and alternative_V > 0)):
                        a_S = -F_v[0, S_cell_V] - D;
                        a_W = - D;
                        a_E = F_u[0, E_cell_U] - D;
                        a_N = - D ;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                    elif ((F_u[0, E_cell_U] < 0 and F_v[0, N_cell_V] < 0) or (alternative_U < 0 and alternative_V < 0)):
                        a_S = - D;
                        a_W = - D;
                        a_E = F_u[0, E_cell_U] - D;
                        a_N = F_v[0, N_cell_V] - D ;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                elif Pe_u > 2 and Pe_v < 2: #upwind only in u-direction
                    if F_u[0, E_cell_U] > 0 or alternative_U > 0:
                        a_S = -F_v[0, S_cell_V]/2 - D;
                        a_W = -F_u[0, W_cell_U] - D;
                        a_E = - D;
                        a_N = F_v[0, N_cell_V]/2 - D;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                    elif F_u[0, E_cell_U] < 0 or alternative_U < 0:
                        a_S = -F_v[0, S_cell_V]/2 - D;
                        a_W = - D;
                        a_E = F_u[0, E_cell_U] - D;
                        a_N = F_v[0, N_cell_V]/2 - D;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                elif (Pe_u < 2 and Pe_v > 2): #upwind only in v-direction
                    if F_v[0, N_cell_V] > 0 or alternative_V > 0:
                        a_S = - F_v[0, S_cell_V] - D;
                        a_W = -F_u[0, W_cell_U]/2 - D;
                        a_E = F_u[0, E_cell_U]/2 - D;
                        a_N = - D;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                    elif F_v[0, N_cell_V] < 0 or alternative_V < 0:
                        a_S = - D;
                        a_W = -F_u[0, W_cell_U]/2 - D;
                        a_E =  F_u[0, E_cell_U]/2 - D;
                        a_N = F_v[0, N_cell_V] - D;
                        a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
                 
                        sol_mat[S, 2] = a_S;
                        sol_mat[W, 2] = a_W;
                        sol_mat[P, 2] = a_P;
                        sol_mat[E, 2] = a_E;
                        sol_mat[N, 2] = a_N;
                elif (Pe_u < 2 and Pe_v < 2): # full CD
                    a_S = -F_v[0, S_cell_V]/2 - D;
                    a_W = -F_u[0, W_cell_U]/2 - D;
                    a_E = F_u[0, E_cell_U]/2 - D;
                    a_N = F_v[0, N_cell_V]/2 - D;
                    a_P = (F_u[0, E_cell_U] - F_u[0, W_cell_U] + F_v[0, N_cell_V] - F_v[0, S_cell_V] - a_E - a_N - a_W - a_S );
        
                    sol_mat[S, 2] = a_S;
                    sol_mat[W, 2] = a_W;
                    sol_mat[P, 2] = a_P;
                    sol_mat[E, 2] = a_E;
                    sol_mat[N, 2] = a_N;