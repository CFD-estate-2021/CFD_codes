#CGS METHOD
#We will now implement a very efficient iterative solution method
#The conjugate gradient S method

#it can be applied to all positive definite systems (the case of physical sciences)
#it works also for non-symmetric matrices

#THE FIRST CODE IS FOR A GENERAL MATRIX (SMALL!!!)
#THE SECOND CODE FOR THE SPARSE MATRICES WE ARE CONSIDERING IN THE CFD SIMULATIONS

import numpy as np

def CGS_solver_1(matrix, source, guess):
    
    #initialize with 0th residual
    residual_0 = source - matrix.dot(guess);
    
    #initialize direction vectors
    dir_vec = residual_0;
    C_dir_vec = residual_0;
    
    #initialize L-2 global residual with residual 0
    R2 = np.sqrt(residual_0.dot(residual_0))
    
    #initialize residual for iteration
    residual = residual_0
    
    #initialize solution
    solution = guess
    #this will be solution at 0th iteration
    
    #set tolerance
    tolerance = 10**-5
    
    while R2 > tolerance:
        
        #alpha^n+1, coefficient needed in solution update
        alpha = (residual_0.dot(residual))/(residual_0.dot(matrix.dot(dir_vec)));
        
        #vector needed to update solution and account for a-symmetric nature of the solution, G^n+1
        G_vec = C_dir_vec - alpha*(matrix.dot(dir_vec));
        
        #update the solution
        solution = solution + alpha*(C_dir_vec + G_vec);
        
        #new residual and new global residual
        residual_n_1 = source - matrix.dot(solution);
        R2 = np.sqrt(residual_n_1.dot(residual_n_1));
        
        #beta^n+1, needed in conjugate direction, search direction vector update
        beta = (residual_0.dot(residual_n_1))/(residual_0.dot(residual));
        
        #update conjugate direction
        C_dir_vec = residual_n_1 + beta*G_vec;
        
        #update search direction
        dir_vec = C_dir_vec + beta*(G_vec + beta*dir_vec);
        
        residual = residual_n_1
        
        print(solution, R2)
        
    return solution

#we create the version of the above algorithm for the sparsely occupied matrix that appears in the CFD problems I condisered

#input vectors are all numpy arrays
#transpositions are needed because the input vectors are arrays as well (1, N), which means that the right dimensions
#must be respected for the np.dot operation to work
#must respect the dimensional rules of linear algebra
def CGS_solver(sol_mat, b_source, guess_0):
    a, b = np.shape(b_source); #size of the system, needed for stopping the algorithm
    columns = int(b);
    
    def matrix_product_5bands(matrix, vector): #specifically for the sparse sol_mat matrix
        rows, columns = np.shape(vector);
        product = np.zeros((1, int(columns))); #list of zeros
        
        for a in range(0, columns):
            t_1 = int(a*5);
            t_2 = int(t_1 +1);
            t_3 = int(t_2 +1);
            t_4 = int(t_3 +1);
            t_5 = int(t_4 +1);
            
            S = int(sol_mat[t_1, 1]);
            W = int(sol_mat[t_2, 1]);
            P = int(sol_mat[t_3, 1]);
            E = int(sol_mat[t_4, 1]);
            N = int(sol_mat[t_5, 1]);
            
            term = sol_mat[t_1, 2]*vector[0, S] + sol_mat[t_2, 2]*vector[0, W] + sol_mat[t_3, 2]*vector[0, P] + sol_mat[t_4, 2]*vector[0, E] + sol_mat[t_5, 2]*vector[0, N];
            
            product[0, a] = term
            
        return product
    
    #initialize with 0th residual
    residual_0 = b_source - matrix_product_5bands(sol_mat, guess_0);
    
    #initialize direction vectors
    dir_vec = residual_0;
    C_dir_vec = residual_0;
    
    #initialize L-2 global residual with residual 0
    trans_residual_0 = residual_0.T;
    R2 = np.sqrt(residual_0.dot(trans_residual_0))
    
    #initialize residual for iteration
    residual = residual_0
    
    #initialize solution
    solution = guess_0
    #this will be solution at 0th iteration
    
    #set tolerance
    tolerance = 10**-5 #ADJUSTED ACCORDING OT THE PROBLEM
    #reference to control the convergence
    ref = 200;
    #number of iterations
    iter_N = 0;
    
    while R2 > tolerance:
        trans_residual = residual.T;
        
        #alpha^n+1, coefficient needed in solution update
        alpha = (residual_0.dot(trans_residual))/(residual_0.dot(matrix_product_5bands(sol_mat, dir_vec).T));
        
        #vector needed to update solution and account for a-symmetric nature of the solution, G^n+1
        G_vec = C_dir_vec - alpha*(matrix_product_5bands(sol_mat, dir_vec));
        
        #update the solution
        solution = solution + alpha*(C_dir_vec + G_vec);
        
        #new residual and new global residual
        residual_n_1 = b_source - matrix_product_5bands(sol_mat, solution);
        trans_residual_n_1 = residual_n_1.T
        R2_temporary = np.sqrt(residual_n_1.dot(residual_n_1.T));
        R2 = R2_temporary[0, 0];
        
        #beta^n+1, needed in conjugate direction, search direction vector update
        beta = (residual_0.dot(trans_residual_n_1))/(residual_0.dot(trans_residual));
        
        #update conjugate direction
        C_dir_vec = residual_n_1 + beta*G_vec;
        
        #update search direction
        dir_vec = C_dir_vec + beta*(G_vec + beta*dir_vec);
        
        residual = residual_n_1
        
        print(solution[0, ref], R2, iter_N)
        
        if iter_N == columns: #to stop the iterations if the system does not converge to the small value of the treshold
            break
        
        iter_N = iter_N + 1
        
    return solution

#the attractive thing about this method is that according to theory the convergence should happen in a number
#of iterations strictly less than the number of equations in the linear system

#THE REASON FOR NO CONVERGENCE IN LARGE SYSTEMS ARE GENERALLY DUE TO NUMERICAL ROUND-OFF ERRORS AND DISCRETISATION/MESH
#ERRORS