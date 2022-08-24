#THE FOLLOWING CODE IS A GAUSS-SIDEDEL SOLUTION ALGORITHM

#Global residual accuracy requirement
treshold = 10**-3;
r = treshold + 1;

#iteration number
iter_N = 0;
iteration_numb = 5000; #maximum number of iterations

while r > treshold:
    
    residual = 0; #must reset the residual at each iteration
    normalisation_factor = 0;
    
    for i in range(0, int(sum(cells_per_row))):
        
        S = i*5;
        W = S + 1;
        P = W + 1;
        E = P + 1;
        N = E + 1;
        
        s_1 = int(sol_mat[S, 1]);
        t_1 = (1/sol_mat[P, 2])*(sol_mat[S, 2]*scalar_CV[0, s_1]);
        w_2 = int(sol_mat[W, 1]);
        t_2 = (1/sol_mat[P, 2])*(sol_mat[W, 2]*scalar_CV[0, w_2]);
        e_4 = int(sol_mat[E, 1]);
        t_4 = (1/sol_mat[P, 2])*(sol_mat[E, 2]*scalar_CV[0, e_4]);
        n_5 = int(sol_mat[N, 1]);
        t_5 = (1/sol_mat[P, 2])*(sol_mat[N, 2]*scalar_CV[0, n_5]);
        
        S_u = (1/sol_mat[P, 2])*b_source[0, i];
        
        scalar_CV[0, i] = S_u - t_1 - t_2 - t_4 - t_5;
   
    for n in range(0, int(sum(cells_per_row))):
        S = n*5;
        W = S + 1;
        P = W + 1;
        E = P + 1;
        N = E + 1;
        
        s_1_n = int(sol_mat[S, 1]);
        t_1_n = sol_mat[S, 2]*scalar_CV[0, s_1_n];
        w_2_n = int(sol_mat[W, 1]);
        t_2_n = sol_mat[W, 2]*scalar_CV[0, w_2_n];
        p_3_n = int(sol_mat[P, 1]) #same, it is the diagonal
        t_3_n = sol_mat[P, 2]*scalar_CV[0, p_3_n];
        e_4_n = int(sol_mat[E, 1]);
        t_4_n = sol_mat[E, 2]*scalar_CV[0, e_4_n];
        n_5_n = int(sol_mat[N, 1]);
        t_5_n = sol_mat[N, 2]*scalar_CV[0, n_5_n];
        
        S_u_n = b_source[0, n];
        
        residual = residual + abs(t_1_n + t_2_n + t_3_n + t_4_n + t_5_n - S_u_n);
        normalisation_factor = normalisation_factor + abs(t_3_n);
    
    r = residual/normalisation_factor;
    
    print(scalar_CV[0, ref], r, iter_N);
    
    if iter_N >= iteration_numb:
        break
        
    iter_N = iter_N + 1