function [s, err_norm, iter_time]=greed_omp_chol_SparseLabWrap(A,Pt,x,m,s_initial,STOPCRIT,STOPTOL,MAXITER,verbose,comp_err,comp_time)
% Wrapper function for SparseLab SolveOMP algorithm
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Make P and Pt functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(A,'float')
    P =@(z) A*z;
    Pt =@(z) A'*z;
elseif isa(A,'object')
    P =@(z) A*z;
    Pt =@(z) A'*z;
elseif isa(A,'function_handle')
    if exist(Pt,'function_handle')
        P=A;
    else
       display('If P is a function handle, please also specify Pt. Exiting.') 
       return  
    end
else
    display('P is of unsupported type. Use matrix, function handle or object. Exiting.') 
    return  
end


    display('Using SolveOMP from SparseLab. Many options currently not available for this solver.')
    
    MyMap=@(mode, m, n, x, I, dim) Map(P,Pt,mode, m, n, x, I, dim);
    
    
    lambdaStop=0;
    
    if strcmp(STOPCRIT,'M')
        MAXITER=STOPTOL;
    elseif strcmp(STOPCRIT,'M')
        lambdaStop=STOPTOL;
    else
        error('Not yet available stoping criterion.')
    end
    
    tic
    [s] = SolveOMP(MyMap, x, N, MAXITER, lambdaStop, 0, verbose);
    TIME=toc;
    if comp_err
        display('Error computation is currently not available. Returning empty array.')
    end
    if comp_time
        display('Full time computation is currently not available. Returning max time.')
        iter_time=TIME;
    else
        iter_time=[];
    end
    err_norm=[];

end

function y=Map(P,Pt,mode, m, n, x, I, dim)

    if mode == 1
        y=P(x);
    else
        y=Pt(x);
    end
    
end