% =====================================================================
% Krylov solver subroutine
% X = GMRES(A,B,RESTART,TOL,MAXIT,M)
% bicgstab(A,b,tol,maxit)
    function dx = krylov(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        tolKrylov   = 1e-2;
        [dx,flag,relres,iter] = bicgstab(@jtjx, r, tolKrylov, 100);
    end
