
% =====================================================================
% Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        
        solMat  = reshape(sol,N);%%% [2087, 2087]
        
        % Laplacian part
%         u           = zeros(2087,2087)%%%
%         lambda      = 1;%%%
%         
%         A           = @(x)(radon(x,thetas));%%%
%         AT          = @(x)(iradon(x,thetas,'linear','None',1.0,output_size));%%%

%         A           = @(x)(A(x)*scaleA);%%%
%         AT          = @(x)(AT(x)*scaleAT);%%%
        
        
        [rows,cols,height] = size(u);
        Dx = [rows,cols,height];
        bTV     = lambda*(Dxt(Dx(solMat))+Dyt(Dy(solMat)));
        
        % Jacobian part
        bJac    = mu*AT(A(solMat));
        
        % Stability term
        bG      = gamma*sol;
        
        b       = bTV(:) + bJac(:) + bG(:);
    end
% =====================================================================
