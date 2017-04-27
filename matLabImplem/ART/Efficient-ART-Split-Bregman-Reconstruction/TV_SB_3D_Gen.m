function [u,err, ctrst] = TV_SB_3D_Gen(A,AT,ANorm,ATNorm,f, N,mu, lambda, gamma, alpha, nInner, nBreg,varargin)
% [u] = TV_SB_3D_Gen(A,AT,f, N,mu, lambda, gamma, alpha, nInner, nBreg)
% [u,err] = TV_SB_3D_Gen(A,AT,f, N,mu, lambda, gamma, alpha, nInner, nBreg,uTarget)
% [u,err] = TV_SB_3D_Gen(A,AT,f, N,mu, lambda, gamma, alpha, nInner, nBreg,uTarget,optSolver)
% [u,err] = TV_SB_3D_Gen(A,AT,f, N,mu, lambda, gamma, alpha, nInner, nBreg,uTarget,optSolver,flagNorm)

% Krylov convergence criterion: decrease to improve precision for solving
% the linear system, increase to go faster (1e-2 to 1e-4)
tolKrylov   = 1e-2;

flagNorm    = 1;
flagDisplay = 0;

% Normalize data
% normFactor  = getNormalizationFactor(f,N(2));
normFactor  = getNormalizationFactor(f(1,:,:),N(2));

%  normFactor  = 1
f           = normFactor*f;

% Stopiing criteria
errMin      = 5e-4; 
errPrevMin  = 1e-4;

switch nargin
    case 13
        uTarget     = varargin{1};
        uTarget     = normFactor*uTarget;
    case 14
        uTarget     = varargin{1};
        uTarget     = normFactor*uTarget;
        optSolver   = varargin{1};
        tolKrylov   = optSolver.tolKrylov;
    case 15
        uTarget     = varargin{1};
        uTarget     = normFactor*uTarget;
        optSolver   = varargin{2};
        if isempty(optSolver) == 0
            tolKrylov   = optSolver.tolKrylov;
        end
        flagNorm    = varargin{3};
    case 16
        uTarget     = varargin{1};
        uTarget     = normFactor*uTarget;
        optSolver   = varargin{2};
        if isempty(optSolver) == 0
            tolKrylov   = optSolver.tolKrylov;
        end
        flagNorm    = varargin{3};
        stopCrit    = varargin{4};
        errMin      = stopCrit.errMin; 
        errPrevMin  = stopCrit.errPrevMin;
end % nargin

errAll      = zeros(nBreg,1);
obj0        = zeros(nBreg,1);
rate        = zeros(nBreg,1);
errPrev     = zeros(nBreg,1);
nnzx        = zeros(nBreg,1);


% Reserve memory for the auxillary variables
rows        = N(2);
cols        = N(3);
profs       = N(1);
f0          = f;
u           = zeros(profs,rows,cols);
uPrev       = u;
x           = zeros(profs,rows,cols);
y           = zeros(profs,rows,cols);
z           = zeros(profs,rows,cols);

bx          = zeros(profs,rows,cols);
by          = zeros(profs,rows,cols);
bz          = zeros(profs,rows,cols);

ctrst       = zeros(ceil(nBreg/100),rows);

% Normalize the forward and adjoint operators to comply with adjoint
% operator definition and |J'Jx|~|x| and to have JJ' with diagonal of the
% order of 1 (so all Hessian terms have values on the same order of
% magnitude)
if flagNorm == 1
    tmp1        = ones(N(2:3));
    Atmp1       = ANorm(tmp1);
    tmp2        = reshape(1:prod(N(2:3)),N(2:3));
    Atmp2       = ANorm(tmp2);
    AAtmp2      = ATNorm(Atmp2);
    a           = abs(tmp1(:)'*AAtmp2(:)/(Atmp1(:)'*Atmp2(:)));
    b           = norm(tmp2(:))/norm(AAtmp2(:));
    scaleA      = sqrt(a*b);
    scaleAT     = sqrt(b/a);
    
    A           = @(x)(A(x)*scaleA);
    AT          = @(x)(AT(x)*scaleAT);
    
    % % Test adjoint definition
    % tmp1        = ones(N);
    % Atmp1       = J(tmp1);
    % tmp2        = reshape(1:prod(N),N);
    % Atmp2       = J(tmp2);
    % AAtmp2      = JT(Atmp2);
    % abs(tmp1(:)'*AAtmp2(:))/(Atmp1(:)'*Atmp2(:))
    % norm(tmp2(:))/norm(AAtmp2(:))
    %for i = 1:N(1)
      uTarget        = uTarget./scaleA;
    %end
end

murf        = mu*AT(f);
outer       = 0;
flagIt      = 1;
uBest       = u;
errBest     = inf;
uDP         = u;
errDPBest   = inf;
iterDP      = 1;
it          = 1;
%  Do the reconstruction
for outer = 1:nBreg
    for inner = 1:nInner;
        % update u
        rhs         = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by)+lambda*Dzt(z-bz)+gamma*u;
        
        u           = reshape(krylov3D(rhs(:)),N);
        
        dx          = Dx(u);
        dy          = Dy(u);
        dz          = Dz(u);
        
        % update x and y
        [x,y,z]       = shrink3(dx+bx,dy+by,dz+bz,alpha/lambda);
        
        % update bregman parameters
        bx          = bx+dx-x;
        by          = by+dy-y;
        bz          = bz+dz-z;
    end   % inner loop
    
    fForw           = A(u);
    f               = f + f0-fForw;
    murf            = mu*AT(f);
    
    nnzx(it)        = nnz(x(:));
    if nargin >= 11
%         if and(outer == 1, inner ==1)
%             it = it + 1;
%         end
        it          = outer;
        
        % Solution error norm, data fidelity and rate of convergence
        errAll(it)    = norm(uTarget(:)-u(:))/norm(uTarget(:));
        %obj0(it)      = 0.5*norm(reshape(A(u)-f0,[],1));
        obj0(it)      = norm(reshape((A(u)-f0)/normFactor,[],1));        
        rate(it)      = norm(u(:)-uTarget(:))/norm(uTarget(:)-uPrev(:));
        errPrev(it)   = norm(u(:)-uPrev(:))/norm(uPrev(:)+1e-2);
        
        if errAll(it) <= errBest
            errBest     = errAll(it);
            uBest       = u;
        end
        
        if flagDisplay == 1 && any([it ==1, it == 10, rem(it, 50) == 0])
            close;
            h=figure;
            subplot(3,2,1);
            imagesc((u)); title(['u, iter. ' num2str(it)]); colorbar;
            subplot(3,2,2);
            imagesc((x));  colorbar; title('x');
            subplot(3,2,3);
            plot(errAll(1:it)); axis tight; title(['Sol. error' ]);
            colormap gray;
            subplot(3,2,4);
            plot(obj0(1:it)); axis tight; title(['obj0' ]);
            subplot(3,2,5);
            plot(rate(1:it)); axis tight; title(['Rate' ]);
            subplot(3,2,6);
            plot(errPrev(1:it)); axis tight; title(['Err prev' ]);
            axis([0 it 0 1]);
            colormap(flipud(colormap));
            drawnow;
        end % rem
    end % nargin
    
    % get a line for contrast
     if mod(outer, 100) == 1
        
        ctrst(outer, :) = u(2,ceil(rows/2),:);
    end
    
    uPrev       = u;
end % outer

if flagDisplay == 1 && nargin >= 13
    close;
    h=figure;
    subplot(3,2,1);
    imagesc(abs(u)); title(['u, iter. ' num2str(it)]); colorbar;
    subplot(3,2,2);
    imagesc(abs(x));  colorbar; title('x');
    subplot(3,2,3);
    plot(errAll(1:it)); axis tight; title(['Sol. error' ]);
    colormap gray;
    subplot(3,2,4);
    plot(obj0(1:it)); axis tight; title(['obj0' ]);
    subplot(3,2,5);
    plot(rate(1:it)); axis tight; title(['Rate' ]);
    subplot(3,2,6);
    plot(errPrev(1:it)); axis tight; title(['Err prev' ]);
    axis([0 it 0 1]);
    drawnow;    
end

if  nargin >= 13
    err.errAll  = errAll;
    err.errPrev = errPrev;
    err.rate    = rate;
    err.obj0    = obj0;
    err.nnzx    = nnzx;
end
err.x       = x;
err.y       = y;



% undo the normalization so that results are scaled properly
if nargin >= 13
    u = u/normFactor;
else
    u = u/normFactor;
end

if flagNorm == 1
    u = u*scaleA;  
end
%u = u*(scaleA)/(normFactor); % 

    function normFactor = getNormalizationFactor(f,n)
        
        %normFactor = 1/norm(f(:)/size(R==1,1));
        normFactor = 1/norm(f(:)/n);
    end


    function d = Dx(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
        d(:,1,:) = u(:,1,:)-u(:,cols,:);
    end

    function d = Dxt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,1:cols-1,:) = u(:,1:cols-1,:)-u(:,2:cols,:);
        d(:,cols,:) = u(:,cols,:)-u(:,1,:);
    end

    function d = Dz(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
        d(1,:,:) = u(1,:,:)-u(rows,:,:);
    end

    function d = Dzt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(1:rows-1,:,:) = u(1:rows-1,:,:)-u(2:rows,:,:);
        d(rows,:,:) = u(rows,:,:)-u(1,:,:);
    end

    function d = Dy(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,:,2:height) = u(:,:,2:height)-u(:,:,1:height-1);
        d(:,:,1) = u(:,:,1)-u(:,:,height);
    end

    function d = Dyt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,:,1:height-1) = u(:,:,1:height-1)-u(:,:,2:height);
        d(:,:,height) = u(:,:,height)-u(:,:,1);
    end

%     function [xs,ys] = shrink2(x,y,lambda)
%         
%         s = sqrt(x.*conj(x)+y.*conj(y));
%         ss = s-lambda;
%         ss = ss.*(ss>0);
%         
%         s = s+(s<lambda);
%         ss = ss./s;
%         
%         xs = ss.*x;
%         ys = ss.*y;
%         
%     end

function [xs,ys,zs] = shrink3(x,y,z,lambda)
        
        s = sqrt(x.*conj(x)+y.*conj(y)+z.*conj(z));
        ss = s-lambda;
        ss = ss.*(ss>0);
        
        s = s+(s<lambda);
        ss = ss./s;
        
        xs = ss.*x;
        ys = ss.*y;
        zs = ss.*z;
        
    end
% % =====================================================================
% % Krylov solver subroutine
% % X = GMRES(A,B,RESTART,TOL,MAXIT,M)
% % bicgstab(A,b,tol,maxit)
%     function dx = krylov(r)
%         %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
%         [dx,flag,relres,iter] = bicgstab(@jtjx, r, tolKrylov, 100);
%     end
% 
% % =====================================================================
% % Callback function for matrix-vector product (called by krylov)
%     function b = jtjx(sol)
%         solMat  = reshape(sol,N);
%         
%         % Laplacian part
%         bTV     = lambda*(Dxt(Dx(solMat))+Dyt(Dy(solMat)));
%         
%         % Jacobian part
%         bJac    = mu*AT(A(solMat));
%         
%         % Stability term
%         bG      = gamma*sol;
%         
%         b       = bTV(:) + bJac(:) + bG(:);
%     end


% =====================================================================
% Krylov solver subroutine
% X = GMRES(A,B,RESTART,TOL,MAXIT,M)
% bicgstab(A,b,tol,maxit)
    function dx = krylov3D(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        [dx,flag,relres,iter] = bicgstab(@jtjx3D, r, tolKrylov, 100);
    end

% =====================================================================
% Callback function for matrix-vector product (called by krylov)
    function b = jtjx3D(sol)
        solMat  = reshape(sol,N);
        
        % Laplacian part
        bTV     = lambda*(Dxt(Dx(solMat))+Dyt(Dy(solMat))+Dzt(Dz(solMat)));
        
        % Jacobian part
        bJac    = mu*AT(A(solMat));
        
        % Stability term
        bG      = gamma*sol;
        
        b       = bTV(:) + bJac(:) + bG(:);
    end

% =====================================================================
end

