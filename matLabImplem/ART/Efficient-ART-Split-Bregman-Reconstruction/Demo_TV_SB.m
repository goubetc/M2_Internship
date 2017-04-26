%Shepp-Logan phantom
N   = 64;
uin = phantom('Modified Shepp-Logan', N);
figure, imagesc(uin), colorbar, colormap gray

Nthetas=100;
thetas = 1:floor(360/Nthetas):360;

numIterART = 200;
relaxParam = 0.5;

%experimental projections
d           = radon(uin, thetas);
indZeros    = randperm(Nthetas);
d(:,indZeros(1:50)) = 0;
figure, imagesc(d), colorbar, colormap gray

NData   = prod(size(d));

uRetro      = iradon(d, thetas,'linear','none',N);
figure, imagesc(uRetro), colorbar, colormap gray

uFBP        = iradon(d, thetas,'spline','Ram-Lak',1,N);
figure, imagesc(uFBP), colorbar, colormap gray

% u   = zeros(N);
u   = uFBP;
% Ru  = radon(u, thetas);
% figure, imagesc(abs(d-Ru)), colorbar, colormap gray

% Define Fwd and adjoint operators
A           = @(x)(radon(x,thetas));
AT          = @(x)(iradon(x,thetas,'linear','None',1.0,N(1)));
N = [N N];
% operator definition and |J'Jx|~|x|
flagNorm    = 1;
if flagNorm == 1
    tmp1        = ones(N);
    Atmp1       = A(tmp1);
    tmp2        = reshape(1:prod(N),N);
    Atmp2       = A(tmp2);
    AAtmp2      = AT(Atmp2);
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
    
    uTarget         = u/scaleA;
end

% SD
% spectralNorm= norm(A(ones(N)))/norm(ones(N));
% % delta       = 0.2*2/(spectralNorm^2);
% delta       = 1e-4
% u           = zeros(N);
% dk          = d;
% numIter     = 20;
% % mu          = 1/delta;
% mu          = 1;
% alpha       = 1;
% lambda      = 1;
% nin         = 2;
% nout        = 5;
% h           = figure;
% err         = inf;
% for it = 1:numIter
% %    du       = -AT(A(u)-dk);
%     v        = u - delta*AT(A(u)-dk);  
%    %objLS(x,d,A)
% %    delta    = 1e-4;
% %    [delta, err] = toastLineSearch (u, du, delta, err, @objLS,d,A);
% %    v        = u + delta*AT(A(u)-dk); 
%    u        = TV_SB_denoising_2D(v,mu,lambda,alpha,nin,nout);
% %    dk       = dk + d -A(u);
%    figure(h); 
%    subplot(2,1,1); imagesc(v); colorbar; axis image;    
%    subplot(2,1,2); imagesc(u); colorbar; axis image;
%    colormap gray; drawnow;
% end

% IMAGE RECONSTRUCTION using the Split Bregman method
% Recover the image
mu          = 1;
lambda      = 1;
gamma       = 1;
alpha       = 0.1;
nInner      = 1;
nBreg       = 200;
arr_idx     = ones(N);
[u,errAll]  = TV_SB_2D_Gen(A,AT,d,[N N],mu, lambda, gamma, alpha, nInner, nBreg,uin);
figure; 
imagesc(u); colorbar; axis image; colormap gray; 




