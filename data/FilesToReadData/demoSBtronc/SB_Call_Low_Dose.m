function [u,errAll, ctrst, timeExec] = SB_Call_Low_Dose(target, numProj, nBreg, noise)

sizeImg= size(target);
output_size = sizeImg(1);

thetas      = (1:numProj)*180/numProj;
p           = radon(target, thetas);
N           = size(p);

if nargin == 4
    %Add poisson noise
    epsilon = 5; % corrects for negative vlaues in the pre logged sinogram data that would otherwise results in infinite projection values
    Noise = noise*exp(-p);
    N_noise = poissrnd(Noise);
    p = -log(N_noise/noise);
    idx = isinf(p);
    p(idx) = -log(epsilon/noise);
%     figure; imagesc(p); colormap gray;colormap(flipud(colormap)); colorbar; 
%     caxis('auto'); 
%     title('Projections with poisson noise '); drawnow;
end

% Define Fwd and adjoint operators
A           = @(x)(radon(x,thetas, N(1)));
AT          = @(x)(iradon(x,thetas,'linear','None',1.0, output_size));%%%%,output_size CHANGED

uretro      = iradon(p,thetas);
%figure; imagesc(uretro); colormap gray;colormap(flipud(colormap)); colorbar; 

flagNorm    = 0;

mu          = 1;
lambda      = 1;
gamma       = 1;
alpha       = 0.1;
nInner      = 1;
%nBreg       = 1000;
%arr_idx     = ones(N);


before1=clock;
[u,errAll, ctrst]  = TV_SB_2D_Gen(A,AT,p,[output_size output_size],mu, lambda, gamma, alpha, nInner, nBreg,target);
timeExec= etime(clock,before1);
