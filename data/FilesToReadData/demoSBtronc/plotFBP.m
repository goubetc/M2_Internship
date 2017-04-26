function [error, pnoise, tv, u] =  plotFBP(target, nbproj, noise)
    theta       = (1:nbproj)*180/nbproj;
    p = radon(target, theta);
    figure; imagesc(p); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('Projections '); drawnow;
    
    if nargin == 3
        %Add poisson noise
        epsilon = 5; % corrects for negative vlaues in the pre logged sinogram data that would otherwise results in infinite projection values
        Noise = noise*exp(-p);
        N_noise = poissrnd(Noise);
        p = -log(N_noise/noise);
        idx = isinf(p);
        p(idx) = -log(epsilon/noise);
        figure; imagesc(p); colormap gray;colormap(flipud(colormap)); colorbar; 
        caxis('auto'); 
        title('Projections with poisson noise '); drawnow;
    end
    u = iradon(p, theta, size(target, 1));
    figure; imagesc(u); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('FBP Target image '); drawnow;
    error = norm(target(:) - u(:))/norm(target(:));
    pnoise = psnr(target, u);
    
    
    sumDiff = target(:,:) - u(:,:);
   
    [gDiff_x, gDiff_y] = gradient(sumDiff);
    
    tv = sum(sum(sqrt(sumDiff.^2 + sumDiff.^2)));
end