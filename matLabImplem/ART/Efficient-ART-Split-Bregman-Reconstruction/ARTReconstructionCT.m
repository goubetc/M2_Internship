function u   = ARTReconstructionCT(thetas,d,relaxParam,numIterART,u)
% u   = ARTReconstruction(A,d,relaxParam,numIterART,u)
%
% Algebraic reconstruction technique, also known as Kaczmarz method
%
% Inputs:
% 
% thetas        = angles of reconstruction
%               = d, data, nr x 1
% relaxParam    = relaxation parameter (between 0 and 1). Choose 0.1 to
% remove noise and 0.9 to get a closer fit of the data
% numIterART    = number of iterations
% u             = initial guess of the solution image, nc x 1. Initialize
% as zeros(nc,1) 
%
% Outputs:
%
% u             = solution image, nc x 1
%
%
% Code downloaded from the repository
% https://github.com/HGGM-LIM/Efficient-ART-Split-Bregman-Reconstruction
%
% If you use this code, please cite Chamorro-Servent et al. Use of Split
% Bregman denoising for iterative reconstruction in fluorescence diffuse
% optical tomography. J Biomed Opt, 18(7):076016, 2013.
% http://dx.doi.org/10.1117/1.JBO.18.7.076016
%
% Juan FPJ Abascal, Judit Chamorro-Servent, Juan Aguirre
% Departamento de Bioingenieria e Ingenieria Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, juchamser@gmail.com, desco@hggm.es

% Norm of rows
%Anorm               = sum(A.*A,2);
img = phantom(64);
bw = double(im2bw(img));

% Loop across iterations
for it = 1:numIterART
    % Loop across angles
    for ir = 1:size(thetas)
        % Loop across u
        for i = 1:64
            %for j = 1:64
            %bw = zeros(64,64);
            %bw(i,:)=1;
            Ai = radon(bw, thetas(ir));
            d2 = radon(u(:,:), thetas(ir));
            %Air = radon(bw(:,:), thetas(ir));
            num         = d(ir)-d2;
            %u(i,j)           = u(i,j) + mtimes((relaxParam*num/((Ai') * Ai)),Ai);
            u           = u + mtimes((.9*num/((bw(i,:)') * bw(i,:))),bw(i,:));
        end %i
        
    end % ir
end % it