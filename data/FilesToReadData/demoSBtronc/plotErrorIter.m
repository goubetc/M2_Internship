function  errors = plotErrorIter(target, proj, err, numProj, nBreg)
    lengthErros = length(nBreg);
    errors = zeros(lengthErros,1);
    for i = 1:lengthErros
        errImg = err{target, proj, i}{1, 1};
        errors(i) = errImg(length(errImg));
    end
    h = figure;
    plot(nBreg(:),errors(:)); axis tight; title(['Error image ' num2str(target) ' nbProj ' num2str(numProj(proj)) ]);
            colormap gray; drawnow;
end
    
% figure;
% plot(1:5000, errors{1,1}); axis tight; 
%             colormap gray; drawnow;