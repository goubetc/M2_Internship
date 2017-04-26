function times = plotTimeIter(target, nbProj, timeTable, numProj, nBreg)
    times = timeTable(target,nbProj, :);
    h = figure;
    plot(nBreg(:),times(:)); axis tight; title(['Execution time target ' num2str(target) ' nbProj ' num2str(numProj(nbProj)) ]);
            colormap gray; drawnow;
end