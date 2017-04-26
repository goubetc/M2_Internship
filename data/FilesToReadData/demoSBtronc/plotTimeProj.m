function times = plotTimeProj(target, it, timeTable, numProj, nBreg)
    times = timeTable(target,:, it);
    h = figure;
    plot(numProj(:),times(:)); axis tight; title(['Execution time target ' num2str(target) ' nbIter ' num2str(nBreg(it)) ]);
            colormap gray; drawnow;
end