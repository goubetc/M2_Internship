function normFactor = getNormalizationFactor(f,n)
        
        %normFactor = 1/norm(f(:)/size(R==1,1));
        normFactor = 1/norm(f(:)/n);
    end