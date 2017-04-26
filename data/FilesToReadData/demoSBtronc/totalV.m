function res = totalV(targetImg, imgRec)
    sumDiff = targetImg(:,:) - imgRec(:,:);
    [gDiff_x, gDiff_y] = gradient(sumDiff);
    res = sum(sum(sqrt(gDiff_x.^2 + gDiff_y.^2)));