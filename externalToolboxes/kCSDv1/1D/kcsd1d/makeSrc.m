function src = makeSrc(X, ext, nSrc)
    xMin = min(X) - ext;
    xMax = max(X) + ext;
    dx = (xMax - xMin)/nSrc;
    src = min(X):dx:max(X);
end