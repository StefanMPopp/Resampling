function s = sMaker(smu, sParam, l, sDist)
switch sDist
    case 1 % "norm"
        % Normal/Uniform
        sPre = normrnd(smu,sParam,l,1);
        sPre(sPre<2) = 1; % Truncating at step length 1
        if sPre(1)==1; sPre(1) = 2; end % 1st step would otherwise be start
        s = cumsum(ceil(sPre));
    case "exp" || 2
        % Exponential
        s = cumsum(ceil(exprnd(smu,l,1)));
    case "poiss" || 3
        % Poisson
        s = cumsum(ceil(poissrnd(smu,l,1)));
    case "levy" || 4
        % Bounded Power Law
        s = cumsum(ceil(levyrnd(l,2,sParam))); % 2 = #of dimensions
        % Add boundedness: if outside bounds, draw again
end