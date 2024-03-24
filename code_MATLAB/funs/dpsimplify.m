function [ps,ix] = dpsimplify(p,tol)
% Recursive Douglas-Peucker Polyline Simplification, Simplify
%
% [ps,ix] = dpsimplify(p,tol)
%
% dpsimplify uses the recursive Douglas-Peucker line simplification 
% algorithm to reduce the number of vertices in a piecewise linear curve 
% according to a specified tolerance. The algorithm is also know as
% Iterative Endpoint Fit. It works also for polylines and polygons
% in higher dimensions.
%
% In case of nans (missing vertex coordinates) dpsimplify assumes that 
% nans separate polylines. As such, dpsimplify treats each line
% separately.
%
% For additional information on the algorithm follow this link
% http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
%
% Input arguments
%
%     p     polyline n*d matrix with n vertices in d 
%           dimensions.
%     tol   tolerance (maximal euclidean distance allowed 
%           between the new line and a vertex)
%
% Output arguments
%
%     ps    simplified line
%     ix    linear index of the vertices retained in p (ps = p(ix))
%
% Examples
%
% 1. Simplify line 
%
%     tol    = 1;
%     x      = 1:0.1:8*pi;
%     y      = sin(x) + randn(size(x))*0.1;
%     p      = [x' y'];
%     ps     = dpsimplify(p,tol);
%
%     plot(p(:,1),p(:,2),'k')
%     hold on
%     plot(ps(:,1),ps(:,2),'r','LineWidth',2);
%     legend('original polyline','simplified')
%
% 2. Reduce polyline so that only knickpoints remain by 
%    choosing a very low tolerance
%
%     p = [(1:10)' [1 2 3 2 4 6 7 8 5 2]'];
%     p2 = dpsimplify(p,eps);
%     plot(p(:,1),p(:,2),'k+--')
%     hold on
%     plot(p2(:,1),p2(:,2),'ro','MarkerSize',10);
%     legend('original line','knickpoints')
%
% 3. Simplify a 3d-curve
% 
%     x = sin(1:0.01:20)'; 
%     y = cos(1:0.01:20)'; 
%     z = x.*y.*(1:0.01:20)';
%     ps = dpsimplify([x y z],0.1);
%     plot3(x,y,z);
%     hold on
%     plot3(ps(:,1),ps(:,2),ps(:,3),'k*-');
%
%
%
% Author: Wolfgang Schwanghart, 13. July, 2010.
% w.schwanghart[at]unibas.ch
if nargin == 0
    help dpsimplify
    return
end
error(nargchk(2, 2, nargin))
% error checking
if ~isscalar(tol) || tol<0
    error('tol must be a positive scalar')
end
% if input track is table, convert
if istable(p)
    p = p{:,{'x' 'y'}};
end
% nr of dimensions
dims    = size(p,2);

% start the recursive algorithm
ixe     = size(p,1);
ixs     = 1;
% logical vector for the vertices to be retained
I   = true(ixe,1);
% call recursive function
p   = simplifyrec(p,tol,ixs,ixe);
ps  = p(I,:);
% if desired return the index of retained vertices
if nargout == 2
    ix  = find(I);
end

% If input was table, give table as output
if istable(p)
    ps = array2table(ps,'VariableNames',{'x','y'});
end
% _________________________________________________________
function p  = simplifyrec(p,tol,ixs,ixe)
 
    % calculate shortest distance of all points to the line from ixs to ixe
    % subtract starting point from other locations
    pt = bsxfun(@minus,p(ixs+1:ixe,:),p(ixs,:));
    % end point
    a = pt(end,:)';
    beta = (a' * pt')./(a'*a);
    b    = pt-bsxfun(@times,beta,a)';
    d    = hypot(b(:,1),b(:,2));
    
    % identify maximum distance and get the linear index of its location
    [dmax,ixc] = max(d);
    ixc  = ixs + ixc; 
    
    % if the maximum distance is smaller than the tolerance remove vertices
    % between ixs and ixe
    if dmax <= tol
        if ixs ~= ixe-1
            I(ixs+1:ixe-1) = false;
        end
    % if not, call simplifyrec for the segments between ixs and ixc (ixc
    % and ixe)
    else   
        p   = simplifyrec(p,tol,ixs,ixc);
        p   = simplifyrec(p,tol,ixc,ixe);
    end
end
end
