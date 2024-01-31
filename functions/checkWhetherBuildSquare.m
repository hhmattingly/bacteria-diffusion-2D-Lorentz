function buildsquare = checkWhetherBuildSquare(xi,past_squares_pos,square_halfL,beta,dt,r)
% buildsquare = checkWhetherBuildSquare(xi,past_squares_pos,square_halfL,beta,dt,r)
%
% Check whether a new square needs to be added to the simulation domain,
% i.e. if the cell approaches the edge of an existing square.
%
% Inputs:
% xi: Current cell position. 1 by d vector, where d is the number of spatial dimensiones
% past_squares_pos: Positions of the centers of previously-placed squares. d by Nsquare matrix.
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar.
% beta: Dimensionless cell swimming speed.
% dt: Time step.
% r: Cell radius.
%
% Outputs:
%
% buildsquare: Whether or not to build a new square. logical scalar.
%
% Henry H. Mattingly, November 2023

% check whether each cell is close to any cube edge
Nsquares = size(past_squares_pos,2);
inside_inds = find(all(abs(repmat(xi',1,Nsquares) - past_squares_pos)<square_halfL));

% refine to just the closest center
[~,mni] = min(sqrt(sum((repmat(xi',1,length(inside_inds))-past_squares_pos(:,inside_inds)).^2,1)));

% now find the closest wall
center_dists = abs(xi'-past_squares_pos(:,inside_inds(mni)));

% output
buildsquare = any(center_dists > (square_halfL-(1+r)-2*beta*dt));


end