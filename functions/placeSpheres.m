function new_circ_pos = placeSpheres(square_halfL,rho,new_square_pos,past_squares_pos)
% new_circ_pos = placeSpheres(square_halfL,rho,new_square_pos,past_squares_pos)
%
% Add new circular obstacles to the simulation domain.
%
% Inputs:
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar.
% rho: Dimensionless number density of obstacles (per area). scalar.
% new_square_pos: Position of the new square. d by 1 vector., where d is the number of spatial dimensiones
% past_squares_pos: Positions of the centers of previously-placed squares. d by Nsquare matrix.
%
% Outputs:
% new_circ_pos: Positions of the centers of newly-placed circular obstacles. These centers lie inside the square centered at new_square_pos, and not inside any previously placed squares. d by Ncircles matrix.
%
% Henry H. Mattingly, November 2023

[d, ~] = size(new_square_pos);
[~, Npast] = size(past_squares_pos);

% mean number of spheres in the new cube
Acube = (2*square_halfL)^d;
mean_Ncirc = rho*Acube;

% generate new spheres
Ncirc_new = poissrnd(mean_Ncirc);
new_circ_pos = 2*square_halfL*(rand(d,Ncirc_new)-1/2)+new_square_pos;

% check if any lie inside other squares
if ~isempty(past_squares_pos)
    rmInds = any(all(abs(repmat(new_circ_pos,[1,1,Npast])-repmat(permute(past_squares_pos,[1,3,2]),[1,Ncirc_new,1]))<square_halfL,1),3);
    new_circ_pos(:,rmInds) = [];
end

end