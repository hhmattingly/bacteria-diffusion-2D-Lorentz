function [past_squares_pos,new_square_pos,circ_pos] = buildSquares(xi,past_squares_pos,circ_pos,square_halfL,rho)
% [past_squares_pos,new_square_pos,circ_pos] = buildSquares(xi,past_squares_pos,circ_pos,square_halfL,rho)
%
% Add a new square to the simulation domain.
%
% Inputs:
% xi: Current cell position. 1 by d vector, where d is the number of spatial dimensiones
% past_squares_pos: Positions of the centers of previously-placed squares. d by Nsquare matrix.
% circ_pos: Positions of the centers of previously-placed circular obstacles. d by Ncircles matrix.
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar.
% rho: Dimensionless number density of obstacles (per area). scalar.
%
% Outputs:
% past_squares_pos: As above, but with an additional column for the new square.
% new_square_pos: Position of the new square. d by 1 vector.
% circ_pos: As above, but with additional columns for new circular obstacles whose centers lie inside the new square (and not in any previously placed squares).
%
% Henry H. Mattingly, November 2023

% set down a cube centered the cell's location
new_square_pos = xi';

% add spheres
new_circ_pos = placeSpheres(square_halfL,rho,new_square_pos,past_squares_pos); % don't need to shift by the latest cube_pos -- already done
circ_pos = [circ_pos, new_circ_pos];

past_squares_pos = [past_squares_pos, new_square_pos];

end