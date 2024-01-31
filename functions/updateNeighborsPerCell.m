function neighbors = updateNeighborsPerCell(cell_pos,circ_pos,circ_rad,cell_rad,max_travel_R)
% neighbors = updateNeighborsPerCell(cell_pos,circ_pos,circ_rad,cell_rad,max_travel_R)
%
% Find obstacles that are close to the cell's current position.
%
% Inputs:
% cell_pos: Cell position. 1 by d matrix, where d is the number of dimensions
% circ_pos: Positions of the centers of the obstacles. Ncircles by d matrix, where d is the number of dimensions
% circ_rad: Obstacle radius. Should always be 1. scalar
% max_travel_R: Maximum distance the cell travels before the next call of this function in the simulation.
%
% Outputs:
% neighbors: Vector of the indices of obstacles that are close to the cell's current position. The obstacle centers are located at circ_pos(neighbors,:). N by 1 vector.
%
% Henry H. Mattingly, November 2023

Ncirc = size(circ_pos,1);

neighbors = find(sqrt(sum((circ_pos-repmat(cell_pos,Ncirc,1)).^2,2))<(max_travel_R+2*(circ_rad+cell_rad))); 


end