function [x,flag] = placeCellsInitially(Ncells,dx,circ_pos,circ_rad,r)
% [x,flag] = placeCellsInitially(Ncells,dx,circ_pos,circ_rad,cell_rad)
%
% Place cells in the void space of the simulation domain.
%
% Inputs:
% Ncells: How many cells to place in the given realization of the environment. scalar
% dx: Side length of the box (square or cube) in which to place cells. scalar
% circ_pos: Positions of the centers of previously-placed circular obstacles. d by Ncircles matrix.
% circ_rad: Circle radius. Should be 1. scalar
% r: Cell radius. scalar
%
% Outputs:
% x: Positions of cells. Ncells by d matrix
% flag: Whether or not there was a problem placing cells. This function
% terminates if placing cells takes longer than 1 second of wall time to
% complete. If it fails to place all cells, flag = 0. logical scalar.
%
% Henry H. Mattingly, November 2023

d = size(circ_pos,1);

circ_pos = circ_pos';
circ_pos = circ_pos(all(abs(circ_pos)<2*dx,2),:);
dx0 = dx;
dx = dx*ones(1,d);
x = [];
n=0;
flag = 0;

t1=tic;
ni = max(100, floor(Ncells/2));
while n<Ncells
    x_i = (rand(ni,d)-1/2).*repmat(dx,ni,1);
    Dnew = pdist2(x_i,circ_pos);
    Dnew = min(Dnew,[],2);
    x_i(Dnew<=(circ_rad+r),:) = [];

    if numel(x_i)==0 && dx(1)<=2/1.1*dx0
        dx = 1.1*dx; % in case the center is fully covered by chance
    end
    
    x = [x; x_i];

    n = size(x,1);

    if toc(t1)>1 % sec
        flag = 1;
        x = [];
        return
    end
end

if flag==0
    x = x(1:Ncells,:);
end

end

