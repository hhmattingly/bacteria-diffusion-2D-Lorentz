function [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos] = generateICs(cell_params,env_params,sim_params,genTheta)
% [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos] = generateICs(cell_params,env_params,sim_params,genTheta)
%
% Generate initial conditions and arrays for a simulation.
%
% Inputs:
% cell_params: A structure with parameters describing the swimmer. 
%   Fields: 
%   beta: normalized cell swimming speed. scalar
% env_params: A structure with parameters describing the environment.
%   Fields:
%   d: number of spatial dimensions. scalar
%   gamma: mean chord length. scalar
%   square_halfL: half the side length of each square, in units of the obstacle radius. scalar.
%   circAreaFrac: the theoretical area fraction occupied by circular obstacles
%   rho: dimensionless number density of obstacles (per area). scalar
% sim_params: A structure with parameters describing the simulation.
%   Fields:
%   T: Duration of the simulation. scalar
%   dt: Simuation time step size. scalar
%   nt: Number of time simulation steps. scalar.
% genTheta: a function handle that randomly samples new heading angles, theta
%
% Outputs:
%
% xt_c: Matrix of cell position over time. 1 by d by nt matrix, for d dimensions
% vt_c: Matrix of cell velocity over time. 1 by d by nt matrix, for d dimensions
% tht_c: Matrix of cell heading angle over time. 1 by nt vector.
% tumbles_c: Matrix of cell run/tumble state over time. 0 indicates run state, 1 indicates tumble state. 1 by nt logical vector.
% contacts_c: Matrix of cell contact state over time. 0 indicates swimming in bulk, 1 indicates contact with one obstacle, 2 indicates contact with two obstacles. 1 by nt vector.
% circ_pos: Matrix containing the positions of the centers of the circular obstacles. d by Nircles matrix, for d dimensions
% squares_pos: Matrix containing the positions of the centers of the squares defining the simulation domain. d by 1 matrix here, for d dimensions
%
% Henry H. Mattingly, November 2023

beta = cell_params.beta;

d = env_params.d;
rho = env_params.rho;
square_halfL = env_params.square_halfL;
gamma = env_params.gamma;

circ_rad = 1; % obstacle radius
r = 0; % cell radius

nt = sim_params.nt;

if gamma<inf

    trapped = 1;
    x0 = [];

    while trapped

        circ_pos = [];
        new_square_pos = zeros(2,1); 
        circ_pos = placeSpheres(square_halfL,rho,new_square_pos,[]); 

        squares_pos = new_square_pos;

        % positions
        disp('Placing cells...')
        % Place 1 cell. R = 1 after rescaling.

        [x0,placeCellFlag] = placeCellsInitially(1,square_halfL/3,circ_pos,1,r); 
        if placeCellFlag
            disp(['WARNING: No cell placed for cell ' num2str(cell_ind) '!!! Simulation not run. \n'])
        else
            disp('done.')
        end

        if ~isempty(x0)
            % check whether the cell is fully surrounded by circles??
            if gamma<3
                trapped = checkWhetherTrapped(x0,circ_pos,circ_rad,square_halfL,squares_pos,rho);
            else
                trapped = 0;
            end
        else
            trapped = 1;
        end

    end
else

    squares_pos = [];
    circ_pos = [];

    x0 = zeros(1,d);
end

%                 figure;hold on
%                 nsph = size(circ_pos,2);
%                 for i = 1:nsph
%                     viscircles(circ_pos',1);
%                 end

%% store ICs
tht_c = nan(1,nt);
th0 = genTheta(1);
thi = th0;
tht_c(1) = thi;

ui = [cos(thi),sin(thi)];

vi = beta*ui;
vt_c = nan(1,d,nt);
vt_c(1,:,1) = vi;

xt_c = nan(1,d,nt);
xt_c(1,:,1) = x0;
% xi = x0;

tumbler = 0;

tumbles_c = nan(1,nt); 
tumbles_c(1) = tumbler;

contacts_c = [];

if gamma<inf
    
    contacts_c = nan(1,nt);
end



end