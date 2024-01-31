function [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos,trapped] = mainSimLoop(xt_c,vt_c,tht_c,tumbles_c,contacts_c,cell_params,env_params,sim_params,circ_pos,squares_pos,RNs_tum_c,genTheta)
% [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos,trapped] = mainSimLoop(xt_c,vt_c,tht_c,tumbles_c,contacts_c,sim_params,cell_params,env_params,circ_pos,squares_pos,RNs_tum_c,genTheta)
%
% Generate initial conditions and arrays for a simulation.
%
% Inputs:
% xt_c: Matrix of cell position over time. 1 by d by nt matrix, for d dimensions
% vt_c: Matrix of cell velocity over time. 1 by d by nt matrix, for d dimensions
% tht_c: Matrix of cell heading angle over time. 1 by nt vector.
% tumbles_c: Matrix of cell run/tumble state over time. 0 indicates run state, 1 indicates tumble state. 1 by nt logical vector.
% contacts_c: Matrix of cell contact state over time. 0 indicates swimming in bulk, 1 indicates contact with one obstacle, 2 indicates contact with two obstacles. 1 by nt vector.
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
% circ_pos: Matrix containing the positions of the centers of the circular obstacles. d by Nircles_initial matrix, for d dimensions
% squares_pos: Matrix containing the positions of the centers of the squares defining the simulation domain. d by 1 matrix, for d dimensions
% RNs_tum_c: Matrix of random numbers sampled from a uniform distribution between 0 and 1 used for generating tumbles. 1 by nt matrix.
% genTheta: a function handle that randomly samples new heading angles, theta
%
% Outputs:
% xt_c: Same as the input, but populated with results of the simulation.
% vt_c: Same as the input, but populated with results of the simulation.
% tht_c: Same as the input, but populated with results of the simulation.
% tumbles_c: Same as the input, but populated with results of the simulation.
% contacts_c: Same as the input, but populated with results of the simulation.
% circ_pos: Same as the input, but populated with results of the simulation. The size of this matrix changes as the simulation domain grows. Final size is 3 by Nircles_final matrix.
% squares_pos: Same as the input, but populated with results of the simulation. The size of this matrix changes as the simulation domain grows. Final size is 3 by Nsquares_final matrix.
% trapped: Indicates whether the cell is fully enclosed by obstacles at any point in the simulation. logical scalar.
%
% Henry H. Mattingly, November 2023

beta = cell_params.beta;

gamma = env_params.gamma;
rho = env_params.rho;
square_halfL = env_params.square_halfL;

dt = sim_params.dt;
nt = sim_params.nt;

r = 0; % cell radius
circ_rad = 1; % obstacle radius

trapped = 0;

% parse initial conditions
xi = xt_c(1,:,1);
thi = tht_c(1);
ui = [cos(thi),sin(thi)];
vi = beta*ui;
tumbler = 0;

neighbors_i = [];
lastNeighborCheck_pos = [];
contact_i = [];
if gamma<inf
    max_travel_R = circ_rad+r+5*beta*dt; 
    
    neighbors_i = updateNeighborsPerCell(xi,circ_pos',1,r,max_travel_R);
    lastNeighborCheck_pos = xi;

    contact_i = contacts_c(1);
end




% simulation loop
for tt = 2:nt

    % tumble or not
    tumbler = (RNs_tum_c(tt)<dt) & (~tumbler); % tumble rate lambda = 1. can't tumble if just tumbled
    tumbles_c(1,tt) = tumbler;

    % reorient tumbling cells
    if tumbler
        % generate angle
        thi = genTheta(1);

        ui = [cos(thi),sin(thi)];

        vi(1,:) = 0;
        vt_c(1,:,tt) = vi(1,:);
    end

    tht_c(1,tt) = thi;

    % obstacles
    if gamma<inf

        [xi,vi,contact_i] = updatePositionVelocity(xi,ui,dt,circ_pos,1+r,neighbors_i,beta,tumbler);
        % if tumbling, vi is set and contact_i is the
        % same as before

        contacts_c(1,tt) = contact_i;

    else
        if ~tumbler
            vi = beta*ui;
            xi = xi + vi*dt;
        end
    end

    % store updates
    xt_c(1,:,tt) = xi;
    vt_c(1,:,tt) = vi;

    % check whether to build new square or update neighbors
    if gamma<inf

        buildsquare = checkWhetherBuildSquare(xi,squares_pos,square_halfL,beta,dt,r);

        if buildsquare
            [squares_pos,~,circ_pos] = buildSquares(xi,squares_pos,circ_pos,square_halfL,rho);

            lastBuildtt = tt;

            % check if cell becomes trapped by this -- expensive
            if gamma<3
                trapped = checkWhetherTrapped(xi,circ_pos,circ_rad,square_halfL,squares_pos,rho);
            end
        end


        % check whether each cell is close to its last neighbor radius
        dr = sqrt(sum((xi-lastNeighborCheck_pos).^2,2));
        update = dr>max_travel_R-2*beta*dt | buildsquare;

        if update
            lastlastNeighborCheck_pos = lastNeighborCheck_pos;
            neighbors_iminus1 = neighbors_i;

            neighbors_i = updateNeighborsPerCell(xi,circ_pos',1,r,max_travel_R);%,maxN);

            lastNeighborCheck_pos = xi;
            lastUpdatett = tt;

        end

        ri = repmat(xi(:),1,length(neighbors_i)) - circ_pos(:,neighbors_i);
        phi_i = sum(ri.^2,1)'-1;

        if any(phi_i<-1000*eps*max(1,norm(xi)))
            keyboard
            disp('Cell entered an obstacle -- error.')
        end

    end

    if trapped
        disp('Cell became trapped -- breaking')
        break
    end

end

