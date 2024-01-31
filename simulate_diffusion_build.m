%% Main script for simulating a run-and-tumble particle in a 2D environment with circular obstacles.
% Henry H. Mattingly, November 2023

addpath(['.' filesep 'functions'])

%% Parameters that don't vary across simulations
Nreps = 1; % # of replicate simulations -- 3 in main text
Ncells = 2; % cells per replicate -- 250 in main text

% function that generates new headings after tumbles
genTheta = @(N) 2*pi*rand(N,1);

% simulation time, in units of the average run duration
T = 20; % 2000 in main text
dt = 1/50; % time step
t = 0:dt:T;
nt = length(t);

sim_params.T = T;
sim_params.dt = dt;
sim_params.nt = nt;

d = 2; % # of dimensions
env_params.d = d;

%% Parameters that do vary across simulations
% circle radius rescaled to 1

% gamma = L/R
gammas = 10.^[1/4:1/4:6/4]; % dimensionless mean chord length
betas = 10.^[-1:1/2:6/2]; % dimensionless swimming speed

rhos = 1./(2*gammas); % obstacle density, average # per area
etas = pi*rhos; % "reduced density"

circAreaFracs = 1-exp(-etas); % obstacle area fraction
voidAreaFracs = exp(-etas); % void area fraction

Ngammas = length(gammas);
Nbetas = length(betas);


%% Main loop
for rep = 1:Nreps
    for beta_ind = 1:Nbetas
        for gamma_ind = 1:Ngammas 

            % get parameters
            beta = betas(beta_ind);
            gamma = gammas(gamma_ind);

            % skipped cases
            if gamma_ind>=6&&beta_ind<3 
                continue
            end

            if gamma_ind==1 && beta_ind>3 
                continue
            end

            if (gamma_ind<2 || gamma_ind>4) && beta_ind>7 
               continue 
            end

            if beta>=500*gamma
                continue
            end
            
            % store parameters
            cell_params.beta = beta;
            env_params.gamma = gamma;

            % for building the environment
            square_halfL = max(5,beta*dt*5);
            env_params.square_halfL = square_halfL;

            circAreaFrac = circAreaFracs(gamma_ind);
            rho = rhos(gamma_ind);

            env_params.circAreaFrac = circAreaFrac;
            env_params.rho = rho;


            %% saving/loading/skipping
            save_dir = ['.', filesep, 'sim_data', filesep, 'gamma=' num2str(round(gamma,2)) '_beta=' num2str(round(beta,2))];

            if ~exist(save_dir,'dir')
                mkdir(save_dir)
            end

            disp(['Gamma = ' num2str(gamma) ', Beta = ' num2str(beta) ', replicate ' num2str(rep)])

%             to skip simulations that were run
            if exist([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'],'file')
                disp('File exists, skipping.')
                continue
            end

            % to load data
%             if exist([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'],'file')
%                 disp('Loading to simulate')
%                 load([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'])
%             end
% 
            %% set up simulation arrays
            xt = nan(Ncells,d,nt);
            tht = nan(Ncells,nt);
            vt = nan(Ncells,d,nt);

            tumbles = zeros(Ncells,length(t));
            contacts = nan(Ncells,length(t));
            contacts(:,1) = 0;

            % random numbers for tumbles
            RNs_tum = rand(Ncells,nt);

            % record the cubes and spheres
            squares = cell(Ncells,1);
            circles = cell(Ncells,1);

            % was environment closed at the beginning/end of the simulation?
            closed = ones(Ncells,1);

            %% run simulations that haven't been run yet, or that ended up in closed pores
            nbreak = 0;
            nsims = 0;

            %%
            while any(closed)
                inds = find(closed);
    
                for i = 1:length(inds)
                    cell_ind = inds(i);

                    disp(['Simulating cell number ' num2str(cell_ind)])

                    % generate realization of environment and ICs
                    disp('Generating environment...')

                    [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos] = generateICs(cell_params,env_params,sim_params,genTheta);

                    %% simulate
                    disp('Simulating...')

                    tic
                    [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circ_pos,squares_pos,trapped] = mainSimLoop(xt_c,vt_c,tht_c,tumbles_c,contacts_c,cell_params,env_params,sim_params,circ_pos,squares_pos,RNs_tum(cell_ind,:),genTheta);
                    toc

                    % store
                    xt(cell_ind,:,:) = xt_c;
                    vt(cell_ind,:,:) = vt_c;
                    tht(cell_ind,:) = tht_c;
                    tumbles(cell_ind,:) = tumbles_c;

                    if gamma<inf
                        contacts(cell_ind,:) = contacts_c;
                        squares{cell_ind} = squares_pos;
                        circles{cell_ind} = circ_pos;
                    end

                    closed(cell_ind) = trapped; 

                    nsims = nsims+1;

                    if trapped==1
                        nbreak = nbreak+1;

                        disp([num2str(nbreak/nsims) ' broken'])
                    end
                end

            end

            %% Visualize outputs

%             for cell_num = cell_ind %find(~closed,1,'first')
%             
%                 options.showRectangles=1;
%                 options.showGrid=1;
%                 visualize_traj(cell_num,nt,circles,squares,square_halfL,xt,tumbles,contacts,vt,options)
% 
%             end
% 
%             drawnow

            %% MSD
            disp('Computing MSD')

            dtmsd = T/2/50; % sub-sample the MSD
            every_n = round(dtmsd/dt);
            tic
            [tmsd,msdr,msdr_N,msdr_std] = compute_MSD(xt, dt, T, every_n);
            toc
            msdr_mean = nanmean(msdr,1);

            %% Plots
            D0 = beta^2/d;

            figure;hold on
            errorbar(tmsd,msdr_mean,nanstd(msdr,[],1)/sqrt(Ncells))
            plot(t,2*d*D0*t*(1-circAreaFrac))
            xlabel('Time')
            ylabel('MSD \langle\Deltax^2\rangle')
            hleg=legend('MSD','D_{liquid} \phi_{void}','Location','northwest');
            hleg.Box = 'off';
            h=gca;h.Box='off';
            h.YLim(1) = 0;

            drawnow

            %% delete unnecessary variable before saving
            clear RNs_tum
            
            %%
            close all

            %%
            disp('Saving data...')
            tic
            save([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'])
            toc

        end
    end
end
