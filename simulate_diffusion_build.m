%% Main script for simulating a run-and-tumble particle in a 2D environment with circular obstacles.
% Henry H. Mattingly, August 2024

addpath(['.' filesep 'functions'])

%% Parameters that don't vary across simulations
% circle radius rescaled to 1
% tumble rate rescaled to 1

Nreps = 10; % # of replicate simulations
Ncells = 50; % cells per replicate

% function that generates new headings after tumbles
genTheta = @(N) 2*pi*rand(N,1);

d = 2; % # of dimensions
env_params.d = d;

%% Parameters that do vary across simulations
gammas = 10.^[1/4:1/4:6/4]; % dimensionless mean chord length
betas = 10.^[-1:1/2:6/2]; % dimensionless swimming speed

rhos = 1./(2*gammas); % dimensionless obstacle density, average # of obstacle centers per area
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

            % cases to skip
            if gamma_ind==1&&beta_ind>=9
                continue
            end

             % store parameters
            cell_params.beta = beta;
            env_params.gamma = gamma;

            %% time
            % simulation time, in units of the average run duration
            T=2000;

            dt = 1/50; % time step
            t = 0:dt:T;
            nt = length(t);
            t1 = ones(Ncells,1);

            sim_params.T = T;
            sim_params.dt = dt;
            sim_params.nt = nt;

            %% for building the environment
            square_halfL = max(5,3*beta*dt); % if obstacles are dense, can be way smaller
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

            % to skip or continue simulations that were run
            continue_flag = 0; % if continue_flag = 1, resume unfinished simulations
            file_name = [save_dir filesep 'simdata_rep' sprintf('%02u',rep) '.mat'];
            
            if exist(file_name,'file')
                disp('File exists...')
                s=load(file_name,'closed','Ncells','nt');

                if ~any(s.closed) && s.Ncells>=Ncells && s.nt>=nt %
                    disp('All simulations run. Skipping.')
                    continue
                elseif s.Ncells<Ncells || s.nt<nt
                    disp('Continuing simulations.')
                    continue_flag = 1;
                elseif any(s.closed)
                    disp('Continuing simulations.')
                    continue_flag = 1;
                else
                    disp('Starting over simulations--data will be overwritten.')
                end
            end


            %% set up simulation arrays
            [xt,tht,vt,tumbles,contacts,squares,circles,circInds,closed,nbreak,nsims,RNs_tum] = allocateArrays(Ncells,d,nt,continue_flag,file_name);


            %% run simulations that haven't been run yet, or that ended up in closed pores
            while any(closed)
                inds = find(closed);
    
                for i = 1:length(inds)
                    cell_ind = inds(i);

                    disp(['Simulating cell number ' num2str(cell_ind)])

                    if t1(cell_ind)==1
                        % generate realization of environment and ICs
                        disp('Generating environment...')
                        [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circInds_c,circ_pos,squares_pos] = generateICs(cell_params,env_params,sim_params,genTheta);
                    else
                        xt_c = xt(cell_ind,:,:);
                        vt_c = vt(cell_ind,:,:);
                        tht_c = tht(cell_ind,:);
                        tumbles_c = tumbles(cell_ind,:);
                        contacts_c = contacts(cell_ind,:);
                        circInds_c = circInds(cell_ind,:,:);
                        circ_pos = circles{cell_ind};
                        squares_pos = squares{cell_ind};
                    end
                    

                    %% simulate
                    disp('Simulating...')

                    tic
                    [xt_c,vt_c,tht_c,tumbles_c,contacts_c,circInds_c,circ_pos,squares_pos,trapped] = mainSimLoop(t1(cell_ind),xt_c,vt_c,tht_c,tumbles_c,contacts_c,circInds_c,cell_params,env_params,sim_params,circ_pos,squares_pos,RNs_tum(cell_ind,:),genTheta);
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
                        circInds(cell_ind,:,:) = circInds_c;
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

            for cell_num = 1:5%0%Ncells%cell_ind %find(~closed,1,'first')
            
                options.showRectangles=0;
                options.showGrid=0;
                options.showDots=0;
                options.quickPlot=1;
                options.showTumbles=0;
                options.showContacts=0;
                visualize_traj(cell_num,nt,circles,squares,square_halfL,xt,tumbles,contacts,vt,options)

            end

            drawnow

            %% MSD
            disp('Computing MSD')

            dtmsd = T/2/50; % sub-sample the MSD
            every_n = round(dtmsd/dt);
            tic
            [tmsd,msdr,msdr_N,msdr_std] = compute_MSD(xt, dt, T, every_n);
            toc
            msdr_mean = nanmean(msdr,1);


            %% Plots
            figure;hold on
            errorbar(tmsd,msdr_mean,nanstd(msdr,[],1)/sqrt(Ncells))
            plot(t,2*d*D0*t*(1-circAreaFrac))
            xlabel('Time')
            ylabel('MSD \langle\Deltar^2\rangle')
            hleg=legend('MSD','2 d D_{liquid} \phi_{void} t','Location','northwest');
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
            save(file_name)
            toc

        end
    end
end