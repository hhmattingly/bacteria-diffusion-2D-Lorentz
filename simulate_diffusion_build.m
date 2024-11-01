%% Main script for simulating a run-and-tumble particle in a 2D environment with circular obstacles.
% Henry H. Mattingly, August 2024

addpath(['.' filesep 'functions'])

%% Parameters that don't vary across simulations
% circle radius rescaled to 1
% tumble rate rescaled to 1

Nreps = 10; % # of replicate simulations -- 3 in first submission
Ncells = 50; % cells per replicate -- 250 in first submission
% 750 total

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


% b100 g10 rep4; b300 g10 rep4
% run b0.32 for longer, esp g5 and g10. maybe increase maximum. fewer total. 10 reps should
% be plenty.

%%
% for beta_ind = 1:Nbetas
% for gamma_ind = 1:Ngammas
% beta = betas(beta_ind);
% gamma = gammas(gamma_ind);
% 
% beta
% gamma
% gamma/beta
% 
% % before
% if beta>10 || (beta==0.1 && gamma>2) || gamma>10 % shorten beta<1, gamma>10. change condition to gamma/beta > X
% T = min(max(2000,2000/beta*gamma/10^(1/4)),20000); % 2000 in main text
% else
% T = min(max(2000,2000*(gamma/10^(1/4)/beta)^2),30000); % 2000 in main text
% end
% 
% T
% 
% % proposed
% if beta>10 || (beta==0.1 && gamma>2) || gamma/beta>10 % shorten beta<1, gamma>10. change condition to gamma/beta > X
% T = min(max(2000,2000/beta*gamma/10^(1/4)),20000); % 2000 in main text
% else
% T = min(max(2000,2000*(gamma/10^(1/4)/beta)^2),30000); % 2000 in main text
% end
% 
% T
% 
% keyboard
% end
% end

% 3 3 is only 2000?
% for gamma_ind = 2:3
%     for beta_ind = 1:5%Nbetas
% 
%         beta = betas(beta_ind)
%         gamma = gammas(gamma_ind)
% 
%         if beta>10 || (beta==0.1 && gamma>2) || gamma>10 % shorten beta<1, gamma>10. change condition to gamma/beta > 10?
%             T = min(max(2000,2000*gamma/10^(1/4)/beta),20000); % 2000 in main text
%         else
%             T = min(max(10000,2000*(gamma/10^(1/4)/beta)^2),30000); % 2000 in main text
%         end
%         T
%     end
% end





%% Main loop
for rep = Nreps:-1:1
    for beta_ind = 2:5
        for gamma_ind = 2:3%Ngammas

            % get parameters
            beta = betas(beta_ind);
            gamma = gammas(gamma_ind);

            if gamma_ind==1&&beta_ind>=9
                continue
            end
%             if beta>300 
%                 continue
%             end

             % store parameters
            cell_params.beta = beta;
            env_params.gamma = gamma;

            %% time
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % simulation time, in units of the average run duration
            if beta>10 || (beta==0.1 && gamma>2) || gamma>10 % shorten beta<1, gamma>10. change condition to gamma/beta > 10?
                T = min(max(2000,2000*gamma/10^(1/4)/beta),20000); % 2000 in main text
            else
                T = min(max(10000,2000*(gamma/10^(1/4)/beta)^2),30000); % 2000 in main text
            end
%             T=2000;

            dt = 1/50; % time step
            t = 0:dt:T;
            nt = length(t);
            t1 = ones(Ncells,1);

            sim_params.T = T;
            sim_params.dt = dt;
            sim_params.nt = nt;

            %% for building the environment
            square_halfL = max(5,3*beta*dt); %%%% 3? 4? if dense, can be way smaller. but how dense is dense enough?
            % smaller for gamma < 2,3?
            % any problem when square_halfL>5? beta > ~83?
            % update neighbor issues?

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
            continue_flag = 0;
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

            T

            % to load data
%             if exist([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'],'file')
%                 disp('Loading to simulate')
%                 load([save_dir filesep 'simdata_rep0' num2str(rep) '.mat'])
%             end
% 
            %% set up simulation arrays
            if continue_flag==0
                xt = nan(Ncells,d,nt);
                tht = nan(Ncells,nt);
                vt = nan(Ncells,d,nt);
            
                tumbles = zeros(Ncells,nt);
                contacts = nan(Ncells,nt);
                contacts(:,1) = 0;
            
                % record the cubes and spheres
                squares = cell(Ncells,1);
                circles = cell(Ncells,1);
                circInds = nan(Ncells,d,nt);
            
                % was environment closed at the beginning/end of the simulation?
                closed = ones(Ncells,1);
            
                nbreak = 0;
                nsims = 0;
            else
                load(file_name,'xt','tht','vt','tumbles','contacts','squares','circles','circInds','closed','nbreak','nsims')
                Nadd=0;
                nt1 = size(xt,3);
            
                if Ncells>size(xt,1)
                    Nadd = Ncells-size(xt,1);
            
                    xt = cat(1,xt,nan(Nadd,d,nt1));
                    tht = [tht; nan(Nadd,nt1)];
                    vt = cat(1,vt,nan(Nadd,d,nt1));
            
                    tumbles = [tumbles; zeros(Nadd,nt1)];
                    contacts = [contacts; nan(Nadd,nt1)];
            
                    % record the cubes and spheres
                    squares = [squares; cell(Nadd,1)];
                    circles = [circles; cell(Nadd,1)];
                    circInds = cat(1, circInds, nan(Nadd,d,nt1));
            
                    % was environment closed at the beginning/end of the simulation?
                    closed = [closed; ones(Nadd,1)];
            
                    Ncells = Ncells + Nadd;
                end
            
                %
                if nt>nt1
                    ntadd = nt-nt1;
                    t1 = nt1*ones(Ncells,1);
            
                    if Nadd>0
                        t1(Ncells-Nadd+1:end) = 1;
                    end
            
                    xt = cat(3,xt,nan(Ncells,d,ntadd));
                    tht = [tht, nan(Ncells,ntadd)];
                    vt = cat(3,vt,nan(Ncells,d,ntadd));
            
                    tumbles = [tumbles, zeros(Ncells,ntadd)];
                    contacts = [contacts, nan(Ncells,ntadd)];
                    circInds = cat(3, circInds, nan(Ncells,d,ntadd));
            
                    % was environment closed at the beginning/end of the simulation?
                    closed = ones(Ncells,1);
                end
            end

            % random numbers for tumbles
            RNs_tum = rand(Ncells,nt);

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

            %% k01 -- sanity check
            N0 = 0;
            N01 = 0;
            theta0 = [];

            for cell_ind = 1:Ncells
                contacts_c = contacts(cell_ind,:);
                circInds_c = squeeze(circInds(cell_ind,:,:));
                circ_pos = circles{cell_ind};
                xt_c = squeeze(xt(cell_ind,:,:));
                tht_c = tht(cell_ind,:);
                ut_c = [cos(tht_c); sin(tht_c)];

                state0inds = find(contacts_c(1:end-1)==0);
                state01inds = find(contacts_c(1:end-1)==0 & contacts_c(2:end)>0);

                N0 = N0 + numel(state0inds);
                N01 = N01 + numel(state01inds);

            end
            k01 = N01/N0/dt / beta

            gamma^-1

            gamma^(-1)/k01

        



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



%% theory 
phiv = voidAreaFracs(gamma_ind);

% parameters -- what do simulations say about these? nu, k10, v12?
nu_0 = 1/2; % large gamma, but also ignoring perpendicular displacement and effect of 0th order forcing 1st theta eqn
nu_0_perp = 1/pi;
nu_1 = 1/2;
nu_1_perp = 1/(3*pi);
k1inf = 2/pi; % large gamma, beta
k1 = k1inf;
k1_1_inf = 2/pi;
k1_1 = k1_1_inf;
k1_1_perp = 1;
v12_1_perp = pi/8 + 1/6;

dilute = 0; % if = 1, constant nu, k10, k12

p10 = 1/2;
p11 = 1-p10;
v12inf = 3/2/pi; v12 = v12inf; v12_1 = v12inf;
AR2 = 0.2215;%0.158;
p21 = 1/2;
p22 = 0.322;%16;
p20 = 1/2-p22;

% catalan's constant
G = 0.91596559417721901505460351493238;

% theory
ginv = 1./gamma;
binv = 1./beta;

% zeroth order quantities
BB = 2/pi*p21/(1-p22) * (v12inf/k1inf + AR2);
Z = pi + ginv .* (pi^2/16 - 1/2 * log(8) + 4*G.* BB);

k1 = 2./Z.*(1 - ginv.*(3/4 - pi/2*BB));
nu_0 = 1/2 + ginv*(log(8)/(4*pi) - 7/(8*pi) + (1-G)*pi/2*BB); % just to gamma^-1

k12p = ginv.*(v12+k1*AR2);
b = k12p./(binv*(1-p22));
a = ginv ./ (binv*p10 + k12p .* (p20/(1-p22)) + k1.*(1-ginv*AR2) );
%     a = ginv ./ (binv*p10 + k12p .* (p20/(1-p22)) + k1.*(1-ginv*AR2) - k1.*(1-exp(-ginv*(pi+6*sqrt(3))/12/2)));
p0_th = 1./(1 + a + a.*b);
p1_th = a./(1 + a + a.*b);
p2_th = a.*b./(1 + a + a.*b);

% first order quantities
k1_1 = 2/pi - ginv*(1/8 + 3/(2*pi) -1/pi^2*log(8));
Lambda1 = binv + k1_1 + v12_1*ginv;
Lambda0 = binv + ginv .* (binv + v12_1*ginv + ginv*AR2*k1_1)./Lambda1;

Lambda1_perp = binv + k1_1_perp + v12_1_perp*ginv;

chi_0 = 1./Lambda0 .* (p0_th + nu_0.*(1-ginv*AR2).*k1_1./Lambda1 .* p1_th);

nu_1 = 1/2 + ginv*(log(8)/(4*pi) - 7/(8*pi));

J0_term = 1/d*chi_0;
J1_parallel_term = 1/d*nu_1./Lambda1.*( ginv./Lambda0.*p0_th + (1+(1-ginv*AR2).*ginv./Lambda0.*k1_1./Lambda1).*nu_0.*p1_th );
J1_perp_term = 1/d*1./Lambda1_perp.*nu_0_perp.*nu_1_perp.*p1_th;

D_eff = phiv*(J0_term + J1_parallel_term + J1_perp_term);

D_sim_units = beta*D_eff;



figure;hold on
errorbar(tmsd,msdr_mean,nanstd(msdr,[],1)/sqrt(Ncells))
plot(t,2*d*D0*t*(1-circAreaFrac))
xlabel('Time')
ylabel('MSD \langle\Deltar^2\rangle')
hleg=legend('MSD','2 d D_{liquid} \phi_{void} t','Location','northwest');
hleg.Box = 'off';
h=gca;h.Box='off';
h.YLim(1) = 0;

plot(tmsd,2*d*tmsd*D_sim_units)