function trapped = checkWhetherTrapped(x0,circ_pos,circ_rad,square_halfL,squares_pos,rho)
% trapped = checkWhetherTrapped(x0,circ_pos,sph_rad,square_halfL,squares_pos,rho)
%
% Check whether a cell is completely enclosed by circular obstacles.
%
% Inputs:
% x0: Current cell position. 1 by d vector, where d is the number of spatial dimensiones
% circ_pos: Positions of the centers of previously-placed circular obstacles. d by Ncircles matrix.
% circ_rad: Obstacle radius. Should equal 1. scalar
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar.
% squares_pos: Positions of the centers of squares comprising the simulation domain. d by Nsquare matrix.
% rho: Dimensionless number density of obstacles (per area). scalar.
%
% Outputs:
% trapped: Whether or not the cell is trapped by obstacles. logical scalar.
%
% Henry H. Mattingly, November 2023

% x0 = x0';
if size(x0,1)~=2
    x0 = x0';
end
nt = size(x0,2);
if nt>1
    x0 = [x0(:,1:round(nt/100):end), x0(:,end)];
end

[d,ncircles] = size(circ_pos);
nsquares = size(squares_pos,2);

% test code
% Npercircle = 1000;
% th = linspace(0,2*pi,Npercircle);
% xy = [cos(th);sin(th)];
% 
% pts = repmat(xy,1,ncircles);
% circleID = zeros(1,size(pts,2));
% for i = 1:ncircles
%     pts(:,1+(i-1)*Npercircle:i*Npercircle) = pts(:,1+(i-1)*Npercircle:i*Npercircle) + repmat(circ_pos(:,i),1,Npercircle);
%     circleID(1+(i-1)*Npercircle:i*Npercircle) = i;
% end
% 
% %%%%%%%%%%%%%%%%%%%
% makefigs = 0;
% if makefigs
%     figure;
%     hold on
%     viscircles(circ_pos',circ_rad,'color','r');
%     for i = 1:nsquares
%         rectangle('Position',[squares_pos(:,i)'-square_halfL*ones(1,2),2*square_halfL*ones(1,2)]);
%     end
%     plot(pts(1,:),pts(2,:),'k.')
% 
%     plot(x0(1),x0(2),'co')
%     
%     pause(0.1)
%     drawnow
% 
% end
% %%%%%%%%%%%%%%%%%%%
% 
% Npts = size(pts,2);
% D = pdist2(pts',circ_pos');
% 
% % A = zeros(Npts);
% % [I,J]=find(D<circ_rad);
% % for i = 1:Npts
% %     ID = circleID(I(i));
% %     A(ID,J(i)) = 1;
% %     A(J(i),ID) = 1;
% % end
% 
% % D = min(D,[],2); % only care about closest circle center for each point
% rmInds = any(D<circ_rad,2);
% pts(:,rmInds) = [];
% circleID(rmInds) = [];





rho_pts = 150*rho;
npts_per_square = ceil(rho_pts*(2*square_halfL)^d);
% average distance between samples:
% prob of no pts in circle is exp(-rho*A) -- survival probability
% -d/dr of that is P(r) = rho*dA/dr * exp(-rho*A) = 2*pi*r*rho * exp(-rho*A)
% rate rho*dA/dr = rho*2*pi*r -- depends on r.
% mean value is 1/\sqrt(rho)
% fraction of samples A that lie within distance r = sqrt(log(1/(1-A))/(pi rho))
A = 1-10^-8;
rthresh = sqrt(log(1/(1-A))/(pi*rho_pts)); 


samples = [];
edge_inds = [];
for i = 1:nsquares
    samples_i = 2*square_halfL*(rand(d,npts_per_square)-1/2);
    % throw out those that fall in other squares

    edge_inds_i = abs(samples_i(1,:)-square_halfL)<=rthresh | abs(samples_i(2,:)-square_halfL)<=rthresh | abs(samples_i(1,:)+square_halfL)<=rthresh | abs(samples_i(2,:)+square_halfL)<=rthresh;

    samples_i =  samples_i + squares_pos(:,i);

    if nsquares>1 && i>1
        other_squares = squares_pos(:,1:i-1);
        
        rmInds = any(all(abs(repmat(samples_i,[1,1,i-1])-repmat(permute(other_squares,[1,3,2]),[1,npts_per_square,1]))<square_halfL,1),3);
        
        samples_i(:,rmInds) = [];
        edge_inds_i(rmInds) = [];
    end
    
    edge_inds_i = find(edge_inds_i);
    if nsquares>1
        other_squares = squares_pos(:,setdiff(1:nsquares,i));
        rmInds = find(any(all(abs(repmat(samples_i(:,edge_inds_i),[1,1,nsquares-1])-repmat(permute(other_squares,[1,3,2]),[1,length(edge_inds_i),1]))<square_halfL,1),3));
        edge_inds_i(rmInds) = [];
    end
    
    edge_inds = [edge_inds; length(samples) + edge_inds_i(:)];
    samples = [samples, samples_i];
end


% throw out samples inside circles
D = pdist2(samples',circ_pos');
D = min(D,[],2); % only care about closest circle center for each point
rmInds = D'<circ_rad;
temp = zeros(1,size(samples,2)); temp(edge_inds) = 1;
samples(:,rmInds) = [];
temp(rmInds) = [];
edge_inds = find(temp); edge_inds = edge_inds(:);

Dpts = pdist2(x0',samples'); % x0 can be a vector. %sqrt(sum((repmat(x0,1,size(samples,2)) - samples).^2,1)); %pdist2(x0',samples');
Dpts = any(Dpts<rthresh,1);

new_pts = find(Dpts); new_pts = new_pts(:);
all_pts = new_pts;
while ~isempty(new_pts) && ~any(ismember(edge_inds,new_pts))
    for i = 1:length(new_pts)
%         Dpts = pdist2(samples(:,new_pts(i))',samples'); % slow step
        Dpts = sqrt(sum((repmat(samples(:,new_pts(i)),1,size(samples,2)) - samples).^2,1));
        Dpts = Dpts<rthresh;
        inds = find(Dpts);
        new_pts = [new_pts; inds(:)];
    end
    new_pts = unique(new_pts);
    new_pts = setdiff(new_pts,all_pts);
    all_pts = [all_pts;new_pts];

%     if any(ismember(edge_inds,new_pts))
%         break
%     end
end

% do x0 and edge samples lie in the same graph?
trapped = all(ismember(edge_inds,all_pts)==0);

%%%%%%%%%%%%%%%%%%%
makefigs = 0;
if makefigs
    figure;
    hold on
    viscircles(circ_pos',circ_rad,'color','k');
    for i = 1:nsquares
        rectangle('Position',[squares_pos(:,i)'-square_halfL*ones(1,2),2*square_halfL*ones(1,2)]);
    end
    plot(samples(1,:),samples(2,:),'.')
    plot(samples(1,edge_inds),samples(2,edge_inds),'g.')
    plot(x0(1),x0(2),'ro')
    plot(samples(1,all_pts),samples(2,all_pts),'c.')
    
    pause(0.1)
    drawnow

end
%%%%%%%%%%%%%%%%%%%
end