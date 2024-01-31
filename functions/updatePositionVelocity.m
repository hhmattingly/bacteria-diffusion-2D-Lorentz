function [xtdt,vtdt,contacts] = updatePositionVelocity(xt,ut,dt,circ_pos,circ_rad,neighbors,beta,tumbler)
% [xtdt,vtdt,contacts] = updatePositionVelocity(xt,ut,dt,sph_pos,sph_rad,neighbors,v0,tumbler)
%
% Updates the cell's position, velocity, and contact state in the next dt.
%
% Inputs:
% xt: The cell's position at the end of the previous time step. 1 by d vector, where d is the number of dimensions
% ut: The cell's heading at the end of the previous time step. 1 by d vector, where d is the number of dimensions
% dt: Simulation time step. scalar
% circ_pos: Positions of the obstacle centers. d by Nircles matrix, where d is the number of dimensions
% circ_rad: Obstacle radius. Should always be 1. scalar
% neighbors: Vector of the indices of obstacles that are close to the cell. N by 1 vector
% beta: Dimensionless cell swimming speed. scalar
% tumbler: Whether the cell is tumbling this time step or not. logical scalar
%
% Outputs:
% xtdt: The cell's position at the end of this time step. 1 by d vector, where d is the number of dimensions
% vtdt: The cell's velocity at the end of this time step. 1 by d vector, where d is the number of dimensions
% contacts: The cell's contact state at the end of this time step. 0 indicates swimming in bulk, 1 indicates contact with one obstacle, 2 indicates contact with two obstacles. scalar
%
% Henry H. Mattingly, November 2023

tol = 1000*eps*max(1,norm(xt)); % tolerance
d=2;
xt = xt';
ut = ut';

circ_pos_i = circ_pos(:,neighbors); % centroids of nearby spheres
nclose = numel(neighbors);

dt_r = dt;
xNext = xt;


makeFigs = 0;
% for plotting %%%%%%%%%%%%%%%%
if makeFigs

figure;hold on
h=viscircles(circ_pos_i',circ_rad);
h.Children(1).Color = 'k';
plot(xt(1),xt(2),'co') % cell position
quiver(xt(1),xt(2),ut(1),ut(2),1,'filled') % cell heading
axis equal
xlabel('x')
ylabel('y')
h=gca;
h.XGrid='on';
h.YGrid='on';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start by checking state: how many obstacles in contact? 
[contacts,circInd] = getContacts(xt,ut,circ_pos_i,circ_rad,tol);

if tumbler
    xtdt = xt';
    vtdt = [0,0];
    return
end

% loop until the end of the time step
while dt_r > 0

    xLast = xNext;

    if contacts==0
        [xNext,vNext,contacts,dt_r,circInd] = update0(xLast,ut,beta,dt_r,circ_pos_i,circ_rad,circInd,tol);
        % if transition to state 1 occurs, contacts becomes 1, dt_r > 0

        % otherwise dt_r = 0

    elseif contacts==1

        [xNext,vNext,contacts,dt_r,circInd] = update1(xLast,ut,beta,dt_r,circInd,circ_pos_i,circ_rad,tol);
        % if no transition occurs, contacts = 1, dt_r = 0
        % if transition to state 0 occurs, contacts = 0, dt_r > 0
        % if transition to state 1' occurs, contacts = 1, dt_r > 0
        % if transition to state 2 occurs, contacts = 2, dt_r > 0

    elseif contacts==2
        % check that state is correct on first iteration (or do before the
        % while loop)...

        xNext = xLast;
        vNext = zeros(1,d);
        dt_r = 0;
        % done; cell doesn't move or change state when in state 2
    end


    % for plotting %%%%%%%%%%%%%%%%
    if makeFigs
        plot(xNext(1),xNext(2),'ko') % cell position
        quiver(xNext(1),xNext(2),ut(1),ut(2),'k')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if dt_r < 0
        keyboard
    end

end

xtdt = xNext;
vtdt = vNext;

rtdt = repmat(xtdt,1,nclose)-circ_pos_i;
phi_tdt = sum(rtdt.^2,1)' - circ_rad^2; % compute distances from cell position to the surfaces of all nearby obstacles

% for plotting %%%%%%%%%%%%%%%%
if makeFigs
    plot(xt(1),xt(2),'co') % initial cell position
    plot(xtdt(1),xtdt(2),'go') % cell position

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if norm(vtdt)>100*beta
    disp('warning: high speed')
    keyboard
end


if any(phi_tdt<-tol) %
    disp('warning: constraint above tolerance')
    keyboard
end

vtdt = vtdt(:)';
xtdt = xtdt(:)';


end




%% functions

%% state 0 update and related
function [xNext,vNext,contacts,dt_r,circInd] = update0(xLast,ut,beta,dt_r,circ_pos_i,circ_rad,circInd,tol)
nclose = size(circ_pos_i,2);
circInd = []; %% 

% compute time to next collision with all obstacles
[dContacts,tContacts] = getLinearTimeToContact(xLast,ut,beta,circ_pos_i,circ_rad,circInd,tol);

% exclude based on normals at contact point
if any(tContacts<=dt_r & tContacts>=0)
    checkInds = find(tContacts<=dt_r & tContacts>=0);
    for i = 1:length(checkInds)
        xContacti = xLast + dContacts(checkInds(i))*ut;
        ni = (xContacti - circ_pos_i(:,checkInds(i))); % outward
        dui = dot(ut,ni);
        if dui>=-tol
            tContacts(checkInds(i)) = Inf;
        end
    end
end

if any(tContacts<=dt_r & tContacts>=0) % trouble after slide-off escape
    [circInd,xNext,vNext,tContact] = getLinearContact(dContacts,tContacts,xLast,ut,beta,circ_pos_i,tol);
    dt_r = dt_r - tContact; % remaining time in dt. for later
    contacts = 1;
else
    vNext = beta*ut;
    xNext = xLast + beta*ut*dt_r;
    contacts = 0;
    dt_r = 0;
    circInd = [];
end

rNext = repmat(xNext,1,nclose) - circ_pos_i;
phiNext = sum(rNext.^2,1)'-circ_rad^2;

end

function [dContacts,tContacts] = getLinearTimeToContact(xt,ut,beta,circ_pos_i,circ_rad,circInd,tol)
nclose = size(circ_pos_i,2);
d=2;

y = repmat(xt,1,nclose)-circ_pos_i;
discrs = dot(repmat(ut,1,nclose),y,1).^2 - (sum(y.^2,1) - circ_rad^2);
d1 = -dot(repmat(ut,1,nclose),y,1) - sqrt(discrs); % must be negative root -- haven't hit yet, so closer intersection will be the - root
d2 = -dot(repmat(ut,1,nclose),y,1) + sqrt(discrs); % must be negative root -- haven't hit yet, so closer intersection will be the - root

% exclude case when cell is on a circle but pointing away
rt = repmat(xt,1,nclose) - circ_pos_i;
nt = rt./repmat(sqrt(sum(rt.^2,1)),d,1); % d

dun = dot(repmat(ut,1,nclose), nt);

d1(dun>-tol) = Inf; % -100*?
d2(dun>-tol) = Inf;
d1(d1<0) = Inf;
d2(d2<0) = Inf;
d1(imag(d1)~=0) = Inf;
d2(imag(d2)~=0) = Inf;
dContacts = min(d1,d2);

tContacts = dContacts/beta;
tContacts(circInd) = Inf;

end


function [circInd,xContact,vContact,tContact] = getLinearContact(ds,ts,xLast,ut,beta,circ_pos_i,tol)

[tContact,circInd] = min(ts);
xContact = xLast + ut*ds(circInd);

ri = xContact - circ_pos_i(:,circInd);
ni = -ri/norm(ri);
vContact = beta*(ut - dot(ut,ni)*ni);

end


%% state 1 update and related
function [xNext,vNext,contacts,dt_r,circInd] = update1(xLast,ut,beta,dt_r,circInd,circ_pos_i,circ_rad,tol)

% if no transition occurs, contacts = 1, dt_r = 0
% if transition to state 0 occurs, contacts = 0, dt_r > 0
% if transition to state 1' occurs, contacts = 1, dt_r > 0
% if transition to state 2 occurs, contacts = 2, dt_r > 0

nclose = size(circ_pos_i,2);

% compute time to next collision with all obstacles
[thContacts,tContacts,dthContacts] = getSurfaceTimeToContact(xLast,ut,beta,circInd,circ_pos_i,tol);
tContacts(isnan(tContacts)) = Inf;
tContacts(circInd) = Inf; %% exclude current sphere -- not a new encounter

% exclude based on normals at contact point -- ie previously contacted
% circle
if any(tContacts<=dt_r & tContacts>=0)
    checkInds = find(tContacts<=dt_r & tContacts>=0);
    for i = 1:length(checkInds)
        xContacti = circ_pos_i(:,circInd) + [cos(thContacts(checkInds(i)));sin(thContacts(checkInds(i)))];
        ni = (xContacti - circ_pos_i(:,checkInds(i))); % outward
        nCurrent = (xContacti - circ_pos_i(:,circInd)); % outward
        vt = ut - dot(ut,nCurrent)*nCurrent; % current v
    
        dui = dot(vt,ni); % including current sphere's force, cell won't hit sphere i
        if dui>=-tol
            tContacts(checkInds(i)) = Inf;
            dthContacts(checkInds(i)) = Inf;
        end
    end
end


[~,mnitc] = min(tContacts); 

% solve dynamics of position on the sphere until potentially tEscape
XLast = xLast - circ_pos_i(:,circInd);
th0 = mod(atan2(XLast(2),XLast(1)),2*pi);
phi = mod(atan2(ut(2),ut(1)),2*pi);
s0 = abs(cos(phi-th0));
tEscape = abs(1/(2*beta)*log((1+s0)/(1-s0)));
tEscape(tEscape<=tol) = Inf;

[~,mni] = min([dt_r,tEscape,tContacts(mnitc)]); % test2 and test5 kind of initiated by tEscape ~ tContacts...

if mni==1 % just slide
    % update position to sliding point
    thNext = phi - 2*acot(exp(dt_r*beta) * cot((phi-th0)/2));

    XNext = [cos(thNext); sin(thNext)];
    xNext = XNext + circ_pos_i(:,circInd);

    dt_r = 0;

    contacts = 1;
    ri = xNext - circ_pos_i(:,circInd);
    ni = ri/norm(ri);
    vNext = beta * (ut - (ut'*ni) * ni)';

elseif mni==2 % escape
    % update position to escape point
    thNext = phi - 2*acot(exp(tEscape*beta) * cot((phi-th0)/2));
    XNext = [cos(thNext); sin(thNext)];
    xNext = XNext + circ_pos_i(:,circInd);

    dt_r = dt_r - tEscape;

%     circInd = [];
    contacts = 0;
    vNext = beta*ut;

elseif mni==3 % collision
    % if there is a collision in that time, find where/when it happens.
    [circInd,xNext,tContact,vNext,contacts] = getSurfaceContact(thContacts,tContacts,dthContacts,ut,beta,circInd,circ_pos_i,tol);
    dt_r = dt_r - tContact;

end

% check for collision with a new sphere before escape/non-escape
rNext = repmat(xNext,1,nclose) - circ_pos_i;
phiNext = sum(rNext.^2,1)'-circ_rad^2;

end



function [thContacts,tContacts,dthContacts] = getSurfaceTimeToContact(xLast,ut,beta,circInd1,circ_pos_i,tol)

% finding a point where two spheres intersect
nclose = size(circ_pos_i,2);
sph_pos_local = circ_pos_i - repmat(circ_pos_i(:,circInd1),1,nclose);
xlocal = xLast - circ_pos_i(:,circInd1);
th0 = mod(atan2(xlocal(2),xlocal(1)),2*pi);
phi = mod(atan2(ut(2),ut(1)),2*pi);
thDir = sign(sin(phi-th0)); % which direction the cell is moving on the circle

% find intersection of two spheres. two solutions for theta
a = sph_pos_local(1,:);
b = sph_pos_local(2,:);
c = 1/2*(a.^2+b.^2);
discr = c.*(2 - c); % solutions only when <= 2
discr(discr<=0)=Inf;

th1 = atan2(b.*c + a.*sqrt(discr), a.*c - b.*sqrt(discr)); % (y,x)
th2 = atan2(b.*c - a.*sqrt(discr), a.*c + b.*sqrt(discr)); % (y,x)

% shift s.t. th0 is at zero
dth1 = th1-th0;
dth2 = th2-th0;

% now want the one that goes in direction thDir and has smallest magnitude.
% correct for direction.
dth1 = thDir*dth1; dth1 = mod(dth1+tol,2*pi)-tol;
dth2 = thDir*dth2; dth2 = mod(dth2+tol,2*pi)-tol;

% remove invalid points
dth1(discr==Inf) = Inf;
dth2(discr==Inf) = Inf;
dth1(dth1>pi/2) = Inf;
dth2(dth2>pi/2) = Inf;

%
dthContacts = min(dth1,dth2);
thContacts = mod(th0 + thDir*dthContacts,2*pi);

% invert solution for s = cos(phi-theta)
s0 = cos(phi-th0);
sContacts = cos(phi-thContacts);
tContacts = abs(1/(2*beta) * log( (s0+1)/(s0-1) * (sContacts - 1)./(sContacts + 1)));
tContacts(thContacts==Inf) = Inf;

end

function [circInd2,xContact,tContact,vContact,contacts] = getSurfaceContact(thContacts,tContacts,dthContacts,ut,beta,circInd1,circ_pos_i,tol)

[~,circInd2] = min(dthContacts); % will regret this?
tContact = tContacts(circInd2);
thContact = thContacts(circInd2);

xContact = [cos(thContact);sin(thContact)] + circ_pos_i(:,circInd1);

% does the cell keep sliding along the next sphere, or does it stop?
n1 = xContact - circ_pos_i(:,circInd1);
n2 = xContact - circ_pos_i(:,circInd2);

dsphs = sqrt(sum((circ_pos_i(:,circInd1) - circ_pos_i(:,circInd2)).^2,1));

n1=-n1;
n2=-n2;

d12 = dot(n1,n2);
du1 = dot(ut,n1);
du2 = dot(ut,n2);

if dsphs<sqrt(2) % NEW
    c = [0;du2];

    % check result
    vt = ut  - [n1,n2]*c;

    dv1 = dot(vt,n1);
    if dv1<0 % for inward pointing normal
        contacts = 1;
    else
        c = [1,d12;d12,1]\[du1;du2]; 
        contacts = 2;
    end
else
    % dot with inward pointing normals
    % v = u - c1 n1 - c2 n2
    % v.n1 = u.n1 - c1 n1.n1 - c2 n2.n1 = 0
    % v.n2 = u.n2 - c1 n1.n2 - c2 n2.n2 = 0
    % u.n1 - c1 - c2 n2.n1 = 0
    % u.n2 - c1 n1.n2 - c2 = 0
    % u.n1 = [1, n2.n1] * [c1;c2]
    % u.n2 = [n2.n1, 1] * [c1;c2]

    c = [1,d12;d12,1]\[du1;du2];

    contacts = 2;
end

% get speed
vt = beta*(ut - [n1,n2]*c);

vContact = vt;

end


%% get state
function [contacts,circInd] = getContacts(xt,ut,circ_pos_i,circ_rad,tol)
nclose = size(circ_pos_i,2);

rt = repmat(xt,1,nclose) - circ_pos_i;
phi_t = sum(rt.^2,1)'-circ_rad^2;
contacts = sum(abs(phi_t)<=tol);
circInd = [];

% verify initial state
if contacts == 1
    circInd = find(abs(phi_t)<=tol);

    % check that the cell really is in state 1 -- ie check that it's
    % pointing into a sphere -- relevant only to the first iteration
    ni = -(xt - circ_pos_i(:,circInd));
%     ni = ri/norm(ri); % unnecessary for unit circles

    % dot with inward pointing normal
    if dot(ut,ni)<=0
        % if <=0, then actually in state 0
        contacts = 0;
        circInd = [];
    end

elseif contacts == 2
    
    sphInds = find(abs(phi_t)<=tol);
    n1 = xt - circ_pos_i(:,sphInds(1));
    n2 = xt - circ_pos_i(:,sphInds(2));

    % dot with inward pointing normals
    n1=-n1;
    n2=-n2;
    du1 = dot(ut,n1);
    du2 = dot(ut,n2);

    if du1<=0 && du2<=0
        contacts = 0; % facing away from spheres

    else
        % remove projection of u onto each sphere normal
        v1 = ut - n1*du1;
        v2 = ut - n2*du2;

        % dot resulting velocity with other sphere's inward normal
        cross_dots = [dot(v1,n2), dot(v2,n1)];

        if any(cross_dots<0)
            % if either dot product <0, cell has a way out and is in state 1
            % find which sphere it slides along next
            if dot(v1,n2)<0
                circInd = sphInds(1); % then cell points away from n2 (inward normal to sphere 2), and thus can move on sphere 1
            else
                circInd = sphInds(2); % opposite case
            end
            contacts = 1;
        else
            contacts = 2;
        end

    end
end

end


