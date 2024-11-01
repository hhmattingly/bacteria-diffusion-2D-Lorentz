function [contacts,circInds,nearbyCircInds,varargout] = getContacts(xt,ut,circ_pos_i,circ_rad,tol,varargin)
% INPUTS

nclose = size(circ_pos_i,2);

xt = xt(:);
ut = ut(:);

if nargin>=6
    phi_t = varargin{1};
else
    rt = repmat(xt,1,nclose) - circ_pos_i;
    phi_t = sum(rt.^2,1)'-circ_rad^2;
end
contacts = sum(abs(phi_t)<=tol);
circInds = [];
nearbyCircInds = [];

% verify initial state
if contacts == 1
    circInds = find(abs(phi_t)<=tol);

    % check that the cell really is in state 1 -- ie check that it's
    % pointing into a sphere -- relevant only to the first iteration
    ni = -(xt - circ_pos_i(:,circInds));
    %     ni = ri/norm(ri); % unnecessary for unit circles

    % dot with inward pointing normal
    if dot(ut,ni)<=0
        % if <=0, then actually in state 0
        contacts = 0;
        circInds = [];
        nearbyCircInds = [];
    end

elseif contacts == 2

    nearbyCircInds = find(abs(phi_t)<=tol);
    circInds = nearbyCircInds;
    n1 = xt - circ_pos_i(:,nearbyCircInds(1));
    n2 = xt - circ_pos_i(:,nearbyCircInds(2));

    % dot with inward pointing normals
    n1=-n1;
    n2=-n2;
    du1 = dot(ut,n1);
    du2 = dot(ut,n2);
    % positive dot product indicates pointing into that obstacle.

    % check rotation directions...
    % if ut lies within n1 and n2, then both of the following are true:
    % 1) rotation direction from n1 to ut is the same as the rotation
    % direction from n1 to n2
    % 2) rotation direction from n2 to ut is the same as the rotation
    % direction from n2 to n1

    % cross products
    n1xut = n1(1)*ut(2) - n1(2)*ut(1); % chosen so that positive means rotate n1 in +theta direction to get to ut
    n1xn2 = n1(1)*n2(2) - n1(2)*n2(1);
    n2xut = n2(1)*ut(2) - n2(2)*ut(1);
    n2xn1 = -n1xn2;

    if du1<=0 && du2<=0
        contacts = 0; % facing away from both spheres

    elseif (du1>=0 && du2<=0)
        % if one dot product is positive, then probably post-tumble. go to
        % state 1

        % pointing into obstacle 1
        vt = ut - n1*du1;

        % does its velocity point away from obstacle 2?
        % if so, state 1, and in contact with obstacle 1
        if dot(vt,n2)<0
            contacts=1;
            circInds=nearbyCircInds(1);
        end

    elseif (du1<=0 && du2>=0)
        % if one dot product is positive, then probably post-tumble. go to
        % state 1

        % pointing into obstacle 2
        vt = ut - n2*du2;

        % does its velocity point away from obstacle 1?
        if dot(vt,n1)<0
            contacts=1;
            circInds=nearbyCircInds(2);
        end
        % if one dot product is positive, then probably post-tumble. go to
        % state 1

    elseif n1xut*n1xn2>0 && n2xut*n2xn1>0
        % but even if both dot products are positive, can still go to state
        % 1!

        % if heading lies between the two normal
        % ok, state 2
        

    elseif n1xut*n1xn2>0
        % finally, this is the 1->1 transition. need to determine which way
        % the cell goes

        % if only n1xut*n1xn2>0
        % implies ut and n2 are on the same side of n1...

        % in this case, cell goes to circle 2
        contacts=1;
        circInds=nearbyCircInds(2);

    elseif n2xut*n2xn1>0
        contacts=1;
        circInds=nearbyCircInds(1);
    end

end

if nargout>=4
    varargout{1} = phi_t;
end

% if contacts>0
% figure;hold on
% viscircles(circ_pos_i',ones(size(circ_pos_i,2),1),'color','k');
% %if contacts>0
%     viscircles(circ_pos_i(:,circInds)',ones(length(circInds),1),'color','r');
% %end
% quiver(xt(1),xt(2),ut(1),ut(2))
% contacts
% keyboard
% end

end



