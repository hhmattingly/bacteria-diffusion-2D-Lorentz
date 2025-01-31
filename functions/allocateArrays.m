function [xt,tht,vt,tumbles,contacts,squares,circles,circInds,closed,nbreak,nsims,RNs_tum] = allocateArrays(Ncells,d,nt,continue_flag,varargin)
% xt stores cell positions
% tht stores cell heading angles
% vt stores cell velocities
% tumbles stores when tumble occur
% contacts stores the cell's surface state
% squares stores the locations of squares in the domain that the cell has
% explored
% circles stores the obstacle locations
% circInds stores which obstacles the cell is in contact with, if any


if nargin>=5
    file_name = varargin{1};
end

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

        % was environment closed at the beginning/end of the
        % simulation? really becomes an indicator of whether
        % the simulation was run or not
        closed = ones(Ncells,1);
    end

end

% random numbers for tumbles
RNs_tum = rand(Ncells,nt);

end