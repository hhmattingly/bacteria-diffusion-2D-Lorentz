function visualize_traj(cell_num,tt,circ_pos,squares_pos,square_halfL,xt,varargin)
% visualize_traj(cell_num,tt,circles,squares,square_halfL,xt,tumbles,contacts,vt,options)
%
% Inputs:
% cell_num: Index of the cell to plot. scalar integer
% tt: Maximum time index to plot. scalar integer
% circ_pos: Matrix containing the positions of the centers of the circular obstacles. d by Nircles_initial matrix, for d dimensions
% squares_pos: Matrix containing the positions of the centers of the squares defining the simulation domain. d by 1 matrix, for d dimensions
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar
% xt: Matrix of all cell positions over time. Ncells by d by nt matrix, for d dimensions
% tumbles: Matrix of cell run/tumble state over time. 0 indicates run state, 1 indicates tumble state. Ncells by nt logical vector.
% contacts: Matrix of cell contact state over time. 0 indicates swimming in bulk, 1 indicates contact with one obstacle, 2 indicates contact with two obstacles. Ncells by nt vector.
% vt: Matrix of cell velocities over time. Ncells by d by nt matrix, for d dimensions. As is, this variable is not used
%
% Optional inputs: 
% options: Structure of plotting options.
%   Fields:
%   showRectangles: Whether or not to plot the squares that form the simulation domain. logical scalar. default value is true
%   showGrid: Whether or not to show major grid lines. logical scalar. default value is true
%   showDots: Whether or not to show dots at the points where the cell's position is recorded.
%
% Outputs:
% None. Generates a figure.
%
% Henry H. Mattingly, November 2023

red = [215,48,39]/255;
gray=0.8*ones(1,3);
circ_rad = 1;

showTumbles=0;
showContacts=0;
if nargin>=7
    tumbles=varargin{1};
end
if nargin>=8
    contacts=varargin{2};
end
if nargin>=9
    vt=varargin{3};
end
if nargin>=10
    options=varargin{4};
end


if ~isempty(options)
    showRectangles=options.showRectangles;
    showGrid=options.showGrid;
    showDots=options.showDots;
    if isfield(options,'showTumbles')
        showTumbles=options.showTumbles;
    else
        showTumbles=1;
    end
    if isfield(options,'showContacts')
        showContacts=options.showContacts;
    else
        showContacts=1;
    end
    if isfield(options,'quickPlot')
        quickPlot=options.quickPlot;
    else
        quickPlot=0;
    end
else
    showRectangles=1;
    showGrid=1;
    showDots=1;
%     showTumbles=1;
%     showContacts=1;
    quickPlot=0;
end

circ_pos = circ_pos{cell_num};
squares_pos = squares_pos{cell_num};
Nsquares = size(squares_pos,2);
Ncirc = size(circ_pos,2);

figure;hold on


if Ncirc>0 && showRectangles
    for i = 1:Nsquares
        h = rectangle('Position',[squares_pos(:,i)'-square_halfL*ones(1,2),2*square_halfL*ones(1,2)]);
    end
end


if Ncirc>0
    if quickPlot
        h=viscircles(circ_pos',circ_rad*ones(size(circ_pos,2),1),'color',gray);
    else
        for i = 1:size(circ_pos,2)
            h = rectangle('Position',[circ_pos(1,i)-1,circ_pos(2,i)-1,2,2],'Curvature',[1,1]);
            h.LineWidth=1;
            h.EdgeColor = gray;
            h.FaceColor = gray;
        end
    end
end

x1 = squeeze(xt(cell_num,:,1:tt));

if showDots
    plot(x1(1,:),x1(2,:),'.-','LineWidth',1)
else
    plot(x1(1,:),x1(2,:),'LineWidth',2)
end

if Ncirc>0

    if showContacts
        contacts_tf_i = contacts(cell_num,1:tt)>0;
        ind1=0;
        ind2=0;

        while ind2<tt
            ind1=find(contacts_tf_i(ind2+1:end)==1,1,'first');
            ind1 = ind1+ind2;
            ind2=find(contacts_tf_i(ind1+1:end)==0,1,'first');
            if isempty(ind2)
                ind2=tt;
            else
                ind2 = ind2+ind1-1;
            end

            plot(x1(1,ind1:ind2),x1(2,ind1:ind2),'k','LineWidth',2)
        end
    end
end

if showTumbles
    tumblesi = logical(tumbles(cell_num,1:tt));
    plot(x1(1,tumblesi),x1(2,tumblesi),'o','Color',red,'MarkerFaceColor',red,'LineWidth',1,'MarkerSize',2) %3
end

axis equal

xlabel('x')
ylabel('y')
zlabel('z')

if showGrid
    h=gca;
    h.XGrid='on';
    h.YGrid='on';
end
end