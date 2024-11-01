function new_circ_pos = placeSpheres(square_halfL,rho,new_square_pos,past_squares_pos)
% new_circ_pos = placeSpheres(square_halfL,rho,new_square_pos,past_squares_pos)
%
% Add new circular obstacles to the simulation domain.
%
% Inputs:
% square_halfL: Half the side length of each square, in units of the obstacle radius. scalar.
% rho: Dimensionless number density of obstacles (per area). scalar.
% new_square_pos: Position of the new square. d by 1 vector., where d is the number of spatial dimensiones
% past_squares_pos: Positions of the centers of previously-placed squares. d by Nsquare matrix.
%
% Outputs:
% new_circ_pos: Positions of the centers of newly-placed circular obstacles. These centers lie inside the square centered at new_square_pos, and not inside any previously placed squares. d by Ncircles matrix.
%
% Henry H. Mattingly, November 2023

[d, ~] = size(new_square_pos);
[~, Npast] = size(past_squares_pos);

% mean number of spheres in the new cube
Acube = (2*square_halfL)^d;
mean_Ncirc = rho*Acube;

% generate new spheres
Ncirc_new = poissrnd(mean_Ncirc);
new_circ_pos = 2*square_halfL*(rand(d,Ncirc_new)-1/2)+new_square_pos;

% check if any lie inside other squares
if ~isempty(past_squares_pos)
    % slow
%     tic
%     rmInds = any(all(abs(repmat(new_circ_pos,[1,1,Npast])-repmat(permute(past_squares_pos,[1,3,2]),[1,Ncirc_new,1]))<square_halfL,1),3);
%     toc

    
%     tic
    rmInds = getRmInds(new_square_pos,past_squares_pos,new_circ_pos,Npast,Ncirc_new,square_halfL);
%     toc

%     tic
%     rmInds = false(1,Ncirc_new);
%     for k = 1:Ncirc_new
%         temp = repmat(new_circ_pos(:,k),1,Npast) - past_squares_pos;
%         rmInds(k) = any(all(abs(temp)<square_halfL,1),2);
%     end
%     toc

    new_circ_pos(:,rmInds) = [];
end

% this could be improved by identifying the minimal set of squares that
% fully cover the space explored

end


function rmInds = getRmInds(new_square_pos,past_squares_pos,new_circ_pos,Npast,Ncirc_new,square_halfL)

    % find close potential squares only
    D = sqrt(sum((repmat(new_square_pos,1,Npast) - past_squares_pos).^2,1));
    close_squares = find(D<2.1*sqrt(2)*square_halfL); %2?
    nclose = length(close_squares);

    %
    rmInds = any(all(abs(repmat(new_circ_pos,[1,1,nclose])-repmat(permute(past_squares_pos(:,close_squares),[1,3,2]),[1,Ncirc_new,1]))<square_halfL,1),3);

%     figure;hold on
%     viscircles(new_circ_pos',1,'color','k');
%     for i = 1:size(past_squares_pos,2)
%         rectangle('Position',[past_squares_pos(:,i)'-square_halfL*ones(1,2),2*square_halfL*ones(1,2)]);
%     end
%     for i = 1:length(close_squares)
%         rectangle('Position',[past_squares_pos(:,close_squares(i))'-square_halfL*ones(1,2),2*square_halfL*ones(1,2)]);
%     end
% 
%     viscircles(new_circ_pos(:,rmInds)',1,'color','r');

%     if any(rmInds==0)
%         keyboard
%     end
end