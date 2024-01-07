function P = Path_generator_esc(Xs, Xf, eta, B, alpha,eps,p_o)
P = Xs; %Initializing P with the start point Xs
current_pos = Xs; % Initializng the current position with the start point Xs
obstacle_indices = unique(B(3,:)); % Indices(Number) of the obstacles given in B
obstacles = []; % Initializing the obstacles struct
for i = obstacle_indices
    % Vertices of the obstacle
    obstacle_points = B(1:2, B(3,:) == i);
    % Saving the vertices in "points" field of the "obstacle" struct
    obstacles{i}.points = obstacle_points;
    % Saving the number of the vertices or actually the number of the sides of the obstacle
    obstacles{i}.numPoints = size(obstacle_points, 2);
end
while norm(current_pos - Xf) >= 0.1
    if size(P,2)>2 % Local minima check
        % If in 2 iteration it goes back to the initial location then
        % it's in local minima
        if current_pos==P(:,end-2)
            % Calculating the potential of the point
            V_current=potential(current_pos, Xf, eta, obstacles, alpha, p_o);
            V_local=V_current;
            % Weighting the repulsive force more than attractive force
            % in order to escape the local minimum
            F_net_local= F_rep_total + 0.05*F_attr;
            pos_local=current_pos;
            eps1=eps;
            P=P(:,1:end-2);
            % If the new point's potential is at least 25% more than the
            % potential of the local minimum, I assume that it has escaped
            % the local minimum.
            while 1.25*V_local>V_current
                % Instead of the testing multiple direction, by weighting
                % the repulsive force more than the attractive force, 
                % I only use that direction for moving.
                pos_local=pos_local+F_net_local/norm(F_net_local)*eps1;
                % If the step size isn't sufficient enough to escape the 
                % local minimum, it's increased in every iteration.
                eps1=eps1*1.25;
                P=[P pos_local];
                V_current=potential(pos_local, Xf, eta, obstacles, alpha, p_o);
            end
            current_pos=P(:,end);
        end
    end
    F_attr = -eta * (current_pos - Xf); % Computing the attractive force.
    F_rep_total = [0; 0]; % Initializing the total repulsive force.
    for i = obstacle_indices
        % Initializing the repulsive force for the i-th obstacle.
        F_rep = [0; 0];
        % Finding the nearest side of the obstacle
        middle=mean(obstacles{i}.points,2);
        [~,~,j]=polyxpoly([obstacles{i}.points(1,:) obstacles{i}.points(1,1)], [obstacles{i}.points(2,:) obstacles{i}.points(2,1)],[current_pos(1) middle(1)],[current_pos(2) middle(2)]);
        j=j(1);
        vertex1=[obstacles{i}.points(:,j);0];
        % Circular order of the vertices
        if j == obstacles{i}.numPoints
            vertex2=[obstacles{i}.points(:, 1);0];
        else
            vertex2=[obstacles{i}.points(:, j+1);0];
        end
        % Finding the shortest distance to the obstacle
        a=vertex1-vertex2;
        b=[current_pos;0]-vertex2;
        if 0<=dot(a,b) && dot(a,b)<=dot(a,a)
            rep_vector=-a*dot(a,b)/(norm(a)^2)+b;
            rep_dist=norm(cross(a,b))/norm(a);
        else
            % In this case the projection point is out of
            % line segment. Thus the shortest distance to
            % the side is from one of the vertices of the side.
            [rep_dist,k]=min([norm(b) norm(b-a)]);
            if k==1
                rep_vector=b;
            else
                rep_vector=b-a;
            end
        end
        if rep_dist < p_o
            rep_force_magnitude = alpha*(1/rep_dist-1/p_o) / (rep_dist^3);
            F_rep = rep_force_magnitude * (rep_vector / rep_dist);
        end
        F_rep_total = F_rep_total + F_rep(1:2);
    end
    F_net = F_attr + F_rep_total;
    new_pos = current_pos + (F_net / norm(F_net)) * eps;
    P = [P new_pos];
    current_pos = new_pos;
end
P=[P Xf];
end
% Potential function
function V = potential(current_pos, Xf, eta, obstacles, alpha, p_o)
V_attr = 0.5*eta * (norm(current_pos - Xf)^2);
V_rep_total = 0;
obstacle_indices = size(obstacles,2);
for i = obstacle_indices
    V_rep=0;
    middle=mean(obstacles{i}.points,2);
    [~,~,j]=polyxpoly([obstacles{i}.points(1,:) obstacles{i}.points(1,1)], [obstacles{i}.points(2,:) obstacles{i}.points(2,1)],[current_pos(1) middle(1)],[current_pos(2) middle(2)]);
    j=j(1);
    vertex1=[obstacles{i}.points(:,j);0];
    if j == obstacles{i}.numPoints
        vertex2=[obstacles{i}.points(:, 1);0];
    else
        vertex2=[obstacles{i}.points(:, j+1);0];
    end
    a=vertex1-vertex2;
    b=[current_pos;0]-vertex2;
    if 0<=dot(a,b) && dot(a,b)<=dot(a,a)
        rep_dist=norm(cross(a,b))/norm(a);
    else
        rep_dist=min([norm(b) norm(b-a)]);
    end
    if rep_dist < p_o
        V_rep=0.5*alpha*(1/rep_dist-1/p_o)^2;
    end
    V_rep_total = V_rep_total + V_rep;
end
V = V_attr + V_rep;
end