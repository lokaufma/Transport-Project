% Transport Project
% 23 March 2015
% Max Payson, Lori Kaufman, Sam Faucher, Molly Wolf

% Function: set_location()
% This function creates a matrix loc that indicates a category for each point
% in the matrix that indicates which edges that point borders.
% Inputs: width (w), length (l), and number of points (n) in matrix
% Outputs: location matrix (loc) defining a type for each matrix point:
% 	loc=0 interior; loc=1 top-left; loc=2 top; loc=3 top-right; loc=4 right 
% 	loc=5 bottom-right; loc=6 bottom; loc=7 bottom-left; loc=8 left.

function[loc] = set_location(w, l, n)
    loc = zeros(n);
    loc(1) = 1; % Top-left point
    loc(n) = 5; % Bottom-right point
    loc(width) = 3; % Top-right point
    loc(n-width+1) = 7; % Bottom-left point
    loc(2:width-1) = 2; % Top points
    loc(n-width+2:n-1) = 6; % Bottom points
    loc(width+1:width:n-2*width+1) = 8; % Left points
    loc(2*width:width:n-width) = 4; % Right points
end

% Function: interior_step()
% This function creates a row to be added to the temperature coefficient matrix 
% corresponding to an interior node.
% Inputs: width (w), length (l), location in matrix (i), Fourier number (Fo)
% Outputs: row to be added to temperature coefficient matrix (new_row)

function[new_row] = interior_step(w, l, i, Fo)
    new_row = zeros(w*l);
    new_row(i) = 1 + 4*Fo;
    new_row(i-1) = -Fo;
    new_row(i+1) = -Fo;
    new_row(i-w) = -Fo;
    new_row(i+w) = -Fo;
end

% Function: edge_step()
% This function creates a row to be added to the temperature coefficient matrix
% corresponding to an edge node.
% Inputs: width (w), length (l), location in matrix (i), Fourier number (Fo),
% location matrix value representing the edges bordered (pos)
% Outputs: row to be added to temperature coefficient matrix (new_row)

function[new_row] = edge_step(w, l, i, Fo, pos)
    new_row = zeros(w*l); % Initialize new row
    new_row(i) = 1+4*Fo; % Coefficient for given node i

    if pos == 2 % If top edge, don't include non-existent top node
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i+w) = -Fo
    elseif pos == 4 % If right edge, don't include non-existent right node
        new_row(i-1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo
    elseif pos == 6 % If bottom edge, don't include non-existent bottom node
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo
    elseif pos == 8 % If left edge, don't include non-existent left node
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo
    end
end

% Function: corner_step()
% This function creates a row to be added to the temperature coefficient matrix
% corresponding to a corner node.
% Inputs: width (w), length (l), location in matrix (i), Fourier number (Fo),
% location matrix value representing the edges bordered (pos)
% Outputs: row to be added to temperature coefficient matrix (new_row)

function[new_row] = corner_step(w, l, i, Fo, pos)
    new_row = zeros(w*l);
    new_row(i) = 1 + 4*Fo;
    
    if pos == 1 % Top-left
        new_row(i+1) = -Fo;
        new_row(i+w) = -Fo;
    elseif pos == 3 % Top-right
        new_row(i-1) = -Fo;
        new_row(i+w) = -Fo;
    elseif pos == 5 % Bottom-right
        new_row(i-1) = -Fo;
        new_row(i-w) = -Fo;
    elseif pos == 7 % Bottom-left
        new_row(i+1) = -Fo;
        new_row(i-w) = -Fo;
    end
end

% Function: corner_stepConst()
% This function finds the relevant constant for the constant vector C
% corresponding to a corner node. The three possible conditions are: 
% 	(1) heat flow from both sides, (2) constant temperature at both sides,
% 	and (3) heat flow from one side and constant temperature at the other.
% Inputs: location matrix value representing the edges bordered (corner),
% Fourier number (Fo), heat flow matrix (q), current temperature matrix (temp)
% Outputs: constant matrix (const)

function[const] = corner_stepConst(corner, Fo, q, temp)
    % Determine whether there is heat flow at each edge for the corner
    if corner == 1 % Top-left
        q1 = q(1); % Top
        q2 = q(4); % Left
    elseif corner == 3 % Top-right
        q1 = q(1); % Top
        q2 = q(2); % Right
    elseif corner == 5 % Bottom-right
        q1 = q(2); % Right
        q2 = q(3); % Bottom
    elseif corner == 7 % Bottom-left
        q1 = q(3); % Bottom
        q2 = q(4); % Left
    end

    if q1 == 0 & q2 == 0 % (1) Constant temperature on both sides
        const = temp + 2*Fo;
    elseif q1 ~= 0 & q2 ~= 0 % (2) Heat flow on both sides
        const = temp + 2*.05*temp;
    else % (3) Heat flow on one side and constant temperature on the other
        const = temp + Fo + .05*temp;
    end
end

function [] = main_func()
    % Define width, length, number of time steps, and 
    % Fourier number Fo (dimensionless time) per step
    l = 20;
    w = 40;
    steps = 100;
    Fo = 0.25;

    % Define the heat flow boundary conditions (q=0 represents no heat flow)
    % q(1) = top, q(2) = right, q(3) = bottom, q(4) = left
    q = [0,0,0,0]

    n=w*l; % Total number of matrix points
    T = zeros(steps,n);  % Initialize temperature matrix
    % T(1,:) = 0.5  % Specify initial conditions for temperature matrix

    % Initialize loc matrix specifying category of point
    loc = set_location(w, l, n);    

    % Initialize MAT, the coefficient matrix, and C, the constant matrix
    MAT = zeros(n,n);
    C = zeros(n,1);  

    for j = 1:steps-1  % Step through time
        for i = 1:n  % Step through the temperature matrix
            if loc(i)==0  % Interior node
                MAT(i,:) = interior_step(w,l, i, Fo); % Add row to MAT
                C(i,1) = T(j,i);
            elseif mod(loc(i), 2) == 0 % Edge node (even number in loc)
                MAT(i,:) = edge_step(w, l, i, Fo, loc(i));
                if q(loc(i)/2) == 0 % No heat flow through this edge
                    C(i,1) = T(j,i) + Fo;
                else % There is heat flow through this edge
                    C(i,1) = T(j,i) + .05*T(j,i); % Increase .05 per time step 
                end
            else % Corner node
                MAT(i,:) = corner_step(w, l, i, Fo, loc);
                C(i,1) = corner_stepConst(loc(i), Fo, q, T(j,i));
            end
        end

        temp = inv(MAT)*C;  % Calculate temperature at next time step
        T(j+1,:) = temp(:);  % Store calculated temperature in T matrix
    end

    % Reformat 2D T matrix to 3D temperature matrix displayT

    displayT = zeros(length, width, steps);
    for k=1:steps % Step through time
        for l=1:length % Step through rows
            displayT(l,:,k) = T(k,(l-1)*width+1:l*width);
        end
    end

    % Plot displayT over time

    for m = 1:steps
        imagesc(displayT(:,:,m));
        tau = (m-1)*Fo;
        tau = num2str(tau);
        text = ['Fo = ' tau];
        title(text,'fontsize',24);
        pause(.1)
    end
end
