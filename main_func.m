% Transport Project
% 23 March 2015
% Max Payson, Lori Kaufman, Sam Faucher, Molly Wolf

% This takes Sam's code and defines a matrix "loc" that categorizes the
% location of each node
% loc=0 interior; loc=1 top left; loc=2 top; loc=3 top right; loc=4 right; 
% loc=5 bottom right; loc=6 bottom; loc=7 bottom left; loc=8 left; 
function[loc] = set_location(w, l, n)
    loc = zeros(n);
    loc(1)=1;  % Top left point
    loc(n)=5;  % Bottom right point
    loc(width)=3;  % Top right point
    loc(n-width+1)=7;  % Bottom left point
    loc(2:width-1)=2; % Top points
    loc(n-width+2:n-1)=6; % Bottom points
    loc(width+1:width:n-2*width+1)=8; % Left points
    loc(2*width:width:n-width)=4; % Right points
end

% This creates a new row to be added to the temperature coefficient matrix corresponding
% to an interior node. The values are from Sam's code and the given
% equation for an interior node. The constant vector is defined in main.
function[new_row] = interior_step(w, l, i, Fo)
    new_row = zeros(w*l);
    new_row(i) = 1+4*Fo;
    new_row(i-1) = -Fo;
    new_row(i+1) = -Fo;
    new_row(i-w) = -Fo;
    new_row(i+w) = -Fo;
end

% This creates a new row to be added to the temperature coefficient matrix corresponding
% to an edge node. Here, the given equation is adapted such that the "node" that exists outside
% the plate is not included and is rather represented in the constant vector defined in main.
% Note pos is location (loc(i)) for the current node
function[new_row] = edge_step(w,l,i, Fo, pos)

    % Initialize new row
    new_row = zeros(w*l);

    % The given node "i" always has the same coefficient
    new_row(i) = 1+4*Fo;

    % If the top edge, don't include the non-existent top node
    if pos == 2
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i+w) = -Fo
    % If right edge, don't include the non-existent right node
    elseif pos == 4
        new_row(i-1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo
    % If the bottom edge, don't include the non-existent bottom node
    elseif pos == 6
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo
    % If the left edge, don't include the non-existent left node
    elseif pos == 8
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo
    end
end

% This creates a new row to be added to the temperature coefficient matrix corresponding
% to a corner node. Here, the given equation is adapted such that the "nodes" that exist outside
% the plate are not included and are rather represented in the constant vector defined in main.
% Note pos is location (loc(i)) for the current node
function[new_row] = edge_step(w,l,i,Fo,pos)
    new_row = zeros(w*l);
    new_row(i) = 1+4*Fo;
    if pos == 1
        new_row(i+1) = -Fo;
        new_row(i+w) = -Fo;
    elseif pos == 3
        new_row(n-1) = -Fo;
        new_row(n+w) = -Fo;
    elseif pos == 5
        new_row(n-1) = -Fo;
        new_row(n-w) = -Fo;
    elseif pos == 7
        new_row(n+1) = -Fo;
        new_row(n-w) = -Fo;
    end
end

% This defines the relevant constant for the constant vector corresponding 
% to a corner node. The constant has three options, either heat flow from
% both sides, constant temperature from both sides, or both
% Note temperature equals the previous temperature at the node (ie T(j,i))
% And corner represents the location (loc(i))
function[const] = corner_stepConst(corner, Fo, q, temp)
    % Determine whether there is heat flow at each edge for the corner
    if corner == 1
        q1 = q(1);
        q2 = q(4);
    elseif corner == 3
        q1 = q(1);
        q2 = q(2);
    elseif corner == 5
        q1 = q(2);
        q2 = q(3);
    elseif corner == 7
        q1 = q(3);
        q2 = q(4);
    end

    % If constant temperature on both sides
    if q1 == 0 & q2 == 0
        const = temp + 2*Fo;
    % Heat flow on both sides
    elseif q1 ~= 0 & q2 ~= 0
        const = temp + 2*.005*temp;
    % heat flow one side constant temperature other side
    else
        const = temp + Fo + .005*temp;
    end
end


function [] = main_func()
    % User can choose width, length, steps, and Fo (dimensionless time) per step.
    l = 20;
    w = 40;
    steps=100;
    Fo=0.25;

    % This defines the q boundary conditions, where a value of 0 indicates no heat
    % flow and a value other than 0 indicates heat flow
    % q(1) = top, q(2) = right, q(3) = bottom, q(4) = left
    q = [0,0,0,0]

    n=w*l;
    T=zeros(steps,n);  % Initialize temperature matrix
    % T(1,:)=0.5  % If you wanted to specify other initial temp, you could here.

    % Initialize matrix that specifies category of point.
    loc = set_location(w, l, n);    
    %  Categorize points as top, bottom, left, right, top left, top right,
    %  bottom left, bottom right

    % Initialize MAT, the coefficient matrix, and C, the constnat matrix, similar to past pset.
    MAT=zeros(n,n);
    C=zeros(n,1);  

    for j=1:steps-1  % For each temperature step
        for i=1:n  % For each point

            %interior node
            if loc(i)==0   
                MAT(i,:) = interior_step(w,l, i, Fo); % Add coefficients to temp coeff matrix
                C(i,1)=T(j,i);  % The constant equals the previous temperature value at that point.

            %edge node since edges are all even locations
            elseif mod(loc(i), 2) == 0
                MAT(i,:) = edge_step(w, l, i, Fo, loc(i)); % Add coefficients to temp coeff matrix

                % no heat flow for relevant edge, the location divided by 2 gives
                % the correct index for the boundary condition vector q
                % the boundary conditions are substituted for neighboring T values
                % note T will always equal 1 if we assume T is the dimensionless T
                if q(loc(i)/2) == 0
                    C(i,1) = T(j,i) + Fo;

                % heat flow for relevant edge
                % If there's heat flow, the heat flow is s.t. it increases normalized temperature
                % by .005 units per time step, thus q*Fo = .005*T(n) 
                else
                    C(i,1) = T(j,i) + .005*T(j,i);
                end
            
            %corner node
            else
                MAT(i,:) = corner_step(w, l, i, Fo, loc); % Add coefficients to temp coefficient matrix
                C(i,1) = corner_stepConst(loc(i), Fo, q, T(j,i)); % Define relevant constant in constant vector
            end
        end
        temp=inv(MAT)*C;  % To calculate the temp at the next time step, multiply inverse of A matrix by C.
        T(j+1,:)=temp(:);  % Store it in the next row of the T matrix.
    end

    % Now we have a T matrix that shows the progression of temperature.


    % Reformat 2D T matrix to 3D T matrix, displayT.

    displayT=zeros(length,width,steps);
    for k=1:steps
        for l=1:length
            displayT(l,:,k)=T(k,(l-1)*width+1:l*width);
        end
    end


    % Plot displayT over time.

    for m=1:steps
        imagesc(displayT(:,:,m));
        tau=(m-1)*Fo;
        tau=num2str(tau);
        text=['Fo = ' tau];
        title(text,'fontsize',24);
        pause(.1)
    end
end
