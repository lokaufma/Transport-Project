% Transport Project
% 23 March 2015
% Max Payson, Lori Kaufman, Sam Faucher, Molly Wolf

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

function[new_row] = interior_step(w, l, i, Fo)
    new_row = zeros(w*l);
    new_row(i) = 1+4*Fo;
    new_row(i-1) = -Fo;
    new_row(i+1) = -Fo;
    new_row(i-w) = -Fo;
    new_row(i+w) = -Fo;
end

function[new_row] = edge_step(w,l,i, Fo, loc)
    new_row = zeros(w*l);
    new_row(i) = 1+4*Fo;
    if loc(i) == 2
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i+w) = -Fo

    elseif loc(i) == 4
        new_row(i-1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo

    elseif loc(i) == 6
        new_row(i-1) = -Fo
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo

    elseif loc(i) == 8
        new_row(i+1) = -Fo
        new_row(i-w) = -Fo
        new_row(i+w) = -Fo
    end
end

function[new_row] = edge_step(w,l,i,Fo,loc)
    new_row = zeros(w*l);
    new_row(i) = 1+4*Fo;
    if loc(i) == 1
        new_row(i+1) = -Fo;
        new_row(i+w) = -Fo;
    elseif loc(i) == 3
        new_row(n-1) = -Fo;
        new_row(n+w) = -Fo;
    elseif loc(i) == 5
        new_row(n-1) = -Fo;
        new_row(n-w) = -Fo;
    elseif loc(i) == 7
        new_row(n+1) = -Fo;
        new_row(n-w) = -Fo;
    end
end

function[const] = corner_stepConst(corner, Fo, q, temp)
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

    if q1 == 0 & q2 == 0
        const = temp + 2*Fo;
    elseif q1 ~= 0 & q2 ~= 0
        const = temp + 2*.005*temp;
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

    q = [0,0,0,0] %if = 0, no heat flow and bc is T = 1

    n=w*l;
    T=zeros(steps,n);  % Initialize temperature matrix
    % T(1,:)=0.5  % If you wanted to specify other initial temp, you could here.

    loc = set_location(w, l, n);  % Initialize matrix that specifies category of point.  
    % cat=0 interior; cat=1 top left; cat=2 top; cat=3 top right; cat=4 right; 
    % cat=5 bottom right; cat=6 bottom; cat=7 bottom left; cat=8 left; 
    %  Categorize points as top, bottom, left, right, top left, top right,
    %  bottom left, bottom right

    % Initialize A and C matrices, similar to past pset.

    MAT=zeros(n,n);
    C=zeros(n,1);  

    for j=1:steps-1  % For each temperature step
        for i=1:n  % For each point

            %interior node
            if loc(i)==0   
                MAT(i,:) = interior_step(w,l, i, Fo);
                C(i,1)=T(j,i);  % The constant equals the previous temperature value at that point.

            %edge node
            elseif mod(loc(i), 2) == 0
                MAT(i,:) = edge_step(w, l, i, Fo, loc);

                %no heat flow for relevant edge
                if q(loc(i)/2) == 0
                    C(i,1) = T(j,i) + Fo;
                %heat flow for relevant edge
                else
                    C(i,1) = T(j,i) + .005*T(j,i);
                end
            
            else
                MAT(i,:) = corner_step(w, l, i, Fo, loc);
                C(i,1) = corner_stepConst(loc(i), Fo, q, T(j,i));
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
