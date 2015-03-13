# Transport-Project
That project we have to do for transport
# Max

width = 6;
length = 4;
extra = 3;
n = width*length-extra;
a = zeros(n);
c = zeros(n,1);
k = 25;
delta_x = .001;
h0 = 1000;
T0 = 1700;
h1 = 200;
T1 = 400;

for i = 1:width*length-extra
    a(i,i) = -4;
    
    if(mod(i,width) == 1) % on left edge
        a(i,i+1) = 2;
    elseif(mod(i,width) == 0) % on right edge
        a(i,i-1) = 2;
    elseif(i==n) % last
        a(i,i-1) = 2;
        c(i) = -2*h1*delta_x/k*T1;
        a(i,i) = -2*(h1*delta_x/k+2);
    else
        a(i,i+1) = 1;
        a(i,i-1) = 1;
    end
    
    if(i<=width) %% on top
        a(i,i+width) = 2;
        c(i) = -2*h0*delta_x/k*T0;
        a(i,i) = -2*(h0*delta_x/k+2);
    elseif(i>=n-width) %% on bottom
        a(i,i-width) = 2;
        if(i<=n-extra) % on boundary
            c(i) = -2*h1*delta_x/k*T1;
            a(i,i) = -2*(h1*delta_x/k+2);
        end
    else
        a(i,i+width) = 1;
        a(i,i-width) = 1;
    end
    
    if(i==n-width) %corner
        a(i,i-1) = 2;
        a(i,i-width) = 2;
        a(i,i+1) = 1;
        a(i,i+width) = 1;
        c(i) = -2*h1*delta_x/k*T1;
        a(i,i) = -2*(h1*delta_x/k+3);
    end
end

t = a\c;

disp = zeros(length,width);
for w = 1:width
    for l = 1:length
        if((l-1)*width+w<=n)
            disp(l,w) = t((l-1)*width+w,1);
        else
            disp(l,w) = NaN;
        end
    end
end

box = zeros(length-1,width-1,4);
for w = 1:width-1
    for l = 1:length-1
        box(l,w,1) = w+0.5;
        box(l,w,2) = l+0.5;
        box(l,w,3) = disp(l,w)-disp(l,w+1) - 1/sqrt(2)*(disp(l+1,w+1)-disp(l,w));
        box(l,w,4) = disp(l,w)-disp(l+1,w) - 1/sqrt(2)*(disp(l+1,w+1)-disp(l,w));
    end
end

hold on
pcolor(disp); shading interp;
quiver(box(:,:,1),box(:,:,2),box(:,:,3),box(:,:,4));
