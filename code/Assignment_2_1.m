%% Assignment 2 Finite Difference Method
% Part 1

clear all;
close all;
set(0,'DefaultFigureWindowStyle','docked');

%%
% set the size of the region being used

L = 3;
W = 2;

% Set Vo arbitrarily to 1
Vo = 1;

% Set the spacing for the mesh and the number of nodes total based on
% region size
dx = 0.1;
dy = 0.1;
nx = L/dx;
ny = W/dy;

%%
% Using the finite difference method, a matrix G is used where GV = F
% V represetnts the voltage at discretized points within thje selescted
% region.
% F represents a matrix that forces the boundary conditions.
% G represents the voltage at each point relative to the points surrounding
% it.
% 


% The G matrix is formed using the number of points calculated above
G = zeros(nx*ny,nx*ny);
Gb = zeros(nx*ny,nx*ny);

% Now we generate the F matrix
F = zeros(nx*ny,1);
Fb = zeros(nx*ny,1);

for i=1:nx
    for j=1:ny
        % set up the mapping variable for G to F
        n = j+(i-1)*ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
            
            Gb(n,:) = 0;
            Gb(n,n) = 1;
            Fb(n) = 0;
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;
            
            Gb(n,:) = 0;
            Gb(n,n) = 1;
            Fb(n) = 0;
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,n+1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;
            
            Gb(n,:) = 0;
            Gb(n,n) = 1;
            Fb(n) = 1;
        elseif j == ny
            G(n,n) = -3;
            G(n,n-1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;
            
            Gb(n,:) = 0;
            Gb(n,n) = 1;
            Fb(n) = 1;
        else
            G(n,n) = -4;
            G(n,n+1) = 1;
            G(n,n-1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;
            
            Gb(n,n) = -4;
            Gb(n,n+1) = 1;
            Gb(n,n-1) = 1;
            Gb(n,n+ny) = 1;
            Gb(n,n-ny) = 1;
        end
    end
end


% Now that the matricies are set up, the solution can be found
V = G\F;
Vb = Gb\Fb;
S = zeros(nx,ny,1);
Sb = zeros(nx,ny,1);

for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        S(i,j) = V(n);
        Sb(i,j) = Vb(n);
    end
end

figure(1);
surf(S');
hold on;
xlabel('x');
ylabel('y');
title('1a) F.D.M. Solution, Left Border V = Vo');
set(gca,'View', [45 45]);

figure(2);
subplot(2,1,1);
surf(Sb);
xlabel('x');
ylabel('y');
title('1b) F.D.M. Solution, Left+Right Border V = Vo');
set(gca,'View', [45 45]);

% The analytical solution for part b must be done to compare to the numeric
a = ny;
b = nx/2;

x2 = linspace(-nx/2,nx/2,30);
y2 = linspace(0,ny,ny);

[i,j] = meshgrid(x2,y2);
Sb2 = zeros(ny,nx);

for n = 1:2:500
    Sb2 = (Sb2+(cosh(n*pi*i/a).*sin(n*pi*j/a))./(n*cosh(n*pi*b/a)));
    subplot(2,1,2);
    surf(x2,y2,(4/pi)*Sb2)
    axis tight
    title('Analytical Solution, Left+Right Border V = Vo');
    xlabel('x');
    ylabel('y');
    set(gca,'View', [45 45]);
    pause(0.01)
end


