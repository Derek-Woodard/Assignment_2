%% Assignment 2 Finite Difference Method
% Part 2a

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

% There are two sigma values for different regions
sig1 = 1;
sig2 = 10^-2;
sigs = zeros(nx,ny);

% We need to set the 'bottleneck' by creating two boxes in the region
boxes = [nx*2/5 nx*3/5 ny*2/5 ny*3/5];

% We need to set the signma values for the full region, giving different
% sigma values to the boxes
for i = 1:nx
    for j = 1:ny
        if i > boxes(1) && i < boxes(2) && (j < boxes(3) || j > boxes(4))
            sigs(i,j) = sig2;
        else
            sigs(i,j) = sig1;
        end
    end
end

%%
% The G matrix is formed using the number of points calculated above
G = zeros(nx*ny);
% Gb = zeros(nx*ny,nx*ny);

% Now we generate the F matrix
F = zeros(1,nx*ny);
% Fb = zeros(nx*ny,1);

for i=1:nx
    for j=1:ny
        % set up the mapping variables
        n = j+(i-1)*ny;
        nxp = j+i*ny;
        nxm = j+(i-2)*ny;
        nyp = j+1+(i-1)*ny;
        nym = j-1+(i-1)*ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;
        elseif j == 1
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp));
        elseif j == ny
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nym));
        else
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
            G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp) + G(n,nym));
        end
    end
end

% plot the sigma for the full region
figure(1)
surf(sigs);
xlabel('x');
ylabel('y');
zlabel('sigma');
title('Plot of sigma values in region');
axis tight
set(gca,'View', [45 45]);

% Now that the matricies are set up, the solution can be found
V = G\F';
% Vb = Gb\Fb;
S = zeros(ny,nx,1);
% Sb = zeros(nx,ny,1);

for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        S(j,i) = V(n);
    end
end

figure(2);
surf(S);
axis tight
xlabel('x');
ylabel('y');
title('Voltage with bottleneck conditions');
set(gca,'View', [45 45]);

% we plot the electric field by taking the gradient of the voltages
[ex ey] = gradient(S);

figure(3)
subplot(2,1,1);
surf(-ex)
axis tight
xlabel('x');
ylabel('y');
zlabel('Electric Field');
title('X-Componenet of the Electric Field');
set(gca,'View', [45 45]);

subplot(2,1,2);
surf(-ey)
axis tight
xlabel('x');
ylabel('y');
zlabel('Electric Field');
title('Y-Componenet of the Electric Field');
set(gca,'View', [45 45]);

% The current density is plotted using the sigma values along with the
% electric fields
Jx = sigs'.*ex;
Jy = sigs'.*ey;
J = sqrt(Jx.^2 + Jy.^2);
figure(4)
surf(J)
axis tight
xlabel('x');
ylabel('y');
zlabel('Current Density');
title('Plot of Current Density for Region');
set(gca,'View', [45 45]);
