%% Assignment 2 Finite Difference Method
% Part 2c

clear all;
close all;
set(0,'DefaultFigureWindowStyle','docked');


%% 
% Test out different bottleneck sizes
for bn = 0.1:0.01:0.9
    nx = 50;
    ny = (3/2)*nx;
    G = zeros(nx*ny);
    F = zeros(1,nx*ny);
    
    % There are two sigma values for different regions
    sigs = zeros(ny,nx);
    sig1 = 1;
    sig2 = 10^-2;
    
    
    % We need to set the 'bottleneck' by creating two boxes in the region
    boxes = [nx*2/5 nx*3/5 ny*bn ny*(1-bn)];
    
    for i = 1:nx
        for j = 1:ny
            n = j+(i-1)*ny;
            if i == 1
                G(n,:) = 0;
                G(n,n) = 1;
                F(n) = 1;
            elseif i == nx
                G(n,:) = 0;
                G(n,n) = 1;
                F(n) = 0;
            elseif j == 1
                if i > boxes(1) && i < boxes(2)
                    G(n,n) = -3;
                    G(n,n+1) = sig2;
                    G(n,n+ny) = sig2;
                    G(n,n-ny) = sig2;
                else
                    G(n,n) = -3;
                    G(n,n+1) = sig1;
                    G(n,n+ny) = sig1;
                    G(n,n-ny) = sig1;
                end
            elseif j == ny
                if i > boxes(1) && i < boxes(2)
                    G(n,n) = -3;
                    G(n,n+1) = sig2;
                    G(n,n+ny) = sig2;
                    G(n,n-ny) = sig2;
                else
                    G(n,n) = -3;
                    G(n,n+1) = sig1;
                    G(n,n+ny) = sig1;
                    G(n,n-ny) = sig1;
                end
            else
                if i > boxes(1) && i < boxes(2) && (j < boxes(3) || j > boxes(4))
                    G(n,n) = -4;
                    G(n,n+1) = sig2;
                    G(n,n-1) = sig2;
                    G(n,n+ny) = sig2;
                    G(n,n-ny) = sig2;
                else
                    G(n,n) = -4;
                    G(n,n+1) = sig1;
                    G(n,n-1) = sig1;
                    G(n,n+ny) = sig1;
                    G(n,n-ny) = sig1;
                end
            end
        end
    end
    
    % We need to set the signma values for the full region, giving different
    % sigma values to the boxes
    for L = 1:nx
        for W = 1:ny
            if L >= boxes(1) && L <= boxes(2)
                sigs(W,L) = sig2;
            else
                sigs(W,L) = sig1;
            end
            if L >= boxes(1) && L <= boxes(2) && W >= boxes(3) && W <= boxes(4)
                sigs(W,L) = sig1;
            end
        end
    end
    
    V = G\F';
    S = zeros(ny,nx,1);
    
    for i = 1:nx
        for j = 1:ny
            n = j+(i-1)*ny;
            S(j,i) = V(n);
        end
    end
    
    [ex ey] = gradient(S);
        
    Jx = sigs.*ex;
    Jy = sigs.*ey;
    J = sqrt(Jx.^2 + Jy.^2);
    
    figure(1)
    hold on
    
    if bn == 0.1
        cur = sum(J,2);
        tot_cur = sum(cur);
        old_cur = tot_cur;
        plot([bn,bn], [old_cur,tot_cur])
    end
    if bn > 0.1
        old_cur = tot_cur;
        cur = sum(J,2);
        tot_cur = sum(cur);
        plot([bn-0.01,bn], [old_cur,tot_cur])
        xlabel('Bottleneck Size')
        ylabel('Current Density')
    end
    title('Effect of Bottleneck size on Current Density')
end