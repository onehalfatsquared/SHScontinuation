%function [r1,r2] = rotateClusters(c1,c2)
    %function to compute optimal rotation of one set of particles onto
    %another
    
    rho = 50
    c1 = getCluster(10,35,6,'LOW')'
    c2 = getCluster(10,17,6,'LOW')'
    
%     [~,b,c,~,~,~] = testSame(c1,c2,1)
%     
%     p1 = b'; c1=p1(:);
%     p2 = c'; c2=p2(:);
    
    figure(1)
    cPlot(c1,40)
    figure(2)
    cPlot(c2,17)
