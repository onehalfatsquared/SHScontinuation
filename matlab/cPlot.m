function [] = cPlot(x,rho)
    %Plot cluster of particles
    
    opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',0,'tfs',36, ...
    'iftight',1,'az',20,'el',18,'zth',-0.5,'pos',[7,5,15,15],'rho',rho);

    plotcluster2b(x,opts);
end

