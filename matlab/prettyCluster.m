function [] = prettyCluster(cluster,n,c,p)
%Plots a given cluster and saves the image how it will be loaded into the
%graphviz script for output. Outputs filename as n(dim)c(cluster num)p(rho).png

     opts = struct('srad',0.1,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
    'scolr',1,'salph',0.7,'ifgrid',0,'iftext',0,...
     'iftight',1,'az',-65,'el',12,'zth',-2.5,'pos',[7,5,15,15]);
    plotcluster2b(cluster',opts);
    filename=strcat('n',num2str(n),'c',num2str(c),'p',num2str(p));
    print(filename,'-deps');
    close;
end

