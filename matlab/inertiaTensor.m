function [M] = inertiaTensor(clust)
%Compute the moment of inertia tesnor for a cluster of particles
particles=c2p(clust); N=length(particles); %seperate into particles
xx=0; yy=0; zz=0; %initialize diaganol elements
xy=0; xz=0; yz=0; %initialize off diagonal 
%Compute the matrix elements
for i=1:N
    p=particles(i,:);
    xx=xx+p(2)^2+p(3)^2;
    yy=yy+p(1)^2+p(3)^2;
    zz=zz+p(2)^2+p(3)^2;
    xy=xy-p(1)*p(2);
    xz=xz-p(1)*p(3);
    yz=yz-p(2)*p(3);
end

%form the matrix
M=[xx,xy,xz;xy,yy,yz;xz,yz,zz];
end