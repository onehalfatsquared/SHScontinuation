n = 11;
potential = "m";
k = ["LOW","MED","HIGH"];
cluster = 139;
rho = 4;

c1 = getCluster2(n,rho,cluster,k(1),potential);
c2 = getCluster2(n,rho,cluster,k(2),potential);
c3 = getCluster2(n,rho,cluster,k(3),potential);


e1 = MP(c2p(c1),rho,1)
e2 = MP(c2p(c2),rho,1)
e2 = MP(c2p(c3),rho,1)