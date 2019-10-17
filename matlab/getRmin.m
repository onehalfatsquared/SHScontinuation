%get the minimum distance of a cluster as fn of rho
clear

n = 9; num = 14;
q = zeros(1,50);

% for rho = 1:49
%     c = getCluster(n,rho,num,"MED");
%     p = pDists(c2p(c));
%     r = min(p(p>0));
%     q(rho) = r;
% end

rho = 49;
count = 0;
for cn = 1:263
    c = getCluster(n,rho,cn,"MED");
    p = pDists(c2p(c));
    r = min(p(p>1.05));
    if (abs(r-sqrt(2))<1e-4)
        count = count + 1
    end
end

q(50) = 1;

plot(1:50,q)
xlabel('\rho')
ylabel('Minimum Distance')



