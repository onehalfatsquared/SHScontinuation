function b=predictMerge(c1, c2)
%predict if clusters c1 and c2 will merge

p1 = c2p(c1); p2 = c2p(c2);

%compute the rmin values
d1=unique(pDists(p1)); d1=d1(d1>1.01); r1=min(d1);
d2=unique(pDists(p2)); d2=d2(d2>1.01); r2=min(d2);

%compute an optimal rotation and permutation of cluster 2 to cluster 1
[~,p1R,p2R,~,~,~]=testSame(c1,c2,1);

%compute adjacency matrices augmented by rmin
A1 = adjMat(p1R, r1); A2 = adjMat(p2R, r2);

%sum over entries of difference. return 0 if same matrix
b = sum(sum(A1 - A2)); 
end