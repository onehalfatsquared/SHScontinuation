function d=euDist(particles,i,j)
    %Compute the euclidean distance between particles i and j
    d=particles(j,:)-particles(i,:);
    d=norm(d);
    %d=DNorm2(d,2);    %bad
    %d=sqrt(sum(d .* d, 2));   %bad
    %d=sqrt(mtimesx(d,d'));
end