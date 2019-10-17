function pgo = point_group_order(xyz)

global tol
tol = 0.01;

% Developed by Yakir Forman, July 2016


%%
[m Na] = size(xyz);
if m == 1
    xyz = reshape(xyz,3,Na/3);
    Na = Na/3;
end
xyz0 = xyz; % save the original configuration

% place the center of mass into the origin
cmass = mean(xyz,2);
xyz = xyz - cmass*ones(1,Na);

% find distances from the origin and sort them
d = sqrt(sum(xyz.^2,1));
[dsort,isort] = sort(d,'descend');
xyz = xyz(:,isort);
d = dsort;


% if there are atoms at the center, remove them
while d(Na) < tol
    d(Na) = [];
    xyz(:,Na) = [];
    Na = Na - 1;
end

% group the distances
% Ngroup = the number of groups
% ngroup(k) = the number of atoms in group k
% ifirst(k) = the first atom in the group k
% ilast(k) = the last atom in the group k
Ngroup = 0; 
k = 0;
while k < Na
    Ngroup = Ngroup + 1;
    k = k + 1;
    ifirst(Ngroup) = k;
    while k < Na & abs(d(ifirst(Ngroup)) - d(k + 1)) < tol
        k = k + 1;
    end
    ilast(Ngroup) = k;
    ngroup(Ngroup) = k - ifirst(Ngroup) + 1; 
end
% fprintf('Ngroup = %d\n',Ngroup);
% for j = 1 : Ngroup
%     fprintf('ifirst = %d, ilast = %d, ngroup = %d\n',ifirst(j),ilast(j),ngroup(j));
% end

%Check if any two atoms overlap
overlap = 0;
kk = 0; %the atom before the beginning of the group
for j = 1 : Ngroup
    Ningroup = ngroup(j);
    for jk = 1 : Ningroup-1
        atomA = xyz(:,kk+jk);
        for jj = jk+1 : Ningroup
            atomB = xyz(:,kk+jj);
            if norm(atomA-atomB) < tol
                overlap = 1;
            end
        end
    end
    kk = kk + Ningroup;
end
if overlap == 1
    fprintf('Warning: two atoms were found to coincide.\nProgram cannot compute point group order.\n')
    pgo = -1;
else

% find the smallest group
% sort groups according to the numbers of atoms in them in the ascending order

[gmin, jmin] = min(ngroup);
[gsort, jsort] = sort(ngroup,'ascend');
ngroup = gsort;
ifirst = ifirst(jsort);
ilast = ilast(jsort);

i1 = ifirst(1); %An atom from the minimum size group
pos1 = xyz(:,i1); %Its position vector
%Now test for an opposite atom
axi = 1; %No opposite atom, so there is 1 axial atom
for i = i1+1 : ilast(1)
    posi = xyz(:,i);
    if norm(posi+pos1) < tol %There is an opposite atom
        %Swap opposite atom to be right next to initial
        xyz(:,i) = xyz(:,i1+1);
        xyz(:,i1+1) = posi;
        axi = 2; %There are now 2 axial atoms
        break
    end
end

%Stabilizer computation:
if gmin-axi > 0 %There are non-axial atoms in the minimum size group
    nonax = i1+axi : ilast(1); %Gives indices of nonaxial atoms
    stab = check_sym(xyz,[0,0,0]',pos1,xyz(:,nonax));
else %We need to take a different group as a non-axial block for stab
    nonax = [];
    ig = 1;
    while ig < Ngroup & isempty(nonax)
        ig = ig+1; %Look at the next group
        dist = point_line_distance(xyz(:,ifirst(ig) : ilast(ig)),[0,0,0]',pos1);
                %Find distances from group to axis
        nonax = find(dist > tol); %Gives indices of nonaxial atoms within group
    end
    if isempty(nonax)
        fprintf('Warning: all atoms are along one axis. Point group is not discrete.\n');
        stab = -1;
    else
        nonax = nonax + ones(size(nonax))*(ifirst(ig)-1); %Gives indices of nonaxial atoms
        stab = check_sym(xyz,[0,0,0]',pos1,xyz(:,nonax));
    end
end

%Orbit computation
orb = 1; %i1 itself is in the orbit
dots = pos1'*xyz(:,nonax); %gives dot products of i1 
                %position vector with position vectors of atoms in 
                %non-axial block; will be used in orbit calculations
n = size(dots,2); %Number of atoms in non-axial block
if axi == 2 %check if opposite atom is in the orbit, if it exists
    %Check inversion
    if dist_conf(xyz, -1*xyz) < tol
        orb = orb + 1;
    else
        %Check reflection
        p = xyz(:,i1)/norm(xyz(:,i1)); %Unit vector along axis
        H = eye(3) - 2*p*p';
        if dist_conf(xyz, H*xyz) < tol
            orb = orb + 1;
        else
            %Check 180-degree rotations
            flag180 = 0;
            for i3 = 1:n
                if abs(dots(i3)) < tol %i3 is perpendicular to i1-axis: Check a rotation that stabilizes i3
                    u = xyz(:,nonax(i3));
                    u = u/norm(u);
                    uu = u*u';
                    R = -1*eye(3) + 2*uu; %The rotation matrix around u by angle 180
                    if dist_conf(xyz, R*xyz) < tol
                        flag180 = 1;
                        break
                    end
                end
                for i4 = i3+1 : n
                    if abs(dots(i3) + dots(i4)) < tol %180-degree rotation should map two other atoms onto each other
                        u = xyz(:,nonax(i3))+xyz(:,nonax(i4));
                        normu = norm(u);
                        if normu > tol %i3 and i4 are not on the same axis
                            u = u/normu;
                            uu = u*u';
                            R = -1*eye(3) + 2*uu; % the rotation matrix around axis u by angle 180
                            if dist_conf(xyz, R*xyz) < tol
                                flag180 = 1;
                                break
                            end
                        else % i3 and i4 are on the same axis
                            u = cross(pos1, xyz(:,nonax(i3))); %Try rotating i1 onto opposite, i3 onto i4
                            u = u/norm(u);
                            uu = u*u';
                            R = -1*eye(3) + 2*uu; %The rotation matrix around u by angle 180
                            if dist_conf(xyz, R*xyz) < tol
                                flag180 = 1;
                                break
                            end
                        end
                    end
                end
                if flag180 == 1
                    break
                end
            end
            if flag180 == 1
                orb = orb+1;
            end
        end
    end
end
for ii = 1:gmin-axi
    posii = xyz(:,i1+axi+ii-1);
    %Check reflection
    u = (pos1-posii)/norm(pos1-posii);
    H = eye(3) - 2*u*u';
    if dist_conf(xyz, H*xyz) < tol
        orb = orb + 1;
    else
        %Check 180 degree rotation
        u = pos1+posii;
        u = u/norm(u);
        uu = u*u';
        R = -1*eye(3) + 2*uu; % the rotation matrix around axis u by angle 180
        if dist_conf(xyz, R*xyz) < tol
            orb = orb + 1;
        else
             %Check rotation -- 
             %     find atom i3 such that if R brings i1 to ii,
             %     then -R brings i1 to im
             for im = 1:gmin-axi
                 if ii ~= im & abs(dots(ii)-dots(im)) < tol
                     posim = xyz(:,i1+axi+im-1);
                     u = cross(pos1 - posii, pos1 - posim);
                     u = u/norm(u); %Axis of rotation -- perpendicular to plane with the three atoms
                     axto1 = pos1 - pos1'*u*u; %Vector from axis to 1
                     axtoii = posii - posii'*u*u; %Vector from axis to ii
                     co = (axto1'*axtoii)/(norm(axto1)*norm(axtoii)); %Cosine of angle of rotation
                     si = sqrt(1 - co*co); %Sine of angle of rotation -- doesn't matter if positive or negative
                     uu = u*u';
                     ux = [0,-u(3),u(2); u(3),0,-u(1); -u(2),u(1),0];
                     R = eye(3)*co + ux*si + (1 - co)*uu; % the rotation matrix around axis p by angle theta
                     if dist_conf(xyz, R*xyz) < tol
                         orb = orb + 1;
                         break
                     end
                 end
             end
        end
    end
end
pgo = stab*orb;

 
 %% Conclusion
% fprintf('pgo = %d\n',pgo);

end %ENDS THE ELSE CONDITION ON OVERLAP

end    
  
%%
function nsym = check_sym(xyz,v1,v2,nonaxi) %Checks all symmetries stabilizing v1-v2 axis
                                           %given a block of equidistant nonaxial atoms
global tol
nsym = 1;
p = v2 - v1; %Vector along axis
[m n] = size(nonaxi);
if n > 1
    p = p/norm(p); % unit vector in the direction of rotation 
    uu = p*p';
    ux = [0,-p(3),p(2); p(3),0,-p(1); -p(2),p(1),0];
    theta = 2*pi/n;
    co = cos(theta); si = sin(theta);
    R = eye(3)*co + ux*si + (1 - co)*uu; % the rotation matrix around axis p by angle theta
    %check rotations
    xyz1 = xyz;
    for j = 1 : n - 1
        xyz1 = R*xyz1;
        if dist_conf(xyz,xyz1) < tol
            nsym = nsym + 1;
        end
    end
end
% check reflections
%An arbitrary atom is chosen as anchor (last in list)
anchor = nonaxi(:,n);
%First, we check reflections in plane that includes v1,v2,anchor
u = cross(p, anchor - v1);
u = u/norm(u); %now gives unit normal to plane including v1,v2,anchor
H = eye(3) - 2*u*u';
xyz1 = H*xyz;
if dist_conf(xyz,xyz1) < tol
    nsym = nsym + 1;
end

%Next, we check if anchor can be reflected to any other atoms
for j = 1 : n-1
    u = nonaxi(:,j)-anchor;
    if abs(p'*u) < tol %Reflection should stabilize v2-v1 axis, so p should be in plane of reflection (perpendicular to u)
        u = u/norm(u);
        H = eye(3) - 2*u*u';
        if dist_conf(xyz,H*xyz) < tol
            nsym = nsym+1;
        end
    end
end

end

%%
function dist = dist_conf(x,y)
[m Na] = size(x);
r = zeros(Na,1);
for k = 1 : Na
    r(k) = min(sqrt(sum((x - y(:,k)*ones(1,Na)).^2,1)));
end
dist = max(r);
end

%%
function dist = point_line_distance(pp,v1,v2)
[m,n] = size(pp);
dist = zeros(1,n);
u = v2 - v1;
uu = u'*u;
for j = 1 : n
    v = pp(:,j) - v1;
    vv = v'*v;
    uv = u'*v;
    dist(j) = sqrt(vv - uv^2/uu);
end
end                 
