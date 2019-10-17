clear; clf; 
n=9
num = 1

%get cluster from the file its written to at given N and rho
filename=strcat('n',num2str(n),'adjust.txt'); %File name nxrhoy.txt
fileID=fopen(filename,'r');             %Open file
formatSpec='%f';                        %Format as float
data=fscanf(fileID,formatSpec);         %read in data
fclose(fileID);                         %Close file

%Seperate column vector into an array of clusters, num clusters by 3n.
num_clusters=length(data)/(3*n);
clusters=reshape(data,[3*n,num_clusters])';

c = clusters(num,:);

figure(1)
cPlot(c')
rmin(c,49)



function rc = rmin(cluster, rho)
    %get rmin value
    comp_val = [1.01, 1.06];
    if rho>30
        comp = 1;
    else
        comp = 2;
    end
    
    di=unique(pDists(c2p(cluster))); di=di(di>comp_val(comp)); %use 1.01 for rho 38, 1.05 for rho 6
    rc = min(di);
end