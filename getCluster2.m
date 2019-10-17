function [c] = getCluster2(n,rho, num, kappa, potential)
    %get cluster from the file its written to at given N and rho - specify
    %kappa as well
    filename=strcat('n',num2str(n),potential,num2str(rho),'Sticky',kappa,'.txt'); %File name nxrhoy.txt
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);         %read in data    
    fclose(fileID);                         %Close file

    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    c = clusters(num,:);
end

