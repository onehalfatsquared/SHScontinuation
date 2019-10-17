function [c] = getCluster(n,rho, num, E)
    %get cluster from the file its written to at given N and rho
    if (rho==50)
        filename=strcat('n',num2str(n),'adjust.txt'); %File name nxrhoy.txt
    else
         filename=strcat('n',num2str(n),'rho',num2str(rho),'Sticky',E,'.txt'); %File name nxrhoy.txt
    end
    fileID=fopen(filename,'r');             %Open file
    formatSpec='%f';                        %Format as float
    data=fscanf(fileID,formatSpec);         %read in data    
    fclose(fileID);                         %Close file

    %Seperate column vector into an array of clusters, num clusters by 3n.
    num_clusters=length(data)/(3*n);
    clusters=reshape(data,[3*n,num_clusters])';
    
    c = clusters(num,:);
end

