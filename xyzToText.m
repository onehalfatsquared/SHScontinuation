n=12;                                    %Number of particles
filename=strcat('minlist',num2str(n),'.txt'); %files are nx.txt s.t. x= num part.
formatSpec='tab';                        %Format as float
data=tdfread(filename,formatSpec);         %read in data    
usable=data.X;

nums=[];
count=1; errors=0;
while errors<1
    try
        if usable(count,1)=='S'
          d=strsplit(usable(count,:));
         nums=[nums str2num(d{2}) str2num(d{3}) str2num(d{4})];
        end
        count=count+1;
    catch
        errors=errors+1;
    end
end

filename=strcat('n',num2str(n),'Morse30.txt'); %File name nxrhoy.txt
fileID=fopen(filename,'w');                     %Open file
fprintf(fileID,'%.16f \n',nums);      %Write cluster
fclose(fileID);