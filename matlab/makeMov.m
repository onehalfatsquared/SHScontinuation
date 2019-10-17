clear;


load('mov1');load('mov2');load('mov3');

figure(1);
set(0,'DefaultAxesFontSize',20);
v = VideoWriter('diff.mp4','MPEG-4');
v.FrameRate = 30; v.Quality = 95;
open(v);
rho = 10; rho_end = 2; step = 0.1;

rho = rho:-step:rho_end;
%plot the data in subplot
for count=1:length(q2)
    clf
    subplot(1,2,1);
    title('Cluster 1');
    cPlot(q1(count,:)',rho(count));
    subplot(1,2,2);
    title('Cluster 2')
    cPlot(q2(count,:)',rho(count));
    str = strcat("\rho = ",num2str(rho(count)));
    text(-4, 5,str, 'fontsize', 18)
    M = getframe(gcf);
    writeVideo(v,M);
    count = count+1;
end

close(v);