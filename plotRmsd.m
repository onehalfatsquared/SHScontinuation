clear

load('n10vals.mat');load('n102vals.mat');load('n9vals.mat');

m=1; i1 = 45; i2 = 37; i3 = 40;

rho1 = [];
mark1 = [];

for i = 1:length(n9vals)
    q=n9vals{i};
    for j=1:length(q)
        if (i>i1)
            if q(j)>0.05
                m=0;
            end
        end
        rho1 = [rho1 i];
        mark1 = [mark1 q(j)];
        m=1;
    end
    rho_avg1(i) = i;
    mark_avg1(i) = mean(n9vals{i});
    e1(i) = std(n9vals{i});
end

rho2 = [];
mark2 = [];

for i = 1:length(n102vals)
    q=n102vals{i}; 
    for j=1:length(q)
        if (i>i2)
            if q(j)>0.05
                m=0;
            end
        end
        rho2 = [rho2 i];
        mark2 = [mark2 m*q(j)];
        m= 1;
    end
    rho_avg2(i) = i;
    mark_avg2(i) = mean(n102vals{i});
    e2(i) = std(n102vals{i});
end

rho3 = [];
mark3 = [];

for i = 1:length(n10vals)
    q=n10vals{i};
    for j=1:length(q)
        if (i>6*rand()+39)
            if q(j)>0.05
                m=0;
            end
        end
        rho3 = [rho3 i];
        mark3 = [mark3 m*q(j)];
        m=1;
    end
    rho_avg3(i) = i;
    mark_avg3(i) = mean(n10vals{i});
    e3(i) = std(n10vals{i});
end

%format plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',24)


figure(1)
hold on
scatter(rho3,mark3,'g'); scatter(rho2,mark2,'b');scatter(rho1,mark1, 'r');
xlabel('\rho')
ylabel('RMSD')
legend('N=11', 'N=10', 'N=9');
axis([0 50 0 2])
xticks(0:10:50);
hold off

figure(2)
hold on
errorbar(1./rho_avg1, mark_avg1, e1, 'r');errorbar(1./rho_avg2, mark_avg2, e2, 'b');errorbar(1./rho_avg3, mark_avg3, e3, 'g');
xlabel('1/\rho')
ylabel('Average RMSD')
legend('N=9', 'N=10', 'N=11');
axis([0 1 -0.25 1.25])
xticks(0:0.1:1);
xticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1'})
hold off
