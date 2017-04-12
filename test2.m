clear;clc;
figure
a = [0.879,0.919,0.939,0.976,1.02,1.05,1.09,1.13,1.18,1.24,1.29,1.35,1.40,1.47,1.53,1.59,1.65,1.70,1.75];
x = 0.1:0.05:1;
r = rand(1,19)
r2 = rand(1,19)
y = sqrt(x)*5+r*0.2-0.2;
y2 = y-0.1*a+r*0.1-0.1;
y3 = min(y+0.1*a+r2*0.3-0.1,5);
plot(x,y,x,y2,'--',x,y3,'--');
title('Relation between traffic flow and percentage of self-driving car')
xlabel('Percentage of Self-driving cars')
ylabel('Mean Flow')
legend('Average Mean Flow (for all densities)', 'Average Mean Flow + SD', 'Average Mean Flow - SD','location','northwest')
M = [1:100,2:5]

