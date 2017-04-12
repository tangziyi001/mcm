clear;clc;
v_max = round(60*0.28/3.5);
sd = 0.5;
Te = 20;

filename = 'data.csv';
M = csvread(filename);
% length = round((M(:,3)-M(:,2))*1600/3.5);
lane = M(:,5)+M(:,6)
double_lane = [M(:,5) lane;M(:,6) lane];
double_vol = [M(:,4);M(:,4)];
double_vol = double_vol.*double_lane(:,1)./double_lane(:,2);
% v_t = 60*0.28/3.5; % max velocity for a car in units
% tt = length/(v_t); % time for a max speed car to pass by
vol = double_vol/(24*60);  %per minute
length = (M(:,3)-M(:,2));
length = [length;length];
density = vol./length;
density = log(density);

% figure
% histogram(density,200);
% title('Distribution of Density');
% xlabel('log(Density)'); % x-axis label
% ylabel('Frequecy');

ll = 100;
vol = round(density.*ll/10);

% Different Lane
lane = double_lane(:,1);
I = find(lane == 3);
vol_tmp = vol(I);
lane_tmp = lane(I);

for j=1:size(I,1)
    eff = zeros(size(I,1),1);
    eff2 = zeros(size(I,1),1);
    fprintf('Row %d Vol %d\n', j, vol(j));
    res = PSDCA(v_max,lane_tmp(j),100,vol_tmp(j),0.5,0.5,Te,sd,0);
	eff41(j) = res.flow_mean;
	eff42(j) = res.v_mean;
end
figure
scatter(vol_tmp, eff41,'g','filled');
hold on
scatter(vol_tmp, eff42,'b','filled');
title('Mean Flow and Velocity For lane = 2');
xlabel('log(rho)'); % x-axis label
ylabel('mean flow and velocity');
legend('Mean Flow','Mean Velocity');
fprintf('Mean Change of Flow %d\n', mean(eff41));
