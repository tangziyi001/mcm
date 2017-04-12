function res = PSDCA(v_max,K,L,N,p,ps,Te,sd,opOut)
% PURPOSE: Traffic simulation using the Partial Self-Driving Cellular
% Automata model
% ------------------------------------------------------------
% SYNTAX: res = PSDCA(v_max,K,L,N,p,Te,sd);
% ------------------------------------------------------------
% INPUT: 
%        v_max: 1x1 -> maximum speed (e.g. legal limit)
%        K: 1x1 -> number of lanes
%        L: 1x1 -> length of the (circular) road
%        N: 1x1 -> number of vehicles (N <= L)
%        p: 1x1 -> probability of sudden stop
%        ps: 1x1 -> probability that a non self-driving car changes lane
%                   when needed
%        Te: 1x1 -> effective number of time steps
%        sd: 1x1 -> percentage of self-driving cars
%        opOut: 1x1 -> 1 produces a summary of the simulation in ASCII file
%        PSDCA.out, 0 does not produce any output. Optional. Default=0.
% ------------------------------------------------------------
% OUTPUT: A structures with some of the inputs (v_max, L, N, p, Te) and
%               rho: 1x1 -> density = N/L
%               v: Tex1 -> velocity across time
%               v_mean: 1x1 -> mean velocity
%               XV: Te x k x L -> position and velocity of vehicles across time
% ------------------------------------------------------------
% NOTE: The programs runs Tr=10*L additional time steps to form a
% relaxation interval before reaching the stationary state. The program
% retains only the last Te time steps.
% ------------------------------------------------------------
% Library: select_uniform

% written by:
%  Enrique M. Quilis
%  <equilis@gmail.com>

% Version 3.2 [December 2015]

% ------------------------------------------------------------
% Housekeeping

% Setting clock
t0 = clock;

% Summary in external text file
if (nargin == 5)
    opOut = 0;
end

% Code for missing
keyNA = -99999;
% Separators for output file
sep ='----------------------------------------------------';
sep1='****************************************************';

% ------------------------------------------------------------
% Checks
if (N > L*K)
    error('*** Number of vehicles must be lower or equal than the size of the road ***');
end

% Density
rho = N/(L);

% Middle point of the road
Lc = fix(L/2);

% Relaxation time (burning sample)
Tr = 0*L;
% Total time
T = Tr + Te;

% Number of self-driving cars
SDC = floor(N*sd);

% Number of self-driving cars
NSDC = N-SDC;
% Allocation and initial velocities and positions
v = zeros(T,1);
pos_sdc = zeros(1,SDC);
pos_nsdc = zeros(1,NSDC);
% Matrix of velocities and type of cars
XV = keyNA*ones(T,K,L); % State of the road
XT = zeros(T,K,L); % Type of car

%% ------------------------------------------------------------
% Initialization

% Initial positions for self-driving cars
aux = 1:L*K;
for j=1:SDC    
    u = select_uniform(1,length(aux));
    pos_sdc(j) = aux(u);
    aux(u) = [];
    % Start velocity is 0
    [x,y] = idxMod2D(pos_sdc(j),L);
    % fprintf('Self Val Init %d %d\n', x,y);
    XV(1,x,y) = 0;
    XT(1,x,y) = 1;
end

% Initial positions for non self-driving cars
for j=1:NSDC        
    u = select_uniform(1,length(aux));
    pos_nsdc(j) = aux(u);
    aux(u) = [];
    % Start velocity is 0
    [x,y] = idxMod2D(pos_nsdc(j),L);
    % fprintf('Non Self Init %d %d\n', x,y);
    XV(1,x,y) = 0;
    XT(1,x,y) = 0;
end

% Positions of all cars
pos = [pos_sdc pos_nsdc];

[d_f,d_lf,d_lb,v_f,v_lf,v_lb,type_f,type_lf,type_lb] = updateDistance(N,L,K,XV,XT,pos,1,keyNA);
vis = squeeze(XV(1,:,:));
disp(vis);
%% ------------------------------------------------------------
% Simulation

for t=2:T
    fprintf('\nStep %d\n', t);
    % Extract velocities and types of cars observed at t-1
    vel = squeeze(XV(t-1,:,:));
    types = squeeze(XT(t-1,:,:));
    % Check potential slow-down for each self_driving car
    slow = zeros(K,L);
    for j=1:N
        [lane,loc] = idxMod2D(pos(j),L);
        if types(lane,loc) == 1
            if vel(lane,loc) > d_f(j)
                % If front car is NSDC, change lane or avoid crashes
                if type_f(j) == 0 && checkSTCA(v_max,vel(lane,loc),d_f(j),d_lf(j),d_lb(j)) == 0
                    slow(lane,loc) = 1;
                end
            end
        end
    end
    for k=1:N
        for j=1:N
            [lane,loc] = idxMod2D(pos(j),L);
            if types(lane,loc) == 1 && type_f(j) == 1
                if vel(lane,loc) > d_f(j) && slow(lane, idxMod1D(loc+d_f(j)+1,L)) == 1
                    slow(lane,loc) = 1;
                end
            end
        end
    end
    % Loop on cars
    for j=1:N
        [lane,loc] = idxMod2D(pos(j),L);
        % fprintf('Now Exe %d %d %d\n',lane,loc,types(lane,loc));
        % Acceleration
        v1 = min(vel(lane,loc)+1, v_max);
        % If the current car is a self-driving car
        if types(lane,loc) == 1
            % check the distance to the front car
            if v1 > d_f(j)
                % If front car is NSDC, change lane or avoid crashes
                if type_f(j) == 0
                    % Change lane if STCA is satisfied
                    if checkSTCA(v_max,v1,d_f(j),d_lf(j),d_lb(j)) == 1
                        newlane = idxMod1D(lane-1,K);
                        % fprintf('Newlane %d\n',newlane);
                        XV(t,newlane,idxMod1D(loc+v1,L)) = v1;
                        XT(t,newlane,idxMod1D(loc+v1,L)) = 1;
                        pos(j) = idxMul2D(newlane,loc+v1,L);
                        continue;
                    else
                        v1 = min(v1, d_f(j));
                    end
                % If front car is SDC, keep going or join the
                % cluster of the front car
                else
                    if slow(lane,loc) == 1
                        XV(t,lane,idxMod1D(v1+d_f(j),L)) = v1;
                        XT(t,lane,idxMod1D(v1+d_f(j),L)) = 1;
                        pos(j) = idxMul2D(lane,d_f(j)+v1,L);
                        continue;
                    end
                    % If the current car is potentially surpass front car
                    % then current car will immediately follow front car
                    % with the same speed
                    if v1 > v_f(j)+1+d_f(j)
                        v1 = v_f(j)+1;
                        XV(t,lane,idxMod1D(v1+d_f(j),L)) = v1;
                        XT(t,lane,idxMod1D(v1+d_f(j),L)) = 1;
                        pos(j) = idxMul2D(lane,d_f(j)+v1,L);
                        continue;
                    end
                end
            end
            XV(t,lane,idxMod1D(loc+v1,L)) = v1;
            XT(t,lane,idxMod1D(loc+v1,L)) = 1;
            pos(j) = idxMul2D(lane,loc+v1,L);
        else
            if v1 > d_f(j)
                % Random Lane Change
                u = rand;
                if checkSTCA(v_max,v1,d_f(j),d_lf(j),d_lb(j)) == 1 && u < ps
                        newlane = idxMod1D(lane-1,K);
                        % fprintf('Newlane %d\n',newlane);
                        XV(t,newlane,idxMod1D(loc+v1,L)) = v1;
                        XT(t,newlane,idxMod1D(loc+v1,L)) = 0;
                        pos(j) = idxMul2D(newlane,loc+v1,L);
                        continue;
                else
                    v1 = min(v1, d_f(j));
                end
            end
            % Random sudden deceleration
            u = rand;
            % Updating velocity
            if (u < p)
                v1 = max(v1-1,0);
            end
            XV(t,lane,idxMod1D(loc+v1,L)) = v1;
            XT(t,lane,idxMod1D(loc+v1,L)) = 0;
            pos(j) = idxMul2D(lane,loc+v1,L);
        end
    end
    vis = squeeze(XV(t,:,:));
    disp(vis);
    I = find(vis ~= keyNA);
    vis = vis(I);
    v(t) = mean(vis);
    [d_f,d_lf,d_lb,v_f,v_lf,v_lb,type_f,type_lf,type_lb] = updateDistance(N,L,K,XV,XT,pos,t,keyNA);
end %of temporal loop

% Burning initial Tr observations
v(1:Tr) = [];
XV(1:Tr,:) = [];

% Global mean
v_mean = mean(v);
flow_mean = v_mean*rho;

% ------------------------------------------------------------
% Elapsed time
et = etime(clock,t0);

% ------------------------------------------------------------
% Generating output

switch opOut
    case 0
        % No output
    case 1
        % Output in file NaSch.out
        fid = fopen('NaSch.out','w');
        fprintf(fid,'\n ');
        fprintf(fid,'NaSch MODEL FOR TRAFFIC SIMULATION \n');
        fprintf(fid,'% s\n',sep);
        fprintf(fid,' Speed limit: %4d\n ',v_max);
        fprintf(fid,'Length of the (circular) road: %4d\n ',L);
        fprintf(fid,'Number of vehicles: %4d\n ',N);
        fprintf(fid,'Probability of sudden stop: %8.4f\n ',p);
        fprintf(fid,'Effective time: %4d\n ',Te);
        fprintf(fid,'Relaxation time: %4d\n ',Tr);
        fprintf(fid,'Total time: %4d\n ',T);
        fprintf(fid,'%s \n',sep);
        fprintf(fid,' Density: %8.4f \n ',rho);
        fprintf(fid,'Mean velocity: %8.4f\n ',v_mean);
        fprintf(fid,'%s \n',sep);
        fprintf(fid,' Elapsed time: %8.4f\n ',et);
        fprintf(fid,'%s \n',sep);
        fprintf(fid,'\n ');
        fclose(fid);
    otherwise
        error('*** opOut must be 0 or 1 ***');
end

% ------------------------------------------------------------
%% Loading structure
res.L = L;
res.N = N;
res.v_max = v_max;
res.p = p;
res.Te = Te;
res.rho = rho;
res.v = v;
res.v_mean = v_mean;
res.flow_mean = flow_mean;
res.XV = XV;
end

%% ------------------------------------------------------------
% Helper Functions
function flag = checkSTCA(v_max, v1, d_f, d_lf, d_lb)
    if d_lf > d_f && d_lb > 1+v_max-v1
        flag = 1;
    else
        flag = 0;
    end
end

function [d_f,d_lf,d_lb,v_f,v_lf,v_lb,type_f,type_lf,type_lb] = updateDistance(N,L,K,XV,XT,pos,t,keyNA)
    % Initial distances (number of cells between vehicles)
    % Distances to the front car
    d_f = zeros(1,N);
    % Velocity of the front car
    v_f = zeros(1,N);
    % Type of the front car
    type_f = zeros(1,N);
    % Distances to the left front car
    d_lf = zeros(1,N);
    % Velocity of the left front car
    v_lf = zeros(1,N);
    % Type of the left front car
    type_lf = zeros(1,N);
    % Distances to the left back car
    d_lb = zeros(1,N);
    % Velocity of the left back car
    v_lb = zeros(1,N);
    % Type of the left back car
    type_lb = zeros(1,N);

    for i=1:N
        [lane,loc] = idxMod2D(pos(i),L);
        %fprintf('\nDistances for %d: %d\n',i,pos(i));
        % Distance to the closest front Car
        front_min = inf;
        for j=1:L
            if XV(t,lane,j) == keyNA
                continue;
            end
            if mod(j-loc+L, L) < front_min
                % Cannot be self
                if mod(j-loc+L, L) == 0
                    continue;
                end
                front_min = mod(j-loc+L, L);
                v_f(i) = XV(t,lane,j);
                type_f(i) = XT(t,lane,j);
            end
        end
        d_f(i) = front_min;
        %fprintf('Front Min %d\n', d_f(i));
        % Distance to the closest left front Car
        left_front_min = inf;
        for j=1:L
            if XV(t,idxMod1D(lane-1,K),j) == keyNA
                continue;
            end
            if mod(j-loc+L, L) < left_front_min
                left_front_min = mod(j-loc+L, L);
                v_lf(i) = XV(t,idxMod1D(lane-1,K),j);
                type_lf(i) = XT(t,idxMod1D(lane-1,K),j);
            end
        end
        d_lf(i) = left_front_min;
        %fprintf('Left Front Min %d\n', d_lf(i));
        % Distance to the closest left back Car
        left_back_min = inf;
        for j=1:L
            if XV(t,idxMod1D(lane-1,K),j) == keyNA
                continue;
            end
            if mod(loc-j+L, L) < left_back_min
                left_back_min = mod(loc-j+L, L);
                v_lb(i) = XV(t,idxMod1D(lane-1,K),j);
                type_lb(i) = XT(t,idxMod1D(lane-1,K),j);
            end
        end
        d_lb(i) = left_back_min;
        %fprintf('Left Back Min %d\n', d_lb(i));
    end
    % T move from gradient d to distance we have to subtract 1
    d_f = d_f - 1;
    d_lf = d_lf - 1;
    d_lb = d_lb - 1;
end

function [x,y] = idxMod2D(a,l)
    x = floor((a-1)/l)+1;
    y = mod(a-1,l)+1;
end
function re2 = idxMul2D(x,y,l)
    re2 = mod((x-1),l)*l + mod((y-1),l) + 1;
end
function re1 = idxMod1D(a,k)
    re1 = mod(a-1,k)+1;
end
