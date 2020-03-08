
% Commands
clear;
close all;
clc;
format long;
graphs = false;
nPast = 10;
semi = 40;
deltaT = 86400/2; % s
miniDeltaT = deltaT/semi;

% Constants
AU = 149.5978707E6; % km
ra_M = 1.66601*AU; % km
rp_M = 1.3814*AU; % km
r_E = 1*AU; % km
r_M = (ra_M + rp_M)/2; % km
e_M = (ra_M/rp_M-1)/(ra_M/rp_M+1); % -
mu_S = 1.327124421E11; % km^3/s^2
v_E = sqrt(mu_S/r_E); % km/s
v_M = sqrt(mu_S/r_M); % km/s

% e_min/max & e_vec
e_Int = 0.005;
e_Hoh = (2*r_M)/(r_M+r_E)-1; % -
e_min = 0.2; % -
e_max = 0.9; % -
e = (e_min:e_Int:e_max)'; % -

% Semi-Major Axis
a = r_E./(1-e); % km

% V
V = sqrt(mu_S*((2/r_E)-(1./a))); % km/s
 
% DeltaV 
deltaV = V - v_E; % km/s

% Data Struct
Data = struct( 'State', cell(semi,size(deltaV,1)), 'n', cell(semi,size(deltaV,1)), ...
               'e', cell(semi,size(deltaV,1)), 'a', cell(semi,size(deltaV,1)), ...
               'M', cell(semi,size(deltaV,1)), 'OutputCart', cell(semi,size(deltaV,1)), ...
               'KeplerState', cell(semi,size(deltaV,1)), 'OutputKepl', cell(semi,size(deltaV,1)) );
      
% Actual Fun
for i = 1:size(Data,2)
    for k = 0:semi-1
        Data(k+1,i).n = 1/sqrt(a(i).^3/mu_S); % rad/s        
        Data(k+1,i).M = 0+k*miniDeltaT*rad2deg(Data(k+1,i).n); % rad   
        [x,y,~,vx,vy,~,th] = Kepl2Cart( a(i), e(i), 0, 0, 0, ...
                                          wrapTo360( Data(k+1,i).M(1) ) );
        Data(k+1,i).State(1,:) = [0+k/semi, x, y, vx, vy]; % km - km/s
        Data(k+1,i).KeplerState(1,:) = [0+k/semi, a(i), e(i), th];
        Data(k+1,i).e = e(i); % -
        Data(k+1,i).a = a(i); % km
        t = 0+miniDeltaT*k;
        j = 1;
        fprintf('Run %i\n\n',i);
        while( true )
%             fprintf( 'Run %i,%i\n\n',i,j) ;   
            t = t + deltaT;
            j = j + 1;
            Data(k+1,i).M(j) = Data(k+1,i).M(j-1) + rad2deg(Data(k+1,i).n)*deltaT ;
            [x,y,~,vx,vy,~,th] = Kepl2Cart( Data(k+1,i).a, Data(k+1,i).e, 0, 0, 0, ...
                                          wrapTo360( Data(k+1,i).M(j) ) );
            Data(k+1,i).State(j,:) = [t/deltaT, x, y, vx, vy];
            Data(k+1,i).KeplerState(j,:) = [t/deltaT, a(i), e(i), th];
            if( e(i) < e_Hoh )
                if( Data(k+1,i).M(j) > 180 )
                    Data(k+1,i).KeplerState(j,:) = [];
                    Data(k+1,i).State(j,:) = [];
                    break;
                end
            else
                if( vecnorm(Data(k+1,i).State(j,2:3),2,2) > r_M )||(Data(k+1,i).M(j) > 180)
                    break;
                end
            end
        end
    end
end
clear x y z vx vy vz j i t;

% Plots
if graphs
    for i = 1:size(Data,2)
        for k = 1:semi
            figure( 'Name', sprintf( 'Orbit e=%.2f - %i', Data(k,i).e, k ) );
            grid minor;
            hold on;
            xlabel( 'x [km]' );
            ylabel( 'y [km]' );
            % Sun
            plot( 0, 0, 'o', 'Color', '#EDB120', 'MarkerSize', 20 );
            % Orbits
            ang=0:0.01:2*pi; 
            x=r_E*cos(ang);
            y=r_E*sin(ang);
            plot( x, y, 'LineWidth', 1.5, 'Color', '#0072BD' );
            ang=0:0.01:2*pi; 
            x=r_M*cos(ang);
            y=r_M*sin(ang);
            plot( x, y, 'LineWidth', 1.5, 'Color', '#A2142F' );
            % Trajectory
            plot( Data(k,i).State(:,2), Data(k,i).State(:,3), 'LineWidth', 1.5 );
            set( gca, 'FontSize', 15 );
            axis equal;
        end
    end
end
    
% Data Treatment
for i = 1:size(Data,2)
    for k = 1:semi
        if size(Data(k,i).State,1) > nPast+1
            for j = size(Data(k,i).State,1):-1:nPast+1
                Data(k,i).OutputCart(j,:) = reshape(Data(k,i).State(j-nPast:j,2:5)',[1,(size(Data(k,i).State,2)-1)*(nPast+1)]);
                Data(k,i).OutputKepl(j,:) = reshape(Data(k,i).KeplerState(j-nPast:j,4)',[1,(nPast+1)])./180;
            end
            Data(k,i).OutputCart(1:nPast,:) = [];
            Data(k,i).OutputKepl(1:nPast,:) = [];
            Data(k,i).OutputCart(:,4*(nPast+1)+1) = e(i);
            Data(k,i).OutputKepl(:,nPast+2) = e(i);    
            Data(k,i).OutputKepl(1,nPast+3) = k;
            Data(k,i).OutputKepl(2:size(Data(k,i).OutputKepl,1),nPast+3) = 0;
        end
    end
end
clear i j k;

OvelhaNegra = 40;
DataOvelhaNegra = Data(5,OvelhaNegra);
Data(:,OvelhaNegra) = [];

% Output Data to File

writematrix(vertcat(Data(:).OutputCart),'orbitCart.dat','Delimiter','tab')
writematrix(vertcat(Data(:,:).OutputKepl),'orbitKepl.dat','Delimiter','tab')
writematrix(vertcat(DataOvelhaNegra(:).OutputCart),'orbitCartOvelha.dat','Delimiter','tab')
writematrix(vertcat(DataOvelhaNegra(:).OutputKepl),'orbitKeplOvelha.dat','Delimiter','tab')


 