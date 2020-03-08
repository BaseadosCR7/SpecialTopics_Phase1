function [x, y, z, u, v, w, th] = Kepl2Cart(a, e, i, Om, om, M)

%Converts from Keplerian parameters to Cartesian coordinates (angles in degrees)

%Necessary to present the required number of significant digits
format long;

%Used constants:
miu = 1.327124421E11; %km^3/s^2

%Eccentric anomaly (Equation (11)):
E_aux = (M*pi/180);           %Initial guess
while 1 < 2                   %Iterative process
    E = E_aux + ((M*pi/180) - E_aux + e*sin(E_aux))/(1 - e*cos(E_aux));
    if (abs(E-E_aux)) < 1E-14 %Stop iterating when the absolute diference is  
        break;                %smaller than 1E-15 radian
    end
    E_aux = E;
end

E = E/pi*180; %Converting the result to degrees

%True anomaly (Equation (9)):
th = 2*atand(sqrt((1+e)/(1-e)) * tand(E/2));
    if th < 0  %Guarantees that the result is positive
        th = th + 180;
    end
    if E > 180 %Guarantees that E and theta are both on the 1st/2nd or 3rd/4th quadrants simultaneously 
        th = th + 180;        
    end
    
%Auxiliary constants (Equation (13)):
l_1 = cosd(Om)*cosd(om) - sind(Om)*sind(om)*cosd(i);
l_2 = -cosd(Om)*sind(om) - sind(Om)*cosd(om)*cosd(i);
m_1 = sind(Om)*cosd(om) + cosd(Om)*sind(om)*cosd(i);
m_2 = -sind(Om)*sind(om) + cosd(Om)*cosd(om)*cosd(i);
n_1 = sind(om)*sind(i);
n_2 = cosd(om)*sind(i);

%Auxiliary matrix (Equation (14)):
A = [l_1 l_2;
     m_1 m_2;
     n_1 n_2];
 
%Distance to central body (Equation (12)):
r = (a*(1-e^2))/(1+e*cosd(th));

%Auxiliary coordinates (Equation (14)):
B = [r*cosd(th);
     r*sind(th)];

%Cartesian coordinates (Equation (14)):
C = A*B;
x = C(1);
y = C(2);
z = C(3);

%Specific angular momentum (Equation (15)):
H = sqrt(miu*a*(1-e^2));

%Cartesian velocities (Equation (16)):
u = (miu/H)*(-l_1*sind(th) + l_2*(e + cosd(th)));
v = (miu/H)*(-m_1*sind(th) + m_2*(e + cosd(th)));
w = (miu/H)*(-n_1*sind(th) + n_2*(e + cosd(th)));


% %Displaying results:
% disp('x = ');
% disp(x);
% disp('y = ');
% disp(y);
% disp('z = ');
% disp(z);
% disp('u = ');
% disp(u);
% disp('v = ');
% disp(v);
% disp('w = ');
% disp(w);

end

