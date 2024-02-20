clear all
clc
tic
%Inputs
g = 9.806;
Hmax = ncread("H.nc","hmax"); %Maximum wave height
H = Hmax(:)';
H(isnan(H)) = 0;
T = ncread("T.nc","tmax"); %Maximum period
T(isnan(T)) = 0;
T(isinf(T)) = 0;
Tp = T(:)';
d = ncread('d.nc','wmb'); %Depth
d = d(:)';
d(isnan(d)) = 0;

%Wavelength
x = ones(1073,50)-0.9;
lambda = zeros(1,1073);
f = linspace(0.05,0.5,50);
w = 2*pi.*f;

for a = 1:1073
    for b = 1:50
        y(a,b) = (w(b).^2.*d(a))/9.81;
        for i = 1:50
            x(a,b) = y(a,b)/tanh(x(a,b));
        end
    end
end

for a = 1:1073
        d(1,a);
        k(a,:) = x(a,:)/d(1,a);
        k(isinf(k)) = 0;
        k(isnan(k)) = 0;
        lambda = 2*pi./k;
end

% plot(k,w)
% title('Wavenumber variation', 'FontSize',15)
% ylabel('Angular frequency (rad/s)','FontSize',12)
% xlabel('Wavenumber (1/m)','FontSize',12)

k = k';

%JONSWAP spectrum
Hs = ncread("Hs.nc","swh");
Hs = Hs(:)';
Gamma = 3.3;

for n = 1:1073
    fp(n)=Tp(n)^-1;
    fp(isinf(fp)) = 0;
end

for n = 1:1073
    Spectrum = jonswap10b(Hs(n),1/fp(n),Gamma,f);
    S_JS(:,n) = Spectrum(:,2);
    S_JS(isnan(S_JS)) = 0;
end
 
% plot(f,S_JS(:,1000));
% title('JONSWAP spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Wave Spectral Density (m^2/Hz)')

%elevation
X = 0;
t = linspace(-25,25,50); %-25 25
w = 2*pi.*f;
w(isinf(w)) = 0;  
dn = S_JS*(w(2)-w(1));

Eta = zeros(1,1073);
alpha = zeros(1,1073);
Z2 = zeros(1,1073);

for i = 1:1073
    while H(i) > (max(Eta(:,i)) - min(Eta(:,i)))
        for a = 1:50
            A(a,:) = dn(a,i).*cos(-w(a).*t); %Distribution of A with frequency
        end
        Z = sum(dn,1);
        Z1 = sum(A,1);
        for b = 1:50
            Eta(b,i) = (alpha(i)/Z(i))*Z1(b); %Distribution of Etas with time
        end
        alpha(i) = alpha(i) + 0.1;
        Z2(i) = min(Eta(:,i));
    end
end

% plot(t,Eta(:,500))

% Force = reshape(alpha, [37 29]);
% y = surf(double(Force));
% colorbar
% title('eta')

%Velocity

for a = 1:1073
    z(a,:) = linspace(d(a),0,10);
    for b = 1:10
        for c = 1:50
            for i = 1:50
                F(a,b,i,c) = (cosh(k(c,a).*(-z(a,b)))/sinh(k(c,a)*d(a)))*(w(c).*dn(c,a).*cos(-w(c).*t(i)));
                F1(a,b,i,c) = (cosh(k(c,a).*(-z(a,b))))/sinh(k(c,a)*d(a))*(w(c).*dn(c,a).*sin(-w(c).*t(i)));
                F(isinf(F)) = 0;
                F(isnan(F)) = 0;
                F1(isinf(F1)) = 0;
                F1(isnan(F1)) = 0;
            end
        end
    end
end

A1 = sum(F,4);
A2 = sum(F1,4);

for a = 1:1073
    for f = 1:10
        for e = 1:50
            v(e,a,f) = (alpha(a)/Z(a))*A1(a,f,e);
            v_dot(e,a,f) = (alpha(a)/Z(a)).*A2(a,f,e);
        end
    end
end

% Force = reshape(d, [37 29]);
% y = surf(double(Force));
% title('Force spatial variability');

% plot(t,squeeze(v),'-b')

% pcolor(t,v(1),z)

%Morison Equation

Cd = 1.05;
Cm = 1.2;
rho = 1025;
D = 5;

for i = 1:1073
    for a = 1:50
        for b = 1:10
            L(i) = z(i,1)-z(i,2);
            Force(a,i,b) = 0.5*Cd*rho*D*L(i)*v(a,i,b).*abs(v(a,i,b))+Cm*rho*0.25*pi*D^2*L(i)*v_dot(a,i,b);
        end
    end
end

% plot(t,squeeze(Force),'-b')

Force = @(x,t)0.5*Cd*rho*D*L.*v.*abs(v)+Cm*rho*0.25*pi*D^2*(L).*v_dot;
Ftotal = integral(Force,-10,0,'ArrayValued',true,'AbsTol',1e-10,'RelTol',1e-6);

for a = 1:1073
Fmax(a) = max(Ftotal(:,a));
end

toc

% F = 0.5*Cd*rho*D*d.*v.*abs(v) + Cm*rho*0.25*pi*D^2*d.*v_dot;

Force = reshape(Fmax, [37 29]);
y = surf(double(Force));
title('Force spatial variability');
colorbar
xlabel('Distance (km)')
ylabel('Distance (km)')
colormap("cool");
shading FLAT;
% 
% % hold on
% % yyaxis left
% % ylabel('v')
% % plot(t,v)
% % yyaxis right
% % plot(t,F)
% % ylabel('Force (N)')
% % xlabel('Time (s)')
% 
% % plot(F,z)
% % xlabel('Force (N)')
% % ylabel('z (m)')
% % title('Peak force')
 
% %delta - standard deviation of crest elevation 
% %w0 - mean wave frequency 
% %Theta - direction relative to mean wave direction

% for a = 1
%     for b = 1:50
%         for e = 1:200
%             Z3 = sum(B,3);
%             Z4 = sum(B1,3);
%             Z5 = sum(F,3);
%             Z = sum(dn,1);
%             v(e,a,b) = (alpha(a)/Z(a)).*Z3(a,e).*Z5(a,b);
% %             v_dot(e,d,f) = (alpha(d)/Z(d)).*Z4(e,d)*Z5(f);
%         end
%     end
% %     Velocity(d) = max(v(:,d,50));
% end