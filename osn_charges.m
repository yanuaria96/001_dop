function osn_charges()
global E
E = 1.0; % electric field

%начальные у�?лови�?
r=1; N=100;%sphere radus
R = linspace(r,r,N)';

%�?оздание �?феры
fig = figure(1);
MySphere(r,0,0,0); hold on;

%генериование �?лучайных точек на �?фере (�?оздание ма�?�?ива)
[theta phi thetad phid] = GeneratePoints(N); %teta and phi array

[x y z] = Sph2Cart(R,theta,phi);
[xd yd zd] = Sph2Cart(R,thetad,phid);

%по�?троить точки на графике 
s(1) = scatter3(x,y,z,30,[0 0 1],'filled'); %точки
% s(2) = scatter3(xd,yd,zd,30,[1 0 0],'filled'); %двойники

dtheta=0.00001*pi/180;
dphi=0.00001*pi/180;

dt = 0.0001*N^(-2);
etta = 1;

%figure(f);
for n = 1:400
%     pause(0.02)
    for k = 1:N
        for r = 1:N
            if r==k
                continue;
            end
            %ча�?тица-ча�?тица
            [theta(k) step] = Moving_theta(dt,etta,...
                theta(k),dtheta,phi(k),R(k),theta(r),phi(r),R(r));
                 
            [phi(k) step] = Moving_phi(dt,etta,...
                phi(k),dphi,theta(k),R(k),phi(r),theta(r),R(r));
            
%             %ча�?тица-двойник
            [theta(k) step] = Moving_theta(dt,etta,...
                theta(k),dtheta,phi(k),R(k),thetad(r),phid(r),R(r));
            
            [phi(k) step] = Moving_phi(dt,etta,...
                phi(k),dphi,theta(k),R(k),phid(r),thetad(r),R(r));
            
%             %двойник-ча�?тица
            [thetad(k) step] = Moving_theta(dt,etta,...
                thetad(k),dtheta,phid(k),R(k),theta(r),phi(r),R(r)); 
            
            [phid(k) step] = Moving_phi(dt,etta,...
                phid(k),dphi,thetad(k),R(k),phi(r),theta(r),R(r));
            
%             %двойник-двойник
            [thetad(k) step] = Moving_theta(dt,etta,...
                thetad(k),dtheta,phid(k),R(k),thetad(r),phid(r),R(r));
            
            [phid(k) step] = Moving_phi(dt,etta,...
                phid(k),dphi,thetad(k),R(k),phid(r),thetad(r),R(r));
        end
    end

    if(rem(n,10)==0)
        delete(s);
        [x y z] = Sph2Cart(R,theta,phi);
        [xd yd zd] = Sph2Cart(R,thetad,phid);   
        s(1) = scatter3(x,y,z,30,[0 0 1],'filled');
        s(2) = scatter3(xd,yd,zd,30,[1 0 0],'filled');
        
        pause(0.01);
        disp('iter');
     end
end

% close(v);
% imwrite(im,map,'test.gif','DelayTime',0,'LoopCount',inf)
% save test

% x_pdf = [-pi:.1:pi];
% y = pdf(pd,x_pdf);
 
% figure
% histogram(phi,'Normalization','pdf')
% line(x_pdf,y)
end


function [a step] = Moving_theta(dt,etta,a,da,b,c,aj,bj,cj)
step = (1/da).*(PSI_theta(a+da,b,c,aj,bj,cj)-PSI_theta(a,b,c,aj,bj,cj));
a = a - (dt/etta).*step;
end

function [a step] = Moving_phi(dt,etta,a,da,b,c,aj,bj,cj)
step = (1/da).*(PSI_phi(a+da,b,c,aj,bj,cj)-PSI_phi(a,b,c,aj,bj,cj));
a = a - (dt/etta).*step;
end

function psi = PSI_theta(thetai,phii,Ri,thetaj,phij,Rj)
%input: E - stifness, const
%       thetai - theta for particle i
%       thetaj - theta for particle j
%       phii - phi for particle i
%       phij - phi for particle j
%output: psi - energy for couple i, j
%global E
myu=[0;1;0];
k=5;
ri=[Ri.*sin(thetai).*cos(phii);Ri.*sin(thetai).*sin(phii);Ri.*cos(thetai)];
rj=[Rj.*sin(thetaj).*cos(phij);Rj.*sin(thetaj).*sin(phij);Rj.*cos(thetaj)];
rij=norm(ri-rj);
psi=1/(rij^2)/Mises_Fisher(k,myu,ri)/Mises_Fisher(k,myu,rj);


% psi_field = ...
%     q1*Ri*sin(thetai)*sin(phii) +...
%     q2*Rj*sin(thetaj)*sin(phij);
% psi = psi + E * psi_field;
end

function psi = PSI_phi(phii,thetai,Ri,phij,thetaj,Rj)
myu=[0;1;0];
k=5;
ri=[Ri.*sin(thetai).*cos(phii);Ri.*sin(thetai).*sin(phii);Ri.*cos(thetai)];
rj=[Rj.*sin(thetaj).*cos(phij);Rj.*sin(thetaj).*sin(phij);Rj.*cos(thetaj)];
rij=norm(ri-rj);
psi=1/(rij^2)/Mises_Fisher(k,myu,ri)/Mises_Fisher(k,myu,rj);
end


function [theta, phi, thetad, phid] = GeneratePoints(N)
theta = pi*rand(N,1);
phi = 2*pi*rand(N,1);

for k = 1:N
    if(theta(k)>=0 && theta(k)<pi/2)
        thetad(k) = -theta(k)+pi;
    elseif(theta(k)>=pi/2 && theta(k)<pi)
        d = pi-theta(k);
        thetad(k) = theta(k)-pi+d*2;
    end
    if(phi(k)<=pi && phi(k)>0)
        phid(k) = phi(k)+pi;
    elseif(phi(k)>=pi && phi(k)<2*pi)
        phid(k) = phi(k)-pi;
    end
end
thetad = thetad';phid = phid';
end


function MySphere(R,x,y,z)
[x,y,z]=sphere;
surf(x,y,z);
colormap([0.9  0.9 0.9]);
end

function [R theta phi] = Cart2Sph(x,y,z)
R = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z/R);
phi = atan(y./x);
end

function [x y z] = Sph2Cart(R,theta,phi)
x=R.*sin(theta).*cos(phi);
y=R.*sin(theta).*sin(phi);
z=R.*cos(theta);
end

function D = GenerateDistribution(N)
sigma = 0.05;
myu = 0;
x = linspace(0, 1, N );
D = NormalDistribution(x, myu, sigma);
D = D / max( D );
end

function D = NormalDistribution(x, myu, sigma)
D = 1./sigma./sqrt(2*pi)*exp(-(x-myu).^2/2./sigma.^2);
end
