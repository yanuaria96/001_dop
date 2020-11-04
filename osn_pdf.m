function osn_pdf()
% PDF - probability distribudion function
t=tic;
global SPHERE_RADIUS D_THETA D_PHI ETTA D_T
N = 20;
c1_const = 0.0001;
SPHERE_RADIUS = 1.0;
D_THETA = c1_const*pi/180;
D_PHI = c1_const*pi/180;
D_T = .04*N^(-2);
ETTA = 1;

[particles, duplicates] = GeneratePoints(N);
h = MySphere(SPHERE_RADIUS,0,0,0);
h = VisualizePoints([particles duplicates],h);

for iter = 1:10000
    % Find potentials for each particle and duplicate
    maxVelocity = 0;
    sumPsi = 0;
    for k = 1:N
        % Current objects:
        particle = particles(k);
        duplicate = duplicates(k);
        
        % For the "other objects" choose all objects except the particle and
        % the duplicate number "k"
        otherIndices = ones(1,N);
        otherIndices(k) = 0;
        otherIndices = logical(otherIndices);
        otherObjects = [particles(otherIndices) duplicates(otherIndices)];
        
        %
        [particles(k).velocity, sumPsi1] = FindVelocity(particle, otherObjects);
        [duplicates(k).velocity, sumPsi2] = FindVelocity(duplicate, otherObjects);
        
        sumPsi = sumPsi + sumPsi1 + sumPsi2;        
        v2 = particles(k).velocity.theta.^2 + particles(k).velocity.phi.^2;
        maxVelocity = max([v2 maxVelocity]);
    end
   
    if maxVelocity < 3e-4
        D_T = .04*N^(-2)*(maxVelocity/3e-4)^(-1/2);
    end
    
    if maxVelocity < 3e-10
        break;
    end
    
    % Move each object according to its velocity
    for k = 1:N
       particles(k) = MoveObject(particles(k));
       duplicates(k) = MoveObject(duplicates(k));
    end
    
    if rem(iter,1) == 0
        h = VisualizePoints([particles duplicates], h);
        pause(0.01);
    end
    
    if rem(iter,1) == 0
        disp([...
            'iter: ' num2str(iter) ...
            ', maxVelocity: ' num2str(maxVelocity) ...
            ', sumPsi: ' num2str(sumPsi)]);
    end
    
    if rem(iter,1) == 0
        SaveData(particles);
    end
end
tt = toc(t);
disp(tt);
SaveData(particles);

end

function SaveData(points)
aTheta = [];
aPhi = [];
aX = [];
aY = [];
aZ = [];
for k = 1:length(points)
   aTheta(k) = points(k).theta;
   aPhi(k) = points(k).phi;
   [aX(k), aY(k), aZ(k)] = Sph2Cart(1, aTheta(k), aPhi(k));
end
dlmwrite('200_anis_x_y_z.txt', [aX' aY' aZ'],'delimiter','\t');
save 'osn_pdf_output';
end

function obj = MoveObject(obj)
    v = obj.velocity;
    obj.phi = obj.phi + v.phi;
    obj.theta = obj.theta + v.theta;
end

function [velocity, sumPsi] = FindVelocity(object, otherObjects)
global D_THETA D_PHI ETTA D_T 
N = length(otherObjects);
th1 = object.theta;
ph1 = object.phi;
velocity.theta = 0;
velocity.phi = 0;
sumPsi = 0;
for k = 1:N
    th2 = otherObjects(k).theta;
    ph2 = otherObjects(k).phi;
    
    psiUnmoved = PSI(th1,ph1,th2,ph2);
    
    dPsidTheta = - (D_T/ETTA).*(1/D_THETA).*(PSI(th1+D_THETA,ph1,th2,ph2) - psiUnmoved);
    velocity.theta = velocity.theta + dPsidTheta;
    
    dPsidPhi = - (D_T/ETTA).*(1/D_PHI).*(PSI(th1,ph1+D_PHI,th2,ph2)-psiUnmoved);    
    velocity.phi = velocity.phi + dPsidPhi;
    
    sumPsi = sumPsi + psiUnmoved;
end
end

function psi = PSI(thetai,phii,thetaj,phij)
myu=[0;0;1];
k=2;
ri=[sin(thetai).*cos(phii);sin(thetai).*sin(phii);cos(thetai)];
rj=[sin(thetaj).*cos(phij);sin(thetaj).*sin(phij);cos(thetaj)];
rij=norm(ri-rj);
%psi=1/(rij^2);
%psi=1/(rij^2)/sqrt(Mises_Fisher(k,myu,ri))/sqrt(Mises_Fisher(k,myu,rj));
%psi=1/(rij^2)/sqrt(Valee_Poussin(k,myu,ri))/sqrt(Valee_Poussin(k,myu,rj));
psi=1/(rij^2)/sqrt(PDF_Matrix(ri))/sqrt(PDF_Matrix(rj));
end

function [particles, duplicates] = GeneratePoints(N)
for k = 1:N
    theta = pi * rand();
    phi = 2.0 * pi * rand();  
    
    [thetad, phid] = GetMirrored(theta,phi);
    
    particles(k) = struct('theta', theta, 'phi', phi, 'velocity', 0.0, 'isDuplicate', 0);
    duplicates(k) = struct('theta', thetad, 'phi', phid, 'velocity', 0.0, 'isDuplicate', 1);
end  
end

function [theta_mirrored, phi_mirrored] = GetMirrored(theta, phi)
    if(theta >= 0 && theta < pi/2)
        theta_mirrored = -theta + pi;
    elseif(theta >= pi/2 && theta < pi)
        d = pi-theta;
        theta_mirrored = theta - pi + d * 2;
    end
    if(phi <= pi && phi > 0)
        phi_mirrored = phi + pi;
    elseif(phi >= pi && phi < 2 * pi)
        phi_mirrored = phi - pi;
    end
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


function handle = VisualizePoints(points,handle)
delete(handle.scatter);
hold on;
N = length(points);
X = [];
Y = [];
Z = [];
Xd = [];
Yd = [];
Zd = [];
for k = 1:N
    [x,y,z] = Sph2Cart(1,points(k).theta,points(k).phi);
    
    if ~points(k).isDuplicate
        X = [X x];
        Y = [Y y];
        Z = [Z z];
    else
        Xd = [Xd x];
        Yd = [Yd y];
        Zd = [Zd z];
    end
end
subplot(handle.subplot(1));
handle.scatter(1) = scatter3(X,Y,Z,30,[0 0 1],'filled');
handle.scatter(2) = scatter3(Xd,Yd,Zd,30,[1 0 0],'filled');
end
function handles = MySphere(R,x,y,z)
handles.fig = figure(20200225);
handles.subplot(1) = subplot(1,2,1);
[x,y,z]=sphere;
handles.sphere = surf(x,y,z);
colormap([0.9  0.9 0.9]);
hold on;
handles.scatter(1) = scatter3(0,0,0,30,[0 0 1],'filled');
handles.scatter(2) = scatter3(0,0,0,30,[1 0 0],'filled');
end