%BRILLOUIN ZONES
%VORONOI 3d
%Poms Johannes
%Based on QHULL DEMO
%http://www.mathworks.com/products/demos/shipping/matlab/qhulldemo.html
%Modified by Lebing Chen
%
%
%
% calculates the lattice first Brillouin zone given any spinW object 
% 
% ### Syntax
% 
% `fbz = sw_FBZ(obj,plt)`
% 
% ### Description
% 
% `fbz = sw_FBZ(obj,plt)` takes a spinW object with lattice parameters and 
% space group and converts it into a point group which can form a convex 
% hull representing the first Brillouin zone.
%
% ### Examples
% 
% test1 = spinw
% test1.genlattice('lat_const',[5 5 5],'angled',[90 90 90],'spgr',89)
% fbz = sw_FBZ(test1, true)
%
% ### Input Arguments
% 
% `obj`
% : spinW object with crystal structure and space group
% 
% `plt`
% : true or false, plot the FBZ or not, default is false
% 

function fbz = sw_FBZ(obj,varargin)


inpForm.fname  = {'plt' 'RtoP'};
inpForm.defval = {false false};
inpForm.size   = {1 1 };
inpForm.soft   = {false false};

param = sw_readparam(inpForm,varargin{:});
plt = param.plt;
R2P = param.RtoP;



sgroup = obj.lattice.label;

a = obj.lattice.lat_const(1);
b = obj.lattice.lat_const(2);
c = obj.lattice.lat_const(3);
alpha = obj.lattice.angle(1);
beta = obj.lattice.angle(2);
gamma = obj.lattice.angle(3);

if sgroup(1) == 'P'
    T1 = [a 0 0];
    T2 = [b*cos(gamma) b*sin(gamma) 0];
    T3 = [c*cos(beta) c/sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) c/sin(gamma)*sqrt(sin(gamma)^2-cos(alpha)^2-cos(beta)^2+2*cos(alpha)*cos(beta)*cos(gamma))];
elseif sgroup(1) == 'F'
    T1 =[0 b/2 c/2];
    T2 =[a/2 0 c/2];
    T3 =[a/2 b/2 0];
elseif sgroup(1) == 'I'
    T1 =[-a/2 b/2 c/2];
    T2 =[a/2 -b/2 c/2];
    T3 =[a/2 b/2 -c/2];
elseif sgroup(1) == 'A' || sgroup(1) == 'C'
    T1 =[a/2 b/2 0];
    T2 =[-a/2 b/2 0];
    T3 =[0 c*cos(alpha) c*sin(alpha)];
elseif sgroup(1) == 'R'
    if alpha == pi/2
        if R2P == false
            T1 =[-a/2 a/2/sqrt(3) c/3];
            T2 =[a/2 a/2/sqrt(3) c/3];
            T3 =[0 -a/sqrt(3) c/3];
            %T1 =[a/2*sqrt(3) -a/2 c/3];
            %T2 =[-a/2*sqrt(3) -a/2 c/3];
            %T3 =[0 a c/3];
        else
            T1 = [a 0 0];
            T2 = [b*cos(gamma) b*sin(gamma) 0];
            T3 = [c*cos(beta) c/sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) c/sin(gamma)*sqrt(sin(gamma)^2-cos(alpha)^2-cos(beta)^2+2*cos(alpha)*cos(beta)*cos(gamma))/3];
        end
    end
else
    T1 = [a 0 0];
    T2 = [b*cos(gamma) b*sin(gamma) 0];
    T3 = [c*cos(beta) c/sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) c/sin(gamma)*sqrt(sin(gamma)^2-cos(alpha)^2-cos(beta)^2+2*cos(alpha)*cos(beta)*cos(gamma))];
    
end
spat=cross(T1,T2)*T3';

G1=2*pi*cross(T2,T3)/spat;
G2=2*pi*cross(T3,T1)/spat;
G3=2*pi*cross(T1,T2)/spat;

X = zeros(125,3);

l=1;
for i=-2:2;
    for j=-2:2;
        for k=-2:2;
            X(l,:)=i*G1+j*G2+k*G3;
            l=l+1;
        end
    end
end


%cla reset; 
hold on;

% Compute Voronoi diagram.
[c,v] = voronoin(X);

nx = c(v{63},:);
tri = convhulln(nx);
fbz = nx;
% Plot the Voronoi diagram.
if plt == true
    fh=figure;
    plot3(X(:,1),X(:,2),X(:,3),'b.','markersize',10);
    for i = 1:size(tri,1)
        patch(nx(tri(i,:),1),nx(tri(i,:),2),nx(tri(i,:),3),i,'FaceAlpha',0.85);
    end
    % Modify the view.
    view(3);
    title('1st Brillouin zone');
    axis equal% tight off  vis3d
    grid on;
    camzoom(2);
    rotate3d on
end
%}
end
