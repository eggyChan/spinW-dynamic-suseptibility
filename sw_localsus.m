%Calculate \chi''(E) = (\chi_{xx}''(E) + \chi_{yy}''(E) + \chi_{zz}''(E))/3
function localsus = sw_localsus(obj,n,E,nE)


fbzhull = sw_FBZ(obj,false);
%Hmax = max(abs(fbzhull(:,1)));
%Kmax = max(abs(fbzhull(:,2)));
%Lmax = max(abs(fbzhull(:,3)));
Hmax = 1;
Kmax = 1;
Lmax = 1;
sel = 'Sxx+Syy+Szz';

rl = zeros(2*n+1,2*n+1,2*n);
Intensity = zeros(nE,1);
localsus = zeros(nE,2);
pts = 0;
wb=waitbar(0,'Calculating local suseptibility...');
for i = -n:n
    H = Hmax*i/n;
    for j = -n:n
        K = Kmax*j/n;
        [r,Spec] = evalc("obj.spinwave({[H, K, -Lmax], [H, K, Lmax], 2*n},'hermit',false,'formfact',false)");
        FBZHKLA = inhull(Spec.hklA',fbzhull);
        if sum(FBZHKLA)>0
            Spec = sw_neutron(Spec);
            Spec = sw_egrid(Spec, 'Evect',linspace(E/40,E,nE+1),'component',sel);
            [r,Spec] = evalc("sw_instrument(Spec,'dE',E/80,'norm',true)");
            rl(i+n+1,j+n+1,:) = FBZHKLA; 
            SinFBZ = Spec.swConv(:,FBZHKLA);
            Intensity = Intensity+sum(SinFBZ,2);
            pts = pts+sum(FBZHKLA);
        end
        waitbar((i+n+(j+n)/(2*n+1))/(2*n+1),wb,'Calculating local suseptibility...');
    end
end
waitbar(1,wb,'Calculation finished!');

Evec = E/40:E*39/40/(nE+1):E;
localsus(:,1) = Evec(2:nE+1);
localsus(:,2) = Intensity/pts*pi/3*13.77;
Ss = Intensity/pts*3.44;
S = sum(Ss)*(E*39/40)/nE/1000;
%close(wb);
disp('Local suseptibility calculation finished in unit MuB^2/eV/cell, assuming g=2.');
disp(['calculated nS = ',num2str(S)]);


    

