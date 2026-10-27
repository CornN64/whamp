Te = 100; % electron temperature (eV)
ne = 1; % number density (cm^-3)
B = 50; % magnetic field strength (nT)

Units = irf_units;
me = Units.me;
qe = Units.e;
cspeed = Units.c;
epso=Units.eps0;

wce = qe*B*1e-9/me;
wpe = sqrt(ne*1e6*qe^2/me/epso);
fce = wce/(2*pi);
veth = sqrt(2*qe*Te/me);


% Guess starting frequency with cold plasma prediction
kparstart = 0.01;
kparend = 0.01*100;
fstart = (kparstart*wce/veth)^2*cspeed^2*wce/(wpe^2)/wce*1.0; % cold plasma prediction for whistler with normalization
%fstart = wpe/wce*1.0+0.5;

Electrons = struct('m',0,'n',ne,'t',Te,'a',5,'vd',0,'b',0,'d',1);
PlasmaModel = struct('B',B);
PlasmaModel.Species = {Electrons};

InputParameters = struct('fstart',fstart,'kperp',[kparstart 0.01 kparend],'kpar',[kparstart 0.01 kparend],'useLog',0); %Whistler mode
tic
Output = mexwrap(PlasmaModel,InputParameters);
toc

%%
kpar = Output.kpar*wce/veth;
kperp = Output.kperp*wce/veth;
freq = Output.f*wce;


%freq(real(freq) > fstart*wce*5) = NaN;
%freq(imag(freq) < -100) = NaN;
%freq(imag(freq) > 100) = NaN;

Ex = Output.Ex;
Ey = Output.Ey;
Ez = Output.Ez;

Bx = Output.Bx;
By = Output.By;
Bz = Output.Bz;

Emag = sqrt(Ex.*conj(Ex)+Ey.*conj(Ey)+Ez.*conj(Ez));
Bmag = sqrt(Bx.*conj(Bx)+By.*conj(By)+Bz.*conj(Bz));

Eperp=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey));
Epolar=real(-2*imag(Ex.*conj(Ey))./(Eperp.*Eperp));

EparoverE = real(sqrt(Ez.*conj(Ez))./Emag);
cBoverE = real(cspeed*Bmag./Emag*1e-6);

% Remove bad points
freq(Output.flagSolutionFound == 0) = NaN;
freq(Output.flagTooHeavilyDamped == 1) = NaN;
freq(Output.flagNoConvergence == 1) = NaN;
freq(real(freq) < 0) = NaN;


%%

fn=figure;
set(fn,'Position',[10 10 900 700])
    h(1)=axes('position',[0.07 0.55 0.41 0.39]); % [x y dx dy]
    h(2)=axes('position',[0.56 0.55 0.41 0.39]);
    h(3)=axes('position',[0.07 0.05 0.41 0.39]);
    h(4)=axes('position',[0.56 0.05 0.41 0.39]);
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

gammamax = max(max((imag(freq)/wce)));
h(1)=irf_panel('disprel');
surf(h(1),kperp,kpar,real(freq)/wce,imag(freq)/wce,'LineStyle','none')
set(h(1),'ylim',[0 kpar(end-1)]);
set(h(1),'xlim',[0 kperp(end-1)]);
caxis(h(1),[-gammamax gammamax])
ylabel(h(1),'k_{||} (m^{-1})','fontsize',14)
xlabel(h(1),'k_{\perp} (m^{-1})','fontsize',14)
zlabel(h(1),'\omega/\Omega_{ce}','fontsize',14)
title(h(1),'Growth rate: \gamma/\Omega_{ce}','fontsize',14)
set(h(1),'fontsize',14);
c=colorbar('peer',h(1),'ver');

h(2)=irf_panel('pol');
surf(h(2),kperp,kpar,real(freq)/wce,Epolar,'LineStyle','none')
set(h(2),'ylim',[0 kpar(end-1)]);
set(h(2),'xlim',[0 kperp(end-1)]);
caxis(h(2),[-1 1])
ylabel(h(2),'k_{||} (m^{-1})','fontsize',14)
xlabel(h(2),'k_{\perp} (m^{-1})','fontsize',14)
zlabel(h(2),'\omega/\Omega_{ce}','fontsize',14)
title(h(2),'Polarization','fontsize',14)
set(h(2),'fontsize',14);
c=colorbar('peer',h(2),'ver');

h(3)=irf_panel('EparoverE');
surf(h(3),kperp,kpar,real(freq)/wce,EparoverE,'LineStyle','none')
set(h(3),'ylim',[0 kpar(end-1)]);
set(h(3),'xlim',[0 kperp(end-1)]);
caxis(h(3),[0 1])
ylabel(h(3),'k_{||} (m^{-1})','fontsize',14)
xlabel(h(3),'k_{\perp} (m^{-1})','fontsize',14)
zlabel(h(3),'\omega/\Omega_{ce}','fontsize',14)
title(h(3),'|E_{||}|/|E|','fontsize',14)
set(h(3),'fontsize',14);
c=colorbar('peer',h(3),'ver');

h(4)=irf_panel('cBoverE');
surf(h(4),kperp,kpar,real(freq)/wce,log10(cBoverE),'LineStyle','none')
set(h(4),'ylim',[0 kpar(end-1)]);
set(h(4),'xlim',[0 kperp(end-1)]);
caxis(h(4),[0 2])
ylabel(h(4),'k_{||} (m^{-1})','fontsize',14)
xlabel(h(4),'k_{\perp} (m^{-1})','fontsize',14)
zlabel(h(4),'\omega/\Omega_{ce}','fontsize',14)
title(h(4),'c|B|/|E|','fontsize',14)
set(h(4),'fontsize',14);
c=colorbar('peer',h(4),'ver');

rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

colormap(h(1),bgrcmap);
colormap(h(2),bgrcmap);
colormap(h(3),'jet');
colormap(h(4),'jet');

set(gcf,'color','w');
%axis(h(1),[0 4e-3 0 4e-3])
