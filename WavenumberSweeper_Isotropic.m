%% Wavenumber Sweeper Isotropic
% Update 2022/10/05
% written by: Armin Huber, armin.huber@dlr.de
% -------------------------------------------------------------------------
% Evaluate characteristic functions of Lamb and Scholte waves in
% fluid-loaded, viscoelastic, isotropic plates.

%% settings
clear
Materials = MaterialList;

FluidLoading = 1;
Fluid = Materials.Fluid.water;

Material = Materials.Isotropic.Aluminum_Disperse_Viscoelastic;
% Material = Materials.Isotropic.Aluminum_Disperse;
Thickness = 1; % (mm)

Frequency = 500; % (kHz)
ModeType = 'Lamb'; % Lamb/Scholte (necessary only when viscoelastic)
Symmetry = 'S'; % S: symmetric A: antisymmetric

SweepRangeReal = 5394:1:5444;
SweepRangeImag = -4:-1:-54;

%% calculation
AngularFrequency = 2*pi*Frequency*1e3;
WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
Half = Thickness/2e3;
k2 = WaveNumber.^2;
kL2 = (AngularFrequency/Material.LongitudinalVelocity_complex)^2;
kT2 = (AngularFrequency/Material.TransverseVelocity_complex)^2;
x = sqrt(kL2-k2);
y = sqrt(kT2-k2);
if  strcmp(Symmetry,'S')
    R = (y.^2-k2).^2./(y.*tan(x*Half))+4*k2.*x./tan(y*Half);
elseif strcmp(Symmetry,'A')
    R = (y.^2-k2).^2./y.*tan(x*Half)+4*k2.*x.*tan(y*Half);
end
if  FluidLoading
    kF2 = (AngularFrequency/Fluid.Velocity)^2;
    a = Fluid.Density*kT2^2*x./(y*Material.Density.*sqrt(kF2-k2));
    if  ~isreal(Material.C) && strcmp(ModeType,'Scholte')
        a = -a;
    end
    if  strcmp(Symmetry,'S')
        Y = abs(R-1i*a);
    elseif strcmp(Symmetry,'A')
        Y = abs(R+1i*a);
    end
else
    Y = abs(R); %#ok<*UNRCH> 
end
Min = zeros(size(Y,1),size(Y,2));
for l = 2:size(Y,2)-1
    for j = 2:size(Y,1)-1
        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
            Min(j,l) = 1;
        end
    end
end
[b1,b2] = find(Min);
if  ~isempty(b1)
    MIN = [b1 b2];
    disp(['cp:    ',num2str(AngularFrequency/real(WaveNumber(MIN(1),MIN(2)))),' m/s',newline...
          'alpha: ',num2str(imag(WaveNumber(MIN(1),MIN(2)))),' Np/m',newline...
          'creal: ',num2str(SweepRangeReal(MIN(1))),' m/s',newline...
          'cimag: ',num2str(SweepRangeImag(MIN(2))),' m/s',newline,'-------------------'])
end
figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))

%% high quality plot
% c = AngularFrequency./real(WaveNumber)/1e3;
% a = imag(WaveNumber);
% f = figure('Toolbar','none','Units','normalized','OuterPosition',[0 0 .75 1],'color','w');
% surf(a(floor(height(a)/2),:),c(:,floor(width(c)/2)),20*log10(Y))
% ax = gca;
% ax.Box = 'on';
% ax.LineWidth = 1;
% ax.XLabel.String = 'Attenuation (Np/m)';
% ax.YLabel.String = 'Phase velocity (m/ms)';
% ax.ZLabel.String = 'abs($G$) (dB)';
% ax.FontSize = 20;
% ax.XLabel.FontSize = 24;
% ax.YLabel.FontSize = 24;
% ax.XLabel.Interpreter = 'latex';
% ax.YLabel.Interpreter = 'latex';
% ax.ZLabel.Interpreter = 'latex';
% ax.TickLabelInterpreter = 'latex';
% tb = axtoolbar('default');
% tb.Visible = 'on';
% % exportgraphics(f,fullfile('C:\Users\Public\Armin\MATLAB','ComplexFunctionAmplitude.png'),'Resolution',150)