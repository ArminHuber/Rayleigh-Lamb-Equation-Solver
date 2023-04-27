%% Rayleigh-Lamb Equation Solver
% Update 2023/04/27
% written by: Armin Huber, armin.huber@dlr.de
% -------------------------------------------------------------------------
% Calculate Lamb wave dispersion diagrams for a free isotropic plate by
% evaluating the Rayleigh-Lamb equations. The Rayleigh-Lamb equation
% amplitudes are scaled to a grayscale. Modal solutions are represented by
% minima, i.e., by dark shading.

%% settings
Material.Name = 'aluminum';
Material.LongitudinalVelocity = 6320; % (m/s)
Material.TransverseVelocity = 3130; % (m/s)

Thickness = 1; % plate thickness (mm)

PhaseVelocityLimit = 20; % (m/ms)
FrequencyLimit = 10000; % (kHz)
PhaseVelocitySteps = 1000;
FrequencySteps = 1000;

%% computation
Half = Thickness/2e3;
Frequency = 0:FrequencyLimit/FrequencySteps:FrequencyLimit; % generate a range of frequencies
PhaseVelocity = (0:PhaseVelocityLimit/PhaseVelocitySteps:PhaseVelocityLimit)*1e3; % generate a range of phase velocities
AngularFrequency = 2*pi*Frequency*1e3;
k2 = (AngularFrequency./PhaseVelocity').^2; % wavenumber^2 of Lamb waves
kL2 = (AngularFrequency/Material.LongitudinalVelocity).^2; % wavenumber^2 of longitudinal bulk waves
kT2 = (AngularFrequency/Material.TransverseVelocity).^2; % wavenumber^2 of transverse bulk waves
x = sqrt(kL2-k2); % out-of-plane wavenumber component of longitudinal bulk waves
y = sqrt(kT2-k2); % out-of-plane wavenumber component of transverse bulk waves
a1 = (y.^2-k2).^2./(4*k2.*x.*y);
a2 = tan(x*Half)./tan(y*Half);
S = abs(a1+a2); % Rayleigh-Lamb equation for symmetric modes (absolute value)
A = abs(a1+1./a2); % Rayleigh-Lamb equation for antisymmetric modes (absolute value)

%% dispersion diagrams
f = figure('Name','Symmetric mode dispersion diagram','Units','normalized','Toolbar','none','OuterPosition',[0 0 1 1],'color','w');
datacursormode on
imagesc(Frequency,PhaseVelocity/1e3,20*log10(S)) % Rayleigh-Lamb equation amplitude in dB
colormap(gray(1024));
ax = gca;
ax.FontSize = 24;
ax.Title.FontSize = 30;
ax.XLabel.FontSize = 30;
ax.YLabel.FontSize = 30;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.Title.String = ['Symmetric mode dispersion diagram of ',num2str(Thickness),'\,mm ',char(join(split(Material.Name,'_'),'\_'))];
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.YDir = 'normal';
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @MyCursor;

f = figure('Name','Antisymmetric mode dispersion diagram','Units','normalized','Toolbar','none','OuterPosition',[0 0 1 1],'color','w');
datacursormode on
imagesc(Frequency,PhaseVelocity/1e3,20*log10(A)) % Rayleigh-Lamb equation amplitude in dB
colormap(gray(1024));
ax = gca;
ax.FontSize = 24;
ax.Title.FontSize = 30;
ax.XLabel.FontSize = 30;
ax.YLabel.FontSize = 30;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.Title.String = ['Antisymmetric mode dispersion diagram of ',num2str(Thickness),'\,mm ',char(join(split(Material.Name,'_'),'\_'))];
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.YDir = 'normal';
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @MyCursor;

function output_txt = MyCursor(~,event_obj)
    output_txt = {['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
end