function theta = load_pulse_theta(Adjfile, pulsefile, ROI, scaleFA)
% Calculate effective theta for a given pulse
% requires mz blochsim code
% mz407

arguments
    Adjfile     string
    pulsefile   string
    ROI         logical
    scaleFA     double  = 1.0
end

sspulse = true; % for now we only use this for 2D slices

adj = extractFieldsforBloch(Adjfile);

nCha = 8; dt=1e-5;
pulse = load(pulsefile); % only .mat now
totalg = pulse.g_low * 1e-3;
totalRF = reshape(pulse.b_est, [], nCha);

if sspulse
    z = zeros(size(adj.z));
else
    z = adj.z;
end
[Mxy, Mz] = blochSim_CK_3D(totalRF*scaleFA, totalg, dt, adj.b0, ...
                           adj.x,adj.y,z, adj.b1, [0;0;1], 1e10, 1e10, 0);

phout = angle(squeeze(Mxy));
faout = atan2(squeeze(abs(Mxy)), squeeze(Mz)); % correct for [0 pi] OK for excitation
complexout = faout .* exp(1j*(phout-pi/2)); % -pi/2 to match axis in EPG

if nargin < 3
    ROI = logical(adj.W);
end

% Output: use one channel only as output
theta = zeros(nnz(ROI), 1, nCha);
theta(:,1,1) = complexout(ROI);
