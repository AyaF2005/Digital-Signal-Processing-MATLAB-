%% ============================
%  Q2: Minimum-Phase & All-Pass
%  ============================

clear; clc; close all;
disp('Designing and analyzing systems for Q2...');

%% -------------------------------------------------------
%  2.1 Minimum-Phase System (FIR) – guaranteed min-phase
%     Approach: design a linear-phase FIR; convert to minimum-phase
% --------------------------------------------------------
Nlin   = 40;             % linear-phase prototype order (even)
Wc     = 0.30;           % normalized cutoff (× pi rad/sample)
h_lin  = fir1(Nlin, Wc, hamming(Nlin+1), 'noscale');   % linear-phase FIR

% Convert to minimum-phase (zeros strictly inside unit circle)
if exist('firminphase','file') == 2
    b_min = firminphase(h_lin);
else
    % Fallback if firminphase is unavailable:
    % rceps decomposes a *sequence* into min/max-phase components.
    % For our filter prototype it works well in practice.
    [b_min, ~] = rceps(h_lin);
end
a_min = 1;

% Sanity checks: zeros inside the unit circle
z_min = roots(b_min);
fprintf('Min-phase FIR: max |zero radius| = %.4f (should be < 1)\n', max(abs(z_min)));

% Open FVTool for the minimum-phase filter
hMin = fvtool(b_min, a_min, 'Name', 'Minimum-Phase FIR (use Analysis menu)');
% In FVTool: View -> Analysis -> Magnitude / Phase / Group Delay / Pole-Zero

% Also show Pole-Zero with zplane (optional helper)
figure('Name','Min-Phase FIR: Pole-Zero (zplane)'); zplane(b_min, a_min); grid on;

% Print difference equation y[n] = sum_k b_min(k+1) x[n-k]
fprintf('\nDifference equation (Min-Phase FIR):\n');
fprintf('y[n] ='); 
for k = 0:numel(b_min)-1
    coeff = b_min(k+1);
    if k==0, fprintf(' %.6g x[n]', coeff);
    else,    fprintf(' + %.6g x[n-%d]', coeff, k);
    end
end
fprintf('\n');

%% ---------------------------------------
%  2.2 All-Pass System (second-order IIR)
%  ---------------------------------------
r   = 0.7;               % pole radius (<1 for stability)
th  = pi/3;              % pole angle
p   = r*exp(1j*th);      % one complex pole
a_ap = poly([p, conj(p)]);       % denominator: 1, -2r cos(th), r^2
b_ap = fliplr(conj(a_ap));       % numerator: reversed & conj to make all-pass
% For real-coefficient all-pass, conj(.) is redundant since a_ap is real.

% Verify unity magnitude numerically
[H_ap, w] = freqz(b_ap, a_ap, 2048);
mag_dev = max(abs(abs(H_ap)-1));
fprintf('All-pass IIR: max |H| deviation from 1 over grid = %.2e\n', mag_dev);

% Open FVTool for the all-pass filter
hAP = fvtool(b_ap, a_ap, 'Name', 'All-Pass IIR (use Analysis menu)');
% In FVTool: View -> Analysis -> Magnitude / Phase / Group Delay / Pole-Zero

% Also show Pole-Zero with zplane
figure('Name','All-Pass IIR: Pole-Zero (zplane)'); zplane(b_ap, a_ap); grid on;

% Print difference equation:
% y[n] - a1 y[n-1] - a2 y[n-2] = b0 x[n] + b1 x[n-1] + b2 x[n-2]
fprintf('\nDifference equation (All-Pass IIR):\n');
a1 = a_ap(2); a2 = a_ap(3);
b0 = b_ap(1); b1 = b_ap(2); b2 = b_ap(3);
fprintf('y[n] = %.6g x[n] + %.6g x[n-1] + %.6g x[n-2] + %.6g y[n-1] + %.6g y[n-2]\n', ...
        b0, b1, b2, a1, a2);

%% ------------------------------------------------------------
%  (Optional) Quick comparison plots (non-FVTool) for convenience
% ------------------------------------------------------------
[H_min, w] = freqz(b_min, a_min, 2048);
Gd_min = grpdelay(b_min, a_min, 2048);
Gd_ap  = grpdelay(b_ap , a_ap , 2048);

figure('Name','Quick Magnitude & Phase (Optional)');
subplot(2,1,1); plot(w/pi, 20*log10(abs(H_min)), 'LineWidth',1.1); grid on;
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Mag (dB)');
title('Minimum-Phase FIR Magnitude');

subplot(2,1,2); plot(w/pi, unwrap(angle(H_min)), 'LineWidth',1.1); grid on;
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Phase (rad)');
title('Minimum-Phase FIR Phase');

figure('Name','All-Pass Checks (Optional)');
subplot(2,1,1); plot(w/pi, 20*log10(abs(H_ap)), 'LineWidth',1.1); grid on;
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Mag (dB)');
title('All-Pass Magnitude (should be ~0 dB across)');

subplot(2,1,2); plot(w/pi, Gd_ap, 'LineWidth',1.1); grid on;
xlabel('Normalized Frequency (\times \pi rad/sample)'); ylabel('Group Delay (samples)');
title('All-Pass Group Delay');

disp('Script finished.');
