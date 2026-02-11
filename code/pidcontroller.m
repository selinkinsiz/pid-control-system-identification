clear; clc; close all;

fprintf('============================================================\n');
fprintf('MECA311 PROJECT\n');
fprintf('============================================================\n');

%% ===================== Q1: LOAD DATA =====================
fprintf('\n=== Q1: DATA ACQUISITION ===\n');
load('measurement_time_output.mat');   % measurement data: variables time and output_response

t = time(:);             % convert time to a column vector
y = output_response(:);  % convert output to a column vector

fprintf('Loaded N = %d samples\n', numel(t));
fprintf('t range  = %.4f to %.4f s\n', t(1), t(end));

figure(1); clf;
plot(t, y, 'b', 'LineWidth', 1.6); grid on;
xlabel('Time (s)'); ylabel('Output');
title('Measured Open-Loop Step Response');
hold on;

%% ===================== Q2–Q3: IDENTIFICATION =====================
fprintf('\n=== Q2–Q3: SYSTEM IDENTIFICATION (zeta, wn) ===\n');

N = numel(y);
tailN = max(100, round(0.1*N));              % last 10% (at least 100 samples) for steady-state estimate
x_ss = mean(y(end-tailN+1:end));             % steady-state estimate from tail average

[y_peak, idxp] = max(y);                     % peak value and its index
Tp = t(idxp);                                % peak time

OS_pct = 100*(y_peak - x_ss)/x_ss;           % percent overshoot
OS_dec = OS_pct/100;                         % overshoot ratio for formula

if ~isfinite(OS_dec) || OS_dec <= 0
    error('Invalid overshoot ratio. Check data / x_ss estimate.');
end

% standard 2nd-order relationships: OS -> zeta, Tp -> wn
zeta = -log(OS_dec)/sqrt(pi^2 + (log(OS_dec))^2);
wn   = pi/(Tp*sqrt(1 - zeta^2));

fprintf('x_ss   = %.6f\n', x_ss);
fprintf('y_peak = %.6f\n', y_peak);
fprintf('Tp     = %.6f s\n', Tp);
fprintf('OS     = %.2f %%\n', OS_pct);
fprintf('zeta   = %.6f\n', zeta);
fprintf('wn     = %.6f rad/s\n', wn);

plot(Tp, y_peak, 'ro', 'MarkerFaceColor','r');
yline(x_ss,'g--','x_{ss}');
legend('Measured','Peak','Steady-state','Location','best');

%% ===================== Q4: k and c (m=1) =====================
fprintf('\n=== Q4: PHYSICAL PARAMETERS (m=1 kg) ===\n');

m = 1;                   % mass given in the project
k = m*wn^2;              % wn^2 = k/m -> k = m*wn^2
c = 2*zeta*wn*m;         % 2*zeta*wn = c/m -> c = 2*zeta*wn*m

fprintf('m = %.2f kg\n', m);
fprintf('k = %.6f N/m\n', k);
fprintf('c = %.6f Ns/m\n', c);

G = tf(1,[m c k]);       % plant: Force -> displacement

%% ===================== Q5: Time-domain specs from measured data (manual) =====================
fprintf('\n=== Q5: TIME-DOMAIN SPECS (MEASURED DATA, MANUAL) ===\n');

Ts_meas = settling_time_2pct(t, y, x_ss);    % 2% settling time from measured data (manual band method)
Ts_approx = 4/(zeta*wn);                     % theoretical approximation for Ts (2% criterion)

% also show theory steady-state for a unit step force: x_ss = 1/k
xss_theory = 1/k;

fprintf('Measured Tp                 = %.6f s\n', Tp);
fprintf('Measured OS                 = %.2f %%\n', OS_pct);
fprintf('Measured Ts (2%% band)       = %.6f s\n', Ts_meas);
fprintf('Approx Ts = 4/(zeta*wn)      = %.6f s\n', Ts_approx);
fprintf('Theory steady-state (1/k)   = %.6f\n', xss_theory);

%% ===================== Q6: PID CONTROLLER =====================
fprintf('\n=== Q6: PID CONTROLLER DESIGN (GAINS GIVEN) ===\n');

% PID gains (use your working values)
Kp = 190.0;
Ki = 1450.0;
Kd = 25.0;
Tf = 0.0037;             % derivative filter constant

C = pid(Kp, Ki, Kd, Tf); % PID controller (filtered derivative)

L = C*G;                 % loop transfer function
Tref = feedback(L, 1);   % closed-loop: reference -> output

% --- manual extraction of Tp/OS/Ts from closed-loop step response ---
t_step = (0:0.001:2).';                 % time vector for step analysis
[y_cl_step, t_cl_step] = step(Tref, t_step);

yss_cl = mean(y_cl_step(end-200:end));  % steady-state estimate from last samples
[ymax_cl, idx_cl] = max(y_cl_step);
Tp_cl = t_cl_step(idx_cl);
OS_cl = 100*(ymax_cl - yss_cl)/yss_cl;
Ts_cl = settling_time_2pct(t_cl_step, y_cl_step, yss_cl);

fprintf('PID gains:\n');
fprintf('  Kp=%.3f, Ki=%.3f, Kd=%.3f, Tf=%.6f\n', Kp, Ki, Kd, Tf);

fprintf('Closed-loop (Tref) step metrics (MANUAL from step response):\n');
fprintf('  Tp = %.6f s\n', Tp_cl);
fprintf('  OS = %.3f %%\n', OS_cl);
fprintf('  Ts = %.6f s\n', Ts_cl);
fprintf('  y_ss (estimated) = %.6f\n', yss_cl);

figure(2); clf;
plot(t_cl_step, y_cl_step, 'r', 'LineWidth', 1.6); grid on;
yline(yss_cl,'k--','y_{ss}');
title('Closed-Loop Step Response (T_{ref})');
xlabel('Time (s)'); ylabel('Output');

%% ===================== ROOT LOCUS PLOT (SAFE LEGEND) =====================
fprintf('\n=== ROOT LOCUS PLOT ===\n');

% target pole marker for visualization only
Tp_target = 0.145;
OS_target = 10;
zeta_d = -log(OS_target/100)/sqrt(pi^2 + (log(OS_target/100))^2);
wn_d   = pi/(Tp_target*sqrt(1-zeta_d^2));
sigma_d = zeta_d*wn_d;
wd_d    = wn_d*sqrt(1-zeta_d^2);
s_target = -sigma_d + 1i*wd_d;

p_cl = pole(Tref);        % actual closed-loop poles

figure(3); clf;
rlocus(L); grid on; hold on;
plot(real(s_target),  imag(s_target), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(s_target), -imag(s_target), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(p_cl), imag(p_cl), 'ko', 'MarkerSize', 7, 'LineWidth', 1.5);

title('Root Locus of L(s)=C(s)G(s) with Target & Actual Poles');
xlim([-80 10]); ylim([-50 50]);

% dummy handles to avoid legend issues with rlocus objects
hRL = plot(nan,nan,'b-');
hT  = plot(nan,nan,'rx','MarkerSize',10,'LineWidth',2);
hP  = plot(nan,nan,'ko','MarkerSize',7,'LineWidth',1.5);
legend([hRL hT hP], {'Root Locus','Target Poles','Actual CL Poles'}, 'Location','best');

%% ===================== Q7: OPEN-LOOP CUTOFF =====================
fprintf('\n=== Q7: OPEN-LOOP CUTOFF FREQUENCY ===\n');

w = logspace(-4, 6, 60000);             % wide sweep so cutoff is not missed
[magG, phG] = bode(G, w);
magG = squeeze(magG);
phG  = squeeze(phG);

DC_G = 1/k;                              % manual DC gain for G(s)=1/(ms^2+cs+k)
targetG = DC_G/sqrt(2);                  % -3 dB relative to DC

[wc_ol, fc_ol, ph_wc_ol] = cutoff_from_mag_robust(w, magG, phG, targetG);

fprintf('Open-loop DC gain |G(0)|   = %.8f\n', DC_G);
fprintf('Open-loop cutoff wc_ol     = %.6f rad/s\n', wc_ol);
fprintf('Open-loop cutoff fc_ol     = %.6f Hz\n', fc_ol);
fprintf('Phase at wc_ol             = %.2f deg\n', ph_wc_ol);

figure(4); clf;
bode(G); grid on;
title('Bode Plot – Open Loop G(s)');

%% ===================== Q8: CLOSED-LOOP CUTOFF (Tref) =====================
fprintf('\n=== Q8: CLOSED-LOOP CUTOFF FREQUENCY (Tref) ===\n');

[magT, phT] = bode(Tref, w);
magT = squeeze(magT);
phT  = squeeze(phT);

DC_T = 1;                                % with integral action, Tref(0) ~ 1 (manual assumption for tracking)
targetT = DC_T/sqrt(2);                  % -3 dB relative to DC

[wc_cl, fc_cl, ph_wc_cl] = cutoff_from_mag_robust(w, magT, phT, targetT);

if ~isfinite(fc_cl) || fc_cl<=0
    error('Q8: cutoff not found for Tref. Increase sweep or re-check cutoff definition.');
end

fprintf('Closed-loop DC gain |T(0)| = %.8f\n', DC_T);
fprintf('Closed-loop cutoff wc_cl   = %.6f rad/s\n', wc_cl);
fprintf('Closed-loop cutoff fc_cl   = %.6f Hz\n', fc_cl);
fprintf('Phase at wc_cl             = %.2f deg\n', ph_wc_cl);

figure(5); clf;
bode(Tref); grid on;
title('Bode Plot – Closed Loop T_{ref}(s)');

%% ===================== Q9: COMPARISON =====================
fprintf('\n=== Q9: OPEN vs CLOSED LOOP COMPARISON ===\n');

fprintf('Open-loop  fc_ol = %.6f Hz\n', fc_ol);
fprintf('Closed-loop fc_cl = %.6f Hz\n', fc_cl);

if isfinite(fc_ol) && isfinite(fc_cl) && fc_ol>0
    fprintf('Bandwidth ratio fc_cl/fc_ol = %.3f x\n', fc_cl/fc_ol);
end

%% ===================== Q10: GAIN & PHASE MARGINS =====================
fprintf('\n=== Q10: GAIN & PHASE MARGINS (Loop L=C*G) ===\n');

% note: margin() is a toolbox calculation; keep only if your instructor allows it
[gm, pm, wcg, wcp] = margin(L);

if isinf(gm)
    gm_dB = Inf;
else
    gm_dB = 20*log10(gm);
end

fprintf('Phase Margin (PM)  = %.2f deg\n', pm);
fprintf('Gain Margin  (GM)  = %s dB\n', ternary(isinf(gm_dB),'Inf',num2str(gm_dB)));
fprintf('Gain crossover wcp = %.6f rad/s\n', wcp);
fprintf('Phase crossover wcg= %.6f rad/s\n', wcg);

figure(6); clf;
margin(L); grid on;
title('Margin Plot – Loop L(s)=C(s)G(s)');

%% ===================== Q11: SINUSOIDAL INPUT ANALYSIS =====================
fprintf('\n=== Q11: SINUSOIDAL INPUT ANALYSIS ===\n');

fprintf('u(t) = 5 + sin(2*pi*f*t)\n');
fprintf('Open-loop input -> plant force (G)\n');
fprintf('Closed-loop input -> reference r(t) (Tref)\n');

if ~isfinite(fc_ol) || fc_ol<=0, error('Q11: fc_ol invalid. Check Q7.'); end
if ~isfinite(fc_cl) || fc_cl<=0, error('Q11: fc_cl invalid. Check Q8.'); end

t_sin = (0:0.001:10).';     % time vector for sinusoidal tests (column)

u1     = 5 + sin(2*pi*1*t_sin);
u_fcOl = 5 + sin(2*pi*fc_ol*t_sin);
r1     = 5 + sin(2*pi*1*t_sin);
r_fcCl = 5 + sin(2*pi*fc_cl*t_sin);

assert(all(isfinite(u1)) && all(isfinite(u_fcOl)) && all(isfinite(r1)) && all(isfinite(r_fcCl)));

y_ol_1Hz = lsim(G,    u1,     t_sin);
y_ol_fc  = lsim(G,    u_fcOl, t_sin);
y_cl_1Hz = lsim(Tref, r1,     t_sin);
y_cl_fc  = lsim(Tref, r_fcCl, t_sin);

amp = @(sig) (max(sig(end-2000:end)) - min(sig(end-2000:end)))/2;

fprintf('Amplitude (steady-state ~ last 2s):\n');
fprintf('  OL @ 1 Hz     : %.6f\n', amp(y_ol_1Hz));
fprintf('  OL @ fc_ol    : %.6f\n', amp(y_ol_fc));
fprintf('  CL @ 1 Hz     : %.6f\n', amp(y_cl_1Hz));
fprintf('  CL @ fc_cl    : %.6f\n', amp(y_cl_fc));

figure(7); clf;

subplot(2,2,1);
plot(t_sin, u1, 'k--', t_sin, y_ol_1Hz, 'b', 'LineWidth', 1.2); grid on;
title('Open Loop: 1 Hz'); xlabel('t (s)'); ylabel('y');
legend('u(t)','y(t)','Location','best'); xlim([8 10]);

subplot(2,2,2);
plot(t_sin, u_fcOl, 'k--', t_sin, y_ol_fc, 'b', 'LineWidth', 1.2); grid on;
title(sprintf('Open Loop: f_{c,ol}=%.3f Hz', fc_ol)); xlabel('t (s)'); ylabel('y');
legend('u(t)','y(t)','Location','best'); xlim([8 10]);

subplot(2,2,3);
plot(t_sin, r1, 'k--', t_sin, y_cl_1Hz, 'r', 'LineWidth', 1.2); grid on;
title('Closed Loop (Tref): 1 Hz'); xlabel('t (s)'); ylabel('y');
legend('r(t)','y(t)','Location','best'); xlim([8 10]);

subplot(2,2,4);
plot(t_sin, r_fcCl, 'k--', t_sin, y_cl_fc, 'r', 'LineWidth', 1.2); grid on;
title(sprintf('Closed Loop (Tref): f_{c,cl}=%.3f Hz', fc_cl)); xlabel('t (s)'); ylabel('y');
legend('r(t)','y(t)','Location','best'); xlim([8 10]);

%% ===================== FINAL CONSOLE SUMMARY =====================
fprintf('\n============================================================\n');
fprintf('--- FINAL CONSOLE SUMMARY ---\n');
fprintf('============================================================\n');

fprintf('\n--- SYSTEM IDENTIFICATION ---\n');
fprintf('wn  = %.3f rad/s\n', wn);
fprintf('zeta = %.3f\n', zeta);
fprintf('k = %.4f N/m\n', k);
fprintf('c = %.4f Ns/m\n', c);
fprintf('Tp = %.4f s\n', Tp);
fprintf('Overshoot = %.2f %%\n', OS_pct);
fprintf('Ts (measured,2%%) = %.4f s\n', Ts_meas);

fprintf('\n--- PID CONTROL RESULTS (MANUAL from step response) ---\n');
fprintf('Kp=%.2f  Ki=%.2f  Kd=%.2f  Tf=%.5f\n', Kp, Ki, Kd, Tf);
fprintf('Tp = %.4f s\n', Tp_cl);
fprintf('Overshoot = %.2f %%\n', OS_cl);
fprintf('Ts = %.4f s\n', Ts_cl);

fprintf('\n--- FREQUENCY DOMAIN RESULTS ---\n');
fprintf('Open Loop Cutoff Frequency   = %.3f rad/s\n', wc_ol);
fprintf('Closed Loop Cutoff Frequency = %.3f rad/s\n', wc_cl);
fprintf('Phase Margin = %.2f deg\n', pm);
fprintf('Gain Margin  = %s dB\n', ternary(isinf(gm_dB),'Inf',num2str(gm_dB)));

fprintf('\nDONE.\n');

%% ===================== LOCAL FUNCTIONS =====================
function Ts = settling_time_2pct(tvec, yvec, yss)
% 2% settling time from data:
% Ts is taken as the first time after the last sample outside the ±2% band around yss
    if abs(yss) < 1e-12
        Ts = NaN; return;
    end
    band = 0.02*abs(yss);
    last_out = find(abs(yvec - yss) > band, 1, 'last');
    if isempty(last_out)
        Ts = tvec(1);
    elseif last_out < numel(tvec)
        Ts = tvec(last_out + 1);
    else
        Ts = NaN;
    end
end

function [wc, fc, ph_wc] = cutoff_from_mag_robust(w, mag, ph, targetMag)
% robust cutoff finder (relative to DC):
% - find resonance peak
% - search after the peak for first crossing below targetMag
% - interpolate on log(w) for smoother estimate
    wc = NaN; fc = NaN; ph_wc = NaN;

    if any(~isfinite(mag)) || any(~isfinite(w))
        return;
    end

    [~, imax] = max(mag);
    startIdx = max(2, imax);

    idx_rel = find(mag(startIdx:end) <= targetMag, 1, 'first');
    if isempty(idx_rel)
        return;
    end
    idx = startIdx + idx_rel - 1;

    if idx <= 1
        return;
    end

    w1 = w(idx-1); w2 = w(idx);
    m1 = mag(idx-1); m2 = mag(idx);
    p1 = ph(idx-1);  p2 = ph(idx);

    % if target is not bracketed, fall back to the nearest sample
    if (m1-targetMag)*(m2-targetMag) > 0
        wc = w(idx);
        fc = wc/(2*pi);
        ph_wc = ph(idx);
        return;
    end

    logw = interp1([m1 m2], log10([w1 w2]), targetMag, 'linear', 'extrap');
    wc = 10^(logw);
    fc = wc/(2*pi);
    ph_wc = interp1(log10([w1 w2]), [p1 p2], log10(wc), 'linear', 'extrap');
end

function out = ternary(cond, a, b)
% small helper: returns a if cond is true, otherwise returns b
    if cond, out = a; else, out = b; end
end

