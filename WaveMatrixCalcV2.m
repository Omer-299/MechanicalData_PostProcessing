clc
clear all
close all
data = readmatrix('Test5_take4.xlsx');

% figure out the region of interest - the portion with cyclic data
% time_meca = data(:,1);
% force = data(:,2);

time_meca = data(2:end,1);
force = data(2:end,2);

figure; plot(time_meca,force); grid on
title('Raw data: Force vs. Time')
xlabel('Time (s)')
ylabel('Force (N)')

% Separating high and low peaks
minProminence = 0.03 * max(force); % Minimum peak prominence
[highs, locs_highs] = findpeaks(force, 'MinPeakProminence', minProminence);
Time_highs = time_meca(locs_highs);
% highs_spacing = diff(locs_highs);

% figure; plot(diff(highs)); grid on
% title('Diff(F_{max}) - Difference b/w high peaks')
% ylabel('Force high peaks difference in each cycle (N)')
% ylim([0 15])

[lows, locs_lows] = findpeaks(-force, 'MinPeakProminence', minProminence);
Time_lows = time_meca(locs_lows);
lows_spacing = diff(locs_lows);

lows = -lows;

% figure; plot(diff(lows)); grid on
% title('Diff(lows) - Difference b/w low peaks')
% ylabel('Force low peaks difference in each cycle (N)')
% ylim([-10 20])

% plot Fhighs and Flows
figure;
plot(Time_highs,highs, '.r',MarkerSize=10);
title('Raw data: Force vs. Time')
hold on
plot(Time_lows, lows, '.b',MarkerSize=10)
plot(time_meca,force,'-k', LineWidth=0.5)
grid on
legend ('Fmax','Fmin', 'location','northwest')
xlabel('Time (s)')
ylabel('Force (N)')



% Combine and sort extrema
all_locs = [locs_highs; locs_lows; 1; length(force)];
all_vals = [highs; lows; force(1); force(end)];
is_max = [true(size(locs_highs)); false(size(locs_lows)); force(1) > force(2); force(end) > force(end-1)];
[all_locs_sorted, sort_idx] = sort(all_locs);
all_vals_sorted = all_vals(sort_idx);
is_max_sorted = is_max(sort_idx);

% Enforce zig-zag alternation
zigzag_idx = all_locs_sorted(1);
zigzag_vals = all_vals_sorted(1);
last_type = is_max_sorted(1);
for i = 2:length(all_locs_sorted)
    if is_max_sorted(i) ~= last_type
        zigzag_idx(end+1) = all_locs_sorted(i);
        zigzag_vals(end+1) = all_vals_sorted(i);
        last_type = is_max_sorted(i);
    end
end

F_filtered = force(zigzag_idx(2:end-2)); % starting from 2 because the first point is from the ramp step of the testing and it is not in the cyclic loaidng step.
T_filtered = time_meca(zigzag_idx(2:end-2));
T_filtered = T_filtered - T_filtered(1);

figure; plot (T_filtered,F_filtered); grid on
title('Filtered force data - only cyclic loading')
xlabel('Time (s)')
ylabel('Force (N)')

[Fhighs, locs_highs] = findpeaks(F_filtered);
Timehigh_filtered = T_filtered(locs_highs);
% highs_spacing = diff(locs_highs);

N = length(Timehigh_filtered);
x = round(0.8*N);
y = round(0.85*N);
% T_exp = mean(diff(Timehigh_filtered(end-500:end-200)));
T_exp = mean(diff(Timehigh_filtered(x:y))); % Loading time period (T)
LF_exp = 1/T_exp;
fprintf('Loading frequency as per the experimental data: %0.2f Hz\n',LF_exp);

[Flows, locs_lows] = findpeaks(-F_filtered);
Flows = -1 * Flows;
Timelows_filtered = T_filtered(locs_lows);
lows_spacing = diff(locs_lows);


figure;
plot(Timehigh_filtered,Fhighs, '.r');
title('Fmax and Fmin')
hold on
plot(Timelows_filtered, Flows, '.b')
grid on
legend ('Fmax','Fmin', 'location','northwest')
xlabel('Time (s)')
ylabel('Force (N)')


figure;
plot(Timehigh_filtered(2:end),diff(Fhighs), '.r');
title('\Delta F_{max} (N/cycle)')
hold on;
plot(Timehigh_filtered(2:end), movmean(diff(Fhighs), 30, "Endpoints", "fill"), '.k', 'markersize', 10);
grid on
xlabel('Time (s)')
ylabel('\Delta F_{max} (N/cycle)')
ylim([0 4])
legend('Instant \Delta F_{max}','Movmean \Delta F_{max}')


% Fhighs = Fhighs(1:length(Flows));Difference b/w high peaks
R = Flows./Fhighs(1:length(Flows)); 

figure; plot(Timelows_filtered(2:end-1),R(2:end-1), '-k', LineWidth=2.5); grid on    % (2:end-1) for eleminating outliers
title('Force ratio R')
xlabel('Time (s)')
ylabel('Force ratio R')
% ylim([-.2 .2])





%% Force ratio using linear regresion
% close all
window_fraction = 0.9;   % fraction of total time span
offset    = 0; %round(T_filtered(end)*0.0,0);       % offset in seconds, if needed

% Mid time
t_mid = (T_filtered(1) + T_filtered(end)) / 2 + offset;

% Half window in time
half_window_time = window_fraction * (T_filtered(end) - T_filtered(1)) / 2;

% Find indices of points within the window
start_idx = find(T_filtered >= t_mid - half_window_time, 1, 'first');
end_idx   = find(T_filtered <= t_mid + half_window_time, 1, 'last');


% Split alternating high/low values
if F_filtered(start_idx) > F_filtered(start_idx+1)
    F_high    = F_filtered(start_idx:2:end_idx);
    Time_high = T_filtered(start_idx:2:end_idx);
    F_low     = F_filtered(start_idx+1:2:end_idx);
    Time_low  = T_filtered(start_idx+1:2:end_idx);
else
    F_low     = F_filtered(start_idx:2:end_idx);
    Time_low  = T_filtered(start_idx:2:end_idx);
    F_high    = F_filtered(start_idx+1:2:end_idx);
    Time_high = T_filtered(start_idx+1:2:end_idx);
end

% Linear Regression Fits
f_high = fit(Time_high(:), F_high(:), 'poly1');
f_low  = fit(Time_low(:),  F_low(:),  'poly1');

% Evaluate fits
yfit_high = feval(f_high, Time_high);
yfit_low  = feval(f_low,  Time_low);

% Compute Force Ratio
Force_Ratio = f_low.p1 / f_high.p1;

% ---------------------- Create Figure ----------------------
figure('Position',[200 150 900 600]);

% --- Top Plot (Main Axes) ---
ax1 = axes('Position',[0.1 0.35 0.85 0.55]); % [x y width height]
hold(ax1,'on'); grid(ax1,'on');

% Plot original and selected data
plot(ax1, T_filtered, F_filtered, 'k.', 'DisplayName','Original Force');
plot(ax1, Time_high, F_high, 'ro', 'DisplayName','Force high (selected)');
plot(ax1, Time_low,  F_low,  'bo', 'DisplayName','Force low (selected)');

% Plot fits
plot(ax1, Time_high, yfit_high, 'g-', 'LineWidth',2.5, 'DisplayName','Fit high');
plot(ax1, Time_low,  yfit_low,  'm-', 'LineWidth',2,  'DisplayName','Fit low');

legend(ax1,'Location','northwest');
title(ax1,'Force vs. Time with Linear Fits; Test: GR250 Test5 take4');
xlabel(ax1,'Time (s)');
ylabel(ax1,'Force (N)');

% Mid-point line
mid_time_value = t_mid; %(min(T_filtered) + max(T_filtered)) / 2;
xline(ax1, mid_time_value, '--c', 'LineWidth', 1.5, 'HandleVisibility','off');
yl = ylim(ax1);
text(mid_time_value, yl(2), 'Mid time', 'Rotation', 90, ...
     'VerticalAlignment','bottom','HorizontalAlignment','right', ...
     'FontSize',9,'Color','c','Parent',ax1);

% ---------------------- Bottom Text Area ----------------------
ax2 = axes('Position',[0.1 0.05 0.85 0.25]); % lower area for text
axis(ax2,'off'); % hide axes

% Prepare text lines (formatted with LaTeX style)
eq1 = sprintf('$F_{high}(t) = %.2f\\,t + %.2f$', f_high.p1, f_high.p2);
eq2 = sprintf('$F_{low}(t) = %.2f\\,t + %.2f$',  f_low.p1,  f_low.p2);
Rtxt = sprintf('$R = F_{low}/F_{high} = %.3f$', Force_Ratio);
Wtxt = sprintf('Window fraction: $\\pm%.1f\\,\\%%$', window_fraction*100);
Otxt = sprintf('Mid-point offset: $%.2f\\,\\mathrm{s}$', offset);

% Combine in vertically aligned text
text(0.05, 0.75, eq1, 'Interpreter','latex', 'FontSize',14, 'Color','g', 'Parent',ax2);
text(0.05, 0.55, eq2, 'Interpreter','latex', 'FontSize',14, 'Color','m', 'Parent',ax2);
text(0.05, 0.35,  Rtxt, 'Interpreter','latex', 'FontSize',14, 'Color','k', 'Parent',ax2);
text(0.05, 0.15, Wtxt, 'Interpreter','latex', 'FontSize',11, 'Color','k', 'Parent',ax2);
text(0.05, 0.00, Otxt, 'Interpreter','latex', 'FontSize',11, 'Color','k', 'Parent',ax2);
