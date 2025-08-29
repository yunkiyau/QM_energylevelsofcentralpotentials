% Parameters
L = 15;            % Length of the potential well
N = 1000;          % Number of discretization points
du = L / N;        % Step size
ulist = linspace(1e-5, L, N + 1)';  % Start from small positive value to avoid division by zero

% Kinetic energy matrix (finite difference method)
e = ones(N + 1, 1);
KE = -spdiags([e -2*e e], -1:1, N + 1, N + 1) / du^2;

HOeigenvalues = [];
ISWeigenvalueslist = [];

for l = 0:5
    % ----- Harmonic Oscillator -----
    % PE: rho^2 + l(l+1)/rho^2
    PE = spdiags(ulist.^2 + l * (l + 1) ./ ulist.^2, 0, N + 1, N + 1);
    H = KE + PE;

    % Find the eigenvalues (lowest 5 by abs value)
    [V, D] = eigs(H, 5, 'smallestabs');
    eigenvalues = diag(D);

    % Collect first (lowest) eigenvalue for selected l
    if l == 0 || l == 1 || l == 2 || l == 3 || l == 4
        HOeigenvalues = [HOeigenvalues eigenvalues(1)];
    end

    % ----- Infinite Square Well -----
    a = 3.105;    % Well width: V=0 for u<a, V=v0 otherwise
    v_0 = 1000;   % Very high outside potential (effectively infinite)

    V_u = zeros(N + 1, 1);
    V_u(ulist >= a) = v_0;
    PEISW = spdiags(V_u + l * (l + 1) ./ ulist.^2, 0, N + 1, N + 1);
    ISWH = KE + PEISW;

    [V_ISW, D_ISW] = eigs(ISWH, 5, 'smallestabs');
    ISWeigenvalues = diag(D_ISW);

    % Collect a few lowest levels by l
    if l == 0
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1:3)'];
    elseif l == 3 || l == 4 || l == 5
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1)'];
    elseif l == 1 || l == 2
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1:2)'];
    end
end

% Display results
disp(HOeigenvalues)
disp(ISWeigenvalueslist)

%% Plot the eigenvalues (HO vs ISW only)
figure; hold on;

% If you have a custom scale() function, keep these; otherwise remove/replace.
scaledHO  = scale(HOeigenvalues);
scaledISW = scale(ISWeigenvalueslist);

% Plot horizontal lines for HO
for i = 1:length(scaledHO)
    plot([0, 4], [scaledHO(i), scaledHO(i)], 'r'); % HO in red
end

% Plot horizontal lines for ISW
for i = 1:length(scaledISW)
    plot([6, 10], [scaledISW(i), scaledISW(i)], 'b'); % ISW in blue
end

hold off;

xlim([0, 10]);
ylim([0, 8.05]);
set(gca, 'XTick', []);
ylabel('Energy levels');
title('Energy levels: Harmonic Oscillator vs Infinite Square Well');
text(2, 0, {'Harmonic', 'Oscillator'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(8, 0, 'Infinite Square Well',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
