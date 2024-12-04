% Parameters
L = 15;  % Length of the potential well
N = 1000;  % Number of discretization points
du = L / N;  % Step size
ulist = linspace(1e-5, L, N + 1)';  % Start from a small positive value to avoid division by zero

% Kinetic energy matrix (finite difference method)
e = ones(N + 1, 1);
KE = -spdiags([e -2*e e], -1:1, N + 1, N + 1) / du^2;

HOeigenvalues = [];
ISWeigenvalueslist = [];
FSWeigenvalueslist = [];
FSWSOCeigenvalueslist = [];

FSWeigenfunctions = [];

for l = 0:5 

    % Construct the Harmonic Oscillator PE matrix
    PE = spdiags(ulist.^2 + l * (l + 1) ./ ulist.^2, 0, N + 1, N + 1);
    H = KE + PE;

    % Find the eigenvalues and eigenvectors for Harmonic Oscillator
    [V, D] = eigs(H, 5, 'smallestabs');
    eigenvalues = diag(D);

    if l == 0 | l == 1 | l == 2 | l==3 | l == 4
        HOeigenvalues = [HOeigenvalues eigenvalues(1)];
    end

    % PE matrix for infinite square well
    a = 3.105;  % Length of the well (u < a means V = 0, u >= a means V = v_0)
    v_0 = 1000;  % Potential outside the well

    V_u = zeros(N + 1, 1);
    V_u(ulist >= a) = v_0;
    PEISW = spdiags(V_u + l * (l + 1) ./ ulist.^2, 0, N + 1, N + 1);
    
    ISWH = KE + PEISW;
     % Find the eigenvalues and eigenvectors for Infinite Square Well
    [V_ISW, D_ISW] = eigs(ISWH, 5, 'smallestabs');

    ISWeigenvalues = diag(D_ISW);
    
    if l == 0 
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1:3)'];
    elseif l == 3 | l == 4 | l == 5
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1)'];
    elseif l == 1 | 2
        ISWeigenvalueslist = [ISWeigenvalueslist ISWeigenvalues(1:2)'];
    end


    % PE matrix for Finite Square Well
    a = 3.105;  % Length of the well (u < a means V = 0, u >= a means V = v_0)
    v_0 = 30;  % Potential outside the well

    V_u = zeros(N + 1, 1);
    V_u(ulist >= a) = v_0;
    PEFSW = spdiags(V_u + l * (l + 1) ./ ulist.^2, 0, N + 1, N + 1);
    
    FSWH = KE + PEFSW;
     % Find the eigenvalues and eigenvectors for Finite Square Well
    [V_FSW, D_FSW] = eigs(FSWH, 5, 'smallestabs');

    FSWeigenvalues = diag(D_FSW);
    
    if l == 0 
        FSWeigenvalueslist = [FSWeigenvalueslist FSWeigenvalues(1:3)'];
        FSWeigenfunctions = [FSWeigenfunctions V_FSW(:, 1:3)];
    elseif l == 3 | l == 4 | l == 5
        FSWeigenvalueslist = [FSWeigenvalueslist FSWeigenvalues(1)'];
        FSWeigenfunctions = [FSWeigenfunctions V_FSW(:, 1)];
    elseif l == 1 | 2
        FSWeigenvalueslist = [FSWeigenvalueslist FSWeigenvalues(1:2)'];
        FSWeigenfunctions = [FSWeigenfunctions V_FSW(:, 1:2)];
    end

    % PE matrix for finite square well with SOC
    U = 0.1;

    PEminus = spdiags(V_u + l * (l + 1) ./ ulist.^2 - U.*l, 0, N + 1, N + 1);
    HSOCminus = KE + PEminus;

    PEplus = spdiags(V_u + l * (l + 1) ./ ulist.^2 - U.*(-l-1), 0, N + 1, N + 1);
    HSOCplus = KE + PEplus;
    
     % Find the eigenvalues and eigenvectors 
    [V_FSWplus, D_FSWplus] = eigs(HSOCplus, 5, 'smallestabs');
    [V_FSWminus, D_FSWminus] = eigs(HSOCminus, 5, 'smallestabs');

    FSWpluseigenvalues = diag(D_FSWplus);
    FSWminuseigenvalues = diag(D_FSWminus);
    
    if l == 0 
        FSWSOCeigenvalueslist = [FSWSOCeigenvalueslist FSWminuseigenvalues(1:3)'];
    elseif l == 3 | l == 4
        FSWSOCeigenvalueslist = [FSWSOCeigenvalueslist FSWpluseigenvalues(1)' FSWminuseigenvalues(1)'];
    elseif l == 1 | l == 2 
        FSWSOCeigenvalueslist = [FSWSOCeigenvalueslist FSWpluseigenvalues(1:2)' FSWminuseigenvalues(1:2)'];
    elseif l == 5
        FSWSOCeigenvalueslist = [FSWSOCeigenvalueslist FSWminuseigenvalues(1)'];
    end

end



disp(HOeigenvalues)
disp(ISWeigenvalueslist)
disp(FSWeigenvalueslist)
disp(FSWSOCeigenvalueslist)

%% Plot the eigenvalues
% Create a new figure
figure;

% Hold on to allow multiple lines to be plotted
hold on;

scaledHO = scale(HOeigenvalues)
scaledISW = scale(ISWeigenvalueslist)
scaledFSW = scale(FSWeigenvalueslist)
scaledFSWSOC = scale(FSWSOCeigenvalueslist)

% Loop through the array and plot each horizontal line for the scaled HOeigenvalues
for i = 1:length(scaledHO)
    % Plot a horizontal line at the height specified by the scaled array element
    plot([0, 4], [scaledHO(i), scaledHO(i)], 'r'); % 'r' sets the line color to red
end

% Loop through the array and plot each horizontal line for ISWeigenvalueslist
for i = 1:length(scaledISW)
    % Plot a horizontal line at the height specified by the array element
    plot([6, 10], [scaledISW(i), scaledISW(i)], 'b'); % 'b' sets the line color to blue
end

% Loop through the array and plot each horizontal line for FSWeigenvalueslist
for i = 1:length(scaledFSW)
    % Plot a horizontal line at the height specified by the array element
    plot([12, 16], [scaledFSW(i), scaledFSW(i)], 'g'); 
end

for i = 1:length(scaledFSWSOC)
    % Plot a horizontal line at the height specified by the array element
    plot([18, 22], [scaledFSWSOC(i), scaledFSWSOC(i)], 'm'); 
end

% Release the hold on the current figure
hold off;

% Set the axis limits for better visualization
xlim([0, 22]);
ylim([0, 8.05]);

% Remove the x-axis tick numbers
set(gca, 'XTick', []);

% Label the axes
ylabel('Energy levels');

% Add a title
title('Energy levels of different central potentials');

% Add custom text labels for specific sections
text(2, 0, {'Harmonic', 'Oscillator'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(8, 0, 'Infinite Square Well', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(14, 0, 'Finite Square Well', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(20, 0, {'Finite Square Well', 'with spin-orbit coupling'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

%% Plot the eigenfunctions for the finite square well
numCols = size(FSWeigenfunctions, 2);

colors = colormap(hsv(numCols)); 

orderedeigenvalues = sort(scaledFSW);

% Determine the number of rows and columns for subplots
numRows = ceil(sqrt(numCols));
numColsPerRow = ceil(numCols / numRows);

figure; % Create a single figure for all subplots

% Define the margins and spacing
margin = 0.05;
subplotWidth = (1 - (numColsPerRow + 1) * margin) / numColsPerRow;
subplotHeight = (1 - (numRows + 1) * margin) / numRows;

for i = 1:numCols
    row = ceil(i / numColsPerRow);
    col = mod(i - 1, numColsPerRow) + 1;
    left = margin + (col - 1) * (subplotWidth + margin);
    bottom = 1 - row * (subplotHeight + margin);
    
    % Create a subplot for each wave function
    axes('Position', [left, bottom, subplotWidth, subplotHeight]);
    plot(ulist, FSWeigenfunctions(:, i), 'LineWidth', 1, 'Color', colors(i, :)); % Add offset without scaling
    hold on; % Hold on to add more elements to the plot
    
    % Add vertical dotted line at ulist = 3.105
    [~, idx] = min(abs(ulist - 3.105)); % Find the closest index to 3.105
    xline(ulist(idx), '--', 'Color', 'k', 'LineWidth', 1.5); % Add the vertical line

    % Add horizontal line at y = offset(i)
    yline(0, '-', 'Color', 'k', 'LineWidth', 1.5); 
    set(gca, 'YTick', 0); 
    
    xlim([0, 4]);
    xlabel('$\rho$', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel('$u(\rho)$', 'Interpreter', 'latex', 'FontSize', 10);

    % Add the corresponding eigenvalue to the legend
    legend_text = [num2str(orderedeigenvalues(i))];
    h = legend(legend_text, 'Location', 'northeastoutside', 'FontSize', 8);
    title(h, 'Scaled eigenvalue'); % Add legend title
end

% Add a global title for the entire figure using annotation
annotation('textbox', [0 0.9 1 0.1], 'String', 'Finite Square Well Wavefunctions', ...
    'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center', 'EdgeColor', 'none');

