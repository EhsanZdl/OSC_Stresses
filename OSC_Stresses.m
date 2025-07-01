% OSC_Stresses.m
% ------------------------------------------------------------
% This script extract yield stress (τ_y) and flow-point stress (τ_f) from 
% TA-TRIOS Oscillatory amplitude sweep test excel extract file (*.xls) that 
% sits in same folder of this file and generate a plot for each file.
% Outputs are saved in a “result” sub-folder created on the same address of 
% parent folder.
%
% QUICK SETUP  ────────────────────────────────────────────────────────────
%  1. Copy this script into the same directory that contains your .xls
%  files
%     
%  2. Edit the “parentFolder” path below to match that directory.
% 
%  3. Run the script in MATLAB.  It will:
%        • compute τ_y and τ_f for every file,
%        • save a summary table  ➜  result/yield_flow_points.xlsx
%        • save one PNG plot     ➜  result/<filename>_plot.png
%
% METHOD SUMMARY
%  • τ_y  – two-tangent method in log-log space:
%           – fit a flat “plateau” to the elastic region
%           – fit the first part of the drop after the breakpoint
%           – the intersection gives yield stress
%  • τ_f  – stress where storage and loss moduli cross (G′ = G″),
%           found by sign change of (G″–G′) and log-space interpolation.
%
% Ehsan Zadehali (zadeh006@umn.edu)(eh.zadehali@gmail.com)
% Version1: 30 May 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc;

parentFolder = ''; 

resultFolder = fullfile(parentFolder,'result');
if ~exist(resultFolder,'dir') 
    mkdir(resultFolder); 
end

files  = dir(fullfile(parentFolder,'*.xls'));     
if isempty(files)
    error('No .xls files found in %s',parentFolder);
end

S = struct('File',{},'Yield',{},'Flow',{});
for k = 1:numel(files)
    F = files(k).name;
    M = readmatrix(F,'Sheet',2,'Range','A4:K1000');
    G  = M(:,1);           % Storage modulus column in excel file
    GG = M(:,2);           % loss modulus column in excel file
    tau_osc  = M(:,11);          % Oscillation stress column in excel file
    ok = ~isnan(G) & ~isnan(GG) & ~isnan(tau_osc);
    G = G(ok); GG = GG(ok); tau_osc = tau_osc(ok);

    %% --- Yield stress (two‑tangent method) -----------------------------------
    logX_tau_osc = log10(tau_osc); 
    logY_G = log10(G);
    min_first_points = 5; 
    N = numel(logY_G);

    % Initial break-----------------------------------------
    Best_Break = min_first_points; 
    bestErr = inf;
    for Break = min_first_points:(N-min_first_points)

        p1 = polyfit(logX_tau_osc(1:Break),logY_G(1:Break),1);
        y1 = polyval(p1, logX_tau_osc(1:Break));

        p2 = polyfit(logX_tau_osc(Break:end),  logY_G(Break:end),  1);
        y2 = polyval(p2, logX_tau_osc(Break:N));

        SSE = sum((logY_G(1:Break)-y1).^2)+sum((logY_G(Break:N)-y2).^2);

        if SSE < bestErr 
            bestErr = SSE; 
            Best_Break = Break; 
        end
    end

    region1 = 1:Best_Break;
    region2 = Best_Break:N;

    % Region2----------------------------------
    logX_reg2 = logX_tau_osc(region2);
    logY_reg2 = logY_G(region2);
    stress_reg2 = tau_osc(region2);

    window = 6;      
    min_points_reg2 = 4;   

    % Backward step in stress-----------------------------
    first_backward = find(diff(stress_reg2) <= 0, 1);
    if  ~isempty(first_backward) && first_backward > window
        last_point = first_backward;
    else
        last_point = min(window, numel(logX_reg2));
    end
    if last_point < min_points_reg2
       last_point = min(numel(logX_reg2), min_points_reg2);
    end
    logX_reg2_fit = logX_reg2(1:last_point);
    logY_reg2_fit = logY_reg2(1:last_point);

    % Final tangent fit region2------------------------------
    p_reg2 = polyfit(logX_reg2_fit, logY_reg2_fit, 1);
    slope_reg2 = p_reg2(1);
    intercept_reg2 = p_reg2(2);

    % Region1-----------------------------------------------
    plateauLevel = prctile(logY_G(region1),90);
    fullX   = logspace(log10(min(tau_osc)), log10(max(tau_osc)),200).';
    plateau = real(10.^plateauLevel*ones(size(fullX)));
    drop    = real(10.^(slope_reg2*log10(fullX)+intercept_reg2));
    tau_y  = real(10^((plateauLevel-intercept_reg2)/slope_reg2));

    rep.fullX      = fullX;
    rep.plateau    = plateau;
    rep.drop       = drop;
    rep.yield      = tau_y;
    rep.breakIndex = Best_Break;

    %% --- Flow‑point stress (G′ = G″) ----------------------------------
    diffG = GG - G;
    crosspint = find(diffG(1:end-1).*diffG(2:end) < 0, 1);
    if isempty(crosspint)
        warning('No G'' = G" crossover in %s',F); tau_f = NaN;
    else
        x1 = log10(tau_osc(crosspint));   
        x2 = log10(tau_osc(crosspint+1));
        y1 = diffG(crosspint);      
        y2 = diffG(crosspint+1);
        r  = abs(y1)/(abs(y1)+abs(y2));
        tau_f = 10^(x1 + r*(x2-x1));
    end

    %% --- Store results -------------------------------------------------
    S(k).File  = F;
    S(k).Yield = tau_y;
    S(k).Flow  = tau_f;
    %% --- Plot replicate ----------------------------------------------
    fullX = logspace(log10(min(tau_osc)), log10(max(tau_osc)), 200)';
    plateau = 10.^plateauLevel * ones(size(fullX));
    drop    = 10.^(slope_reg2*log10(fullX) + intercept_reg2);

    yl = [1e-2 1e5];
    fig = figure('Units','pixels','Position',[100 100 1000 500]);

    % Left panel – two‑tangent yield stress-----------------
    subplot(1,2,1); 
    hold on; 
    axis square; 
    box on; 
    grid off;
    set(gca,'XScale','log','YScale','log');

    plot(fullX, plateau,'k-','LineWidth',0.7);
    plot(fullX, drop,   'k-','LineWidth',0.7);
    plot(tau_osc, G,'or','MarkerSize',4,'MarkerFaceColor','none');
    plot([tau_y tau_y], yl,'k--','LineWidth',1);

    Gy = 10^plateauLevel;          
    plot(tau_y, Gy,'kx','MarkerSize',10,'LineWidth',1.5);

    text(0.05,0.1,sprintf('yield stress = %.1f Pa',tau_y),'Units','normalized');
    
    xlim([1e-2 1e5]); 
    ylim(yl);
    xlabel('\sigma_{osc} (Pa)'); 
    ylabel('G'' (Pa)');
    

    % Right panel – G′ & G″ crossover for flow point---------------------
    subplot(1,2,2); 
    hold on; 
    axis square; 
    box on; 
    grid off;
    set(gca,'XScale','log','YScale','log');

    plot(tau_osc, G, 'or','MarkerSize',4,'MarkerFaceColor','none','DisplayName','G''');
    plot(tau_osc, GG,'ob','MarkerSize',4,'MarkerFaceColor','none','DisplayName','G"');

    if ~isnan(tau_f)
        % approximate G'' at crossover for marker height-----------------
        [~, idxC] = min(abs(log10(tau_osc) - log10(tau_f)));
        Gc = G(idxC);
        plot([tau_f tau_f], yl,'k--','LineWidth',1);
        plot(tau_f, Gc,'rx','MarkerSize',10,'LineWidth',1.5);
    end

   text(0.05,0.1,sprintf('flow-point stress = %.1f Pa',tau_f),'Units','normalized');


    xlim([1e-2 1e5]); ylim(yl);
    xlabel('\sigma_{osc} (Pa)'); ylabel('Modulus (Pa)');
   

    saveas(fig, fullfile(resultFolder, sprintf('%s_plot.png', strrep(F,'.xls',''))));
    
end

%% --- Write summary table ----------------------------------------------
T = struct2table(S);
writetable(T,'yield_flow_points.xlsx');

fprintf('\n%-30s  %-10s  %-10s\n','File','τ_y (Pa)','τ_f (Pa)');
for k = 1:height(T)
    fprintf('%-30s  %10.3g  %10.3g\n',T.File{k},T.Yield(k),T.Flow(k));
end
fprintf('\nResults saved to yield_flow_points.xlsx and plots/ folder in %s\n',resultFolder);
