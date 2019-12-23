clear all; clc
close all
%addpath C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures
%% data set a9a
for i = 1:5
    figure(i)
    if i == 1
        load('CappL1_Alla9a.mat')
    elseif i == 2
        load('L1L2_Alla9a.mat')
    elseif i == 3
        load('LSP_Alla9a.mat')
    elseif  i == 4
        load('MCP_Alla9a.mat')
    else
        load('SCAD_Alla9a.mat')
    end
    semilogy(TimVRSPA,GapVRSPA,'-r','LineWidth', 2);
    hold on
    semilogy(TimSPDCA,GapSPDCA,'-b','LineWidth', 2);
    xh = xlabel( 'time (s)', 'FontSize', 16, 'FontWeight', 'bold' );
    set(xh,'Interpreter','latex');
    yh = ylabel( '$\|\mathcal{P}(\bar{\omega}_i, \frac{1}{\lambda_0}(\omega_i -\bar{\omega}_{i}), 1)\|$', 'FontSize', 14, 'FontWeight', 'bold' );
    set(yh,'Interpreter','latex');
    
    h = legend( { 'VRSPA','Algorithm 2'},'location','NorthEast', 'FontSize',16, 'FontWeight', 'bold' );
    set(h,'Interpreter','latex');
    set( gca, 'FontSize',  16, 'FontWeight', 'bold','linewidth',1.08 ) ;
end

% data set mnist

for i = 1:5
    figure(i+5)
    if i == 1
        load('CappL1_Allmnist.mat')
    elseif i == 2
        load('L1L2_Allmnist.mat')
    elseif i == 3
        load('LSP_Allmnist.mat')
    elseif  i == 4
        load('MCP_Allmnist.mat')
    else
        load('SCAD_Allmnist.mat')
    end
    semilogy(TimVRSPA,GapVRSPA,'-r','LineWidth', 2);
    hold on
    semilogy(TimSPDCA,GapSPDCA,'-b','LineWidth', 2);
    xh = xlabel( 'time (s)', 'FontSize', 16, 'FontWeight', 'bold' );
    set(xh,'Interpreter','latex');
    yh = ylabel( '$\|\mathcal{P}(\bar{\omega}_i, \frac{1}{\lambda_0}(\omega_i -\bar{\omega}_{i}), 1)\|$', 'FontSize', 14, 'FontWeight', 'bold' );
    set(yh,'Interpreter','latex');
    
    h = legend( { 'VRSPA','Algorithm 2'},'location','NorthEast', 'FontSize',16, 'FontWeight', 'bold' );
    set(h,'Interpreter','latex');
    set( gca, 'FontSize',  16, 'FontWeight', 'bold','linewidth',1.08 ) ;
end