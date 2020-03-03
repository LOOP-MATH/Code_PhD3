clear all; clc
close all
%addpath C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures
%% data set a9a
for i = 1:5
    figure(i)
    if i == 1
        load('CappL1_Alla9a_SS.mat')
    elseif i == 2
        load('L1L2_Alla9a_SS.mat')
    elseif i == 3
        load('LSP_Alla9a_SS.mat')
    elseif  i == 4
        load('MCP_Alla9a_SS.mat')
    else
        load('SCAD_Alla9a_SS.mat')
    end
    semilogy(TimVRSPA,FunVRSPA,'-r','LineWidth', 2);
    hold on
    semilogy(TimSPDCA,FunSPDCA,'-b','LineWidth', 2);
    hold on
    semilogy(TimSSDCSPG,FunSSDCSPG,'+g','LineWidth', 2);
    hold on
    semilogy(TimSSPDCA,FunSSPDCA,'*y','LineWidth', 2);
    
    xh = xlabel( 'time (s)', 'FontSize', 16, 'FontWeight', 'bold' );
    set(xh,'Interpreter','latex');
    yh = ylabel( '$f(\omega_k)\|$', 'FontSize', 14, 'FontWeight', 'bold' );
    set(yh,'Interpreter','latex');
    
    h = legend( { 'VRSPA','SPIDER-Moreau-PDCA', 'SSDCSPG', 'SPIDER-Inexact-PDCA' },'location','NorthEast', 'FontSize',16, 'FontWeight', 'bold' );
    set(h,'Interpreter','latex');
    set( gca, 'FontSize',  16, 'FontWeight', 'bold','linewidth',1.08 ) ;
end

% data set mnist

for i = 1:5
    figure(i+5)
    if i == 1
        load('CappL1_Allmnist_SS.mat')
    elseif i == 2
        load('L1L2_Allmnist_SS.mat')
    elseif i == 3
        load('LSP_Allmnist_SS.mat')
    elseif  i == 4
        load('MCP_Allmnist_SS.mat')
    else
        load('SCAD_Allmnist_SS.mat')
    end
    semilogy(TimVRSPA,FunVRSPA,'-r','LineWidth', 2);
    hold on
    semilogy(TimSPDCA,FunSPDCA,'-b','LineWidth', 2);
    hold on
    semilogy(TimSSDCSPG,FunSSDCSPG,'+g','LineWidth', 2);
    hold on
    semilogy(TimSSPDCA,FunSSPDCA,'*y','LineWidth', 2);
    
    xh = xlabel( 'time (s)', 'FontSize', 16, 'FontWeight', 'bold' );
    set(xh,'Interpreter','latex');
    yh = ylabel( '$f(\omega_k)\|$', 'FontSize', 14, 'FontWeight', 'bold' );
    set(yh,'Interpreter','latex');
    
    h = legend( { 'VRSPA','SPIDER-Moreau-PDCA', 'SSDCSPG', 'SPIDER-Inexact-PDCA' },'location','NorthEast', 'FontSize',16, 'FontWeight', 'bold' );
    set(h,'Interpreter','latex');
    set( gca, 'FontSize',  16, 'FontWeight', 'bold','linewidth',1.08 ) ;
end