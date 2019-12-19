%======================================
clear all; clc
addpath Datasets
addpath Solvers_Figure
addpath Subfunctions
%======================================
load('a9a.mat')
%======================================
X             = sparse(X);
[n,d]         = size(X);
Max_prox_num  = 1e4;
Max_grad_num  = 30*n;

ka    = 1e-4;       % kappa regularizer parameter
x0    = ones(d,1); % initial point setting
for  Gfunindex  = 1:5
    if Gfunindex == 1
        name_pena  = 'CappedL1';
        theta  = 0.01;
        lambda_capp = 1;
        Prox   = @(wk,lambda) Prox_CappedL1(wk,ka,lambda,lambda_capp,theta);
        Mygfun = @(wk) ka*sum(Fun_CappL1(wk,lambda_capp,theta));
    elseif  Gfunindex == 2
        name_pena  = 'L1l2';
        theta  = 1;
        lambda_l1l2  = 1;
        Prox   = @(wk,lambda) Prox_L1L2(wk, ka,lambda,lambda_l1l2,theta);
        Mygfun = @(wk)ka*Fun_L1L2(wk,lambda_l1l2,theta);
    elseif  Gfunindex == 3
        theta = 2.1;
        lambda_scad = 1;
        name_pena  = 'SCAD';
        Prox   = @(wk,lambda) Prox_SCAD(wk, ka,lambda,lambda_scad,theta);
        Mygfun = @(wk) ka*sum(Fun_SCAD(wk,lambda_scad,theta));
    elseif Gfunindex == 4
        theta = 2;
        lambda_mcp = 1;
        name_pena  = 'MCP';
        Prox   = @(wk,lambda) Prox_MCP(wk, ka,lambda,lambda_mcp,theta);
        Mygfun = @(wk) ka*sum(Fun_MCP(wk,lambda_mcp,theta));
    else
        theta = 1;
        lambda_lsp = 0.1;
        name_pena  = 'LSP';
        Prox   = @(wk,lambda) Prox_LSP(wk, ka,lambda,lambda_lsp,theta);
        Mygfun = @(wk) ka*sum(Fun_LSP(wk,lambda_lsp,theta));
    end
    % VRSPA
    [VRSPAoutput, lambda_vrspa]  = VRSPA_Figure(n,X,x0,Max_grad_num,Max_prox_num,Prox);
    [SPDCAoutput,  lambda_spdca] = SPDCA_Figure(n,X,x0,Max_grad_num,Max_prox_num,Prox);
    
    [FunVRSPA, GapVRSPA]        = Err_Fun_Gap(VRSPAoutput(1:d,:), X,lambda_vrspa,Prox,Mygfun,'VRSPA');
    TimVRSPA = VRSPAoutput(d+1,:);
    [FunSPDCA, GapSPDCA] = Err_Fun_Gap(SPDCAoutput(1:d,:), X,lambda_spdca,Prox,Mygfun,'SPDCA');
    TimSPDCA = SPDCAoutput(d+1,:);
    clear SPDCAoutput VRSGAoutput
    if Gfunindex == 1
        save('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Figures\CappL1_Alla9a.mat','TimVRSPA','TimSPDCA', 'FunVRSPA','FunSPDCA','GapVRSPA','GapSPDCA');
    elseif Gfunindex == 2
        save('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Figures\L1L2_Alla9a.mat','TimVRSPA','TimSPDCA', 'FunVRSPA','FunSPDCA','GapVRSPA','GapSPDCA');
    elseif Gfunindex == 3
        save('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Figures\SCAD_Alla9a.mat','TimVRSPA','TimSPDCA', 'FunVRSPA','FunSPDCA','GapVRSPA','GapSPDCA');
    elseif Gfunindex == 4
        save('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Figures\MCP_Alla9a.mat','TimVRSPA','TimSPDCA', 'FunVRSPA','FunSPDCA','GapVRSPA','GapSPDCA');
    else
        save('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Figures\LSP_Alla9a.mat','TimVRSPA','TimSPDCA', 'FunVRSPA','FunSPDCA','GapVRSPA','GapSPDCA');
    end
end