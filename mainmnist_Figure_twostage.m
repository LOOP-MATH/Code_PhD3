%======================================
clear all; clc
addpath Datasets
addpath Solvers_Figure
addpath Subfunctions
%======================================
load('mnist.mat')
%======================================
X              = sparse(X);
[n,d]          = size(X);
Max_nonsmooth_subgradient_num = 1e4;
Max_prox_num          = 1e4;
Max_smooth_grad_num   = 30*n;

ka    = 1e-4;       % kappa regularizer parameter
x0    = ones(d,1); % initial point setting
for  Gfunindex  = 1:5
    if Gfunindex == 1
        name_pena  = 'CappedL1';
        theta  = 0.01;
        lambda_capp = 1;
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_capp);
        subgradient_fun_f2 = @(x) subCapp(x, ka, lambda_capp, theta);
        Mygfun = @(wk) ka*sum(Fun_CappL1(wk,lambda_capp,theta));
        Prox   = @(wk,lambda) Prox_CappedL1(wk,ka,lambda,lambda_capp,theta);
    elseif  Gfunindex == 2
        name_pena  = 'L1l2';
        theta  = 1;
        lambda_l1l2  = 1;
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_l1l2);
        subgradient_fun_f2 = @(x) subl1l2(x, ka, lambda_l1l2, theta);
        Mygfun = @(wk)ka*Fun_L1L2(wk,lambda_l1l2,theta);
        Prox   = @(wk,lambda) Prox_L1L2(wk, ka,lambda,lambda_l1l2,theta);
    elseif  Gfunindex == 3
        theta = 2.1;
        lambda_scad = 1;
        name_pena  = 'SCAD';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_scad);
        subgradient_fun_f2 = @(x) subscad(x, ka, lambda_scad, theta);
        Mygfun = @(wk) ka*sum(Fun_SCAD(wk,lambda_scad,theta));
        Prox   = @(wk,lambda) Prox_SCAD(wk, ka,lambda,lambda_scad,theta);
    elseif Gfunindex == 4
        theta = 2;
        lambda_mcp = 1;
        name_pena  = 'MCP';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_mcp);
        subgradient_fun_f2 = @(x) submcp(x, ka, lambda_mcp, theta);
        Mygfun = @(wk) ka*sum(Fun_MCP(wk,lambda_mcp,theta));
        Prox   = @(wk,lambda) Prox_MCP(wk, ka,lambda,lambda_mcp,theta);
    else
        theta = 1;
        lambda_lsp = 0.1;
        name_pena  = 'LSP';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_lsp);
        subgradient_fun_f2 = @(x) sublsp(x, ka, lambda_lsp, theta);
        Mygfun = @(wk) ka*sum(Fun_LSP(wk,lambda_lsp,theta));
        Prox   = @(wk,lambda) Prox_LSP(wk, ka,lambda,lambda_lsp,theta);
    end
   
    % SSDCSPG
    b1  = 100;
    output_SSDCSPG  =   SSDCSPG_Figure(n,X,x0,b1,subgradient_fun_f1,...
        subgradient_fun_f2,Max_smooth_grad_num,Max_nonsmooth_subgradient_num);
    FunSSDCSPG  = Err_Fun_Gap_SS(output_SSDCSPG(1:d,:), X, Mygfun,'SSDCSPG');
    TimSSDCSPG  = output_SSDCSPG(d+1,:);

    % VRSPA
    [output_VRSPA, lambda_vrspa]  = VRSPA_Figure(n,X,x0,Max_smooth_grad_num,Max_prox_num,Prox);
    FunVRSPA = Err_Fun_Gap_SS(output_VRSPA(1:d,:), X, Mygfun,'VRSPA');
    TimVRSPA  = output_VRSPA(d+1,:);
    % SPDCA
    [output_SPDCA,  lambda_spdca] = SPDCA_Figure(n,X,x0,Max_smooth_grad_num,Max_prox_num,Prox);
    FunSPDCA   = Err_Fun_Gap_SS(output_SPDCA(1:d,:), X,Mygfun,'SPDCA');
    TimSPDCA   = output_SPDCA(d+1,:);
    
    % SSDCPDCA
    output_SSDPDCA  = SSPDCA_Figure(n,X,x0,subgradient_fun_f1,...
        subgradient_fun_f2,Max_smooth_grad_num,Max_nonsmooth_subgradient_num);
    FunSSPDCA   = Err_Fun_Gap_SS(output_SSDPDCA(1:d,:), X,Mygfun,'SSPDCA');
    TimSSPDCA   = output_SSDPDCA(d+1,:);
    
    
    clear output_SSDCSPG output_SSDPDCA output_VRSPA  output_SPDCA
    if Gfunindex == 1
        save('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures\CappL1_Allmnist_SS.mat','TimSSDCSPG','TimSSPDCA', 'FunSSDCSPG','FunSSPDCA','TimVRSPA','TimSPDCA','FunVRSPA','FunSPDCA');
    elseif Gfunindex == 2
        save('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures\L1L2_Allmnist_SS.mat','TimSSDCSPG','TimSSPDCA', 'FunSSDCSPG','FunSSPDCA','TimVRSPA','TimSPDCA','FunVRSPA','FunSPDCA');
    elseif Gfunindex == 3
        save('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures\SCAD_Allmnist_SS.mat','TimSSDCSPG','TimSSPDCA', 'FunSSDCSPG','FunSSPDCA','TimVRSPA','TimSPDCA','FunVRSPA','FunSPDCA');
    elseif Gfunindex == 4
        save('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures\MCP_Allmnist_SS.mat','TimSSDCSPG','TimSSPDCA', 'FunSSDCSPG','FunSSPDCA','TimVRSPA','TimSPDCA','FunVRSPA','FunSPDCA');
    else
        save('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Figures\LSP_Allmnist_SS.mat','TimSSDCSPG','TimSSPDCA', 'FunSSDCSPG','FunSSPDCA','TimVRSPA','TimSPDCA','FunVRSPA','FunSPDCA');
    end
end