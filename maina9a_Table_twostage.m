% Readme

%======================================
clear all; clc
addpath Datasets
addpath Solvers_Table
addpath Subfunctions
%======================================
load('a9a.mat')
%======================================
X              = sparse(X);
[n,d]          = size(X);
Max_nonsmooth_subgradient_num = 1e4; %非光滑函数次梯度计算最大次数
Max_prox_num          = 1e4;         % 临近算子最大计算次数
Max_smooth_grad_num   = 30*n;        %  单一光滑函数梯度计算最大次数

ka    = 1e-4;       % kappa  平衡拟合项和正则项的正则参数
x0    = ones(d,1);  % 初始点
for  Gfunindex  = 5:5
    if Gfunindex == 1
        name_pena   = 'CappedL1';
        theta       = 0.01;
        lambda_capp = 1;
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_capp);
        subgradient_fun_f2 = @(x) subCapp(x, ka, lambda_capp, theta);
        Mygfun   = @(wk) ka*sum(Fun_CappL1(wk,lambda_capp,theta));
        Prox     = @(wk,lambda) Prox_CappedL1(wk,ka,lambda,lambda_capp,theta);
        kalambda = ka*lambda_capp;
    elseif  Gfunindex == 2
        name_pena  = 'L1l2';
        theta  = 1;
        lambda_l1l2  = 1;
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_l1l2);
        subgradient_fun_f2 = @(x) subl1l2(x, ka, lambda_l1l2, theta);
        Mygfun = @(wk) ka*Fun_L1L2(wk,lambda_l1l2,theta);
        Prox   = @(wk,lambda) Prox_L1L2(wk, ka,lambda,lambda_l1l2,theta);
        kalambda = ka*lambda_l1l2;
    elseif  Gfunindex == 3
        theta = 2.1;
        lambda_scad = 1;
        name_pena  = 'SCAD';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_scad);
        subgradient_fun_f2 = @(x) subscad(x, ka, lambda_scad, theta);
        Mygfun = @(wk) ka*sum(Fun_SCAD(wk,lambda_scad,theta));
        Prox   = @(wk,lambda) Prox_SCAD(wk, ka,lambda,lambda_scad,theta);
        kalambda = ka*lambda_scad;
    elseif Gfunindex == 4
        theta = 2;
        lambda_mcp = 1;
        name_pena  = 'MCP';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_mcp);
        subgradient_fun_f2 = @(x) submcp(x, ka, lambda_mcp, theta);
        Mygfun = @(wk) ka*sum(Fun_MCP(wk,lambda_mcp,theta));
        Prox   = @(wk,lambda) Prox_MCP(wk, ka,lambda,lambda_mcp,theta);
        kalambda = ka*lambda_mcp;
    else
        theta = 1;
        lambda_lsp = 0.1;
        name_pena  = 'LSP';
        subgradient_fun_f1 = @(x) subL1(x,ka,lambda_lsp);
        subgradient_fun_f2 = @(x) sublsp(x, ka, lambda_lsp, theta);
        Mygfun = @(wk) ka*sum(Fun_LSP(wk,lambda_lsp,theta));
        Prox   = @(wk,lambda) Prox_LSP(wk, ka,lambda,lambda_lsp,theta);
        kalambda = ka*lambda_lsp;
    end
    
    b1 = 100;  epsilon  = 1e-6;
    
%     output1 = SSDCSPG_Table_SS(n,X,x0,b1,subgradient_fun_f1,...
%         subgradient_fun_f2,Max_smooth_grad_num,Max_nonsmooth_subgradient_num,epsilon, Mygfun);
% %     output2 = VRSPA_Table_SS(n,X,x0,Max_smooth_grad_num,Max_prox_num,Prox,epsilon,Mygfun);
  %   output3 = SPDCA_Table_SS(n,X,x0,Max_smooth_grad_num,Max_prox_num,Prox,epsilon,Mygfun);
    
    epsilon  = 1e-16;
    output4 = SSPDCA_Table_SS(n,X,x0,...
        subgradient_fun_f2,Max_smooth_grad_num,Max_nonsmooth_subgradient_num, epsilon, Mygfun,kalambda);
    
    
    %output1
    %output2
    
   % output3
    output4
end