%======================================
clear all; clc
addpath Datasets
addpath Solvers_Table
addpath Subfunctions
%======================================
load('a9a.mat')
%======================================
X     = sparse(X);
[n,d] = size(X);

ka    = 1e-4;       % kappa regularizer parameter
X0    = [rand(d,1) ones(d,1)];

maxaver = 5;
Max_grad_num = 30*n;
Max_prox_num = 1e4;
Max_subgradient_num = 1e4;
eps = 1e-3;

fid = fopen('C:\Users\Kaitu\Desktop\Phd Papers\Paper 3-stochastic methods_phd3\phd3_SDCA\SDCA\Code_PhD3\Results\Tables\Result_PCA_a9a.txt','w');
fprintf(fid,'Algorithm     Fun,        Gap,       Num_Grad,     Num_Prox,     Time_iter \r\n');
for xlength = 1:2
    x0    = X0(:,xlength);
    for  Gfunindex  = 1:5
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
        
        FunVRSPA=0; GapVRSPA=0; Time_stopVRSPA=0; Time_iterVRSPA=0; Num_GradVRSPA=0;Num_ProxVRSPA=0;
        FunSPDCA=0; GapSPDCA=0; Time_stopSPDCA=0; Time_iterSPDCA=0; Num_GradSPDCA=0; Num_ProxSPDCA=0;
        FunSSPDCA=0; GapSSPDCA=0; Time_stopSSPDCA=0; Time_iterSSPDCA=0; Num_GradSSPDCA=0; Num_SubfSSPDCA=0;
        FunSSDCSPG=0; GapSSDCSPG=0; Time_stopSSDCSPG=0; Time_iterSSDCSPG=0; Num_GradSSDCSPG=0; Num_SubfSSDCSPG=0;
        
        for aver = 1:maxaver
            % VRSPA SPDCA SSPDCA SSDCSPG
            [SSPDCAoutput, Num_Grad_sspdca,Num_subgrad_sspdca] = ...
                SSPDCA_Table(n,X,x0,Max_grad_num,Max_subgradient_num,subgradient_fun_f2,kalambda);
            [VRSPAoutput, lambda_vrspa, Num_Grad_vrspa] = ...
                VRSPA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox);
            [SPDCAoutput,  lambda_spdca,Num_Grad_spdca] = ...
                SPDCA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox);
            b1 = 100;
            [SSDCSPGoutput, Num_Grad_ssdcspg, Num_subgrad_ssdcspg]  = ...
                SSDCSPG_Table(n,X,x0,b1,subgradient_fun_f1,subgradient_fun_f2,Max_grad_num,Max_subgradient_num);
            
            TimVRSPA = VRSPAoutput(d+1,:);
            [funVRSPA1, gapVRSPA1,time_stopVRSPA1,time_iterVRSPA1,num_GradVRSPA1,num_ProxVRSPA1] = ...
                Err_Fun_Gap_Table(VRSPAoutput(1:d,:),TimVRSPA,Num_Grad_vrspa, X,lambda_vrspa,Prox,Mygfun,'VRSPA',eps);
            FunVRSPA  = FunVRSPA  + funVRSPA1;
            GapVRSPA  = GapVRSPA  + gapVRSPA1;
            Time_stopVRSPA = Time_stopVRSPA + time_stopVRSPA1;
            Time_iterVRSPA = Time_iterVRSPA + time_iterVRSPA1;
            Num_GradVRSPA = Num_GradVRSPA + num_GradVRSPA1;
            Num_ProxVRSPA = Num_ProxVRSPA + num_ProxVRSPA1;
            
            TimSPDCA = SPDCAoutput(d+1,:);
            [funSPDCA1, gapSPDCA1,time_stopSPDCA1,time_iterSPDCA1,num_GradSPDCA1,num_ProxSPDCA1] =...
                Err_Fun_Gap_Table(SPDCAoutput(1:d,:),TimSPDCA, Num_Grad_spdca,X,lambda_spdca,Prox,Mygfun,'SPDCA',eps);
            FunSPDCA  = FunSPDCA + funSPDCA1;
            GapSPDCA  = GapSPDCA + gapSPDCA1;
            Time_stopSPDCA = Time_stopSPDCA + time_stopSPDCA1;
            Time_iterSPDCA = Time_iterSPDCA + time_iterSPDCA1;
            Num_GradSPDCA  = Num_GradSPDCA  + num_GradSPDCA1;
            Num_ProxSPDCA  = Num_ProxSPDCA  + num_ProxSPDCA1;
            
            
            TimSSPDCA = SSPDCAoutput(d+1,:);
            [funSSPDCA1, gapSSPDCA1,time_stopSSPDCA1,time_iterSSPDCA1,num_GradSSPDCA1,num_SubfSSPDCA1] =...
                Err_Fun_Gap_Table_SS(SSPDCAoutput(1:d,:),TimSSPDCA, Num_Grad_sspdca,Num_subgrad_sspdca,X,lambda_spdca,Prox,Mygfun,'SSPDCA',eps);
            FunSSPDCA       = FunSSPDCA + funSSPDCA1;
            GapSSPDCA       = GapSSPDCA + gapSSPDCA1;
            Time_stopSSPDCA = Time_stopSSPDCA + time_stopSSPDCA1;
            Time_iterSSPDCA = Time_iterSSPDCA + time_iterSSPDCA1;
            Num_GradSSPDCA  = Num_GradSSPDCA  + num_GradSSPDCA1;
            Num_SubfSSPDCA  = Num_SubfSSPDCA  + num_SubfSSPDCA1;
            
            
            TimSSDCSPG = SSDCSPGoutput(d+1,:);
            [funSSDCSPG1, gapSSDCSPG1,time_stopSSDCSPG1,time_iterSSDCSPG1,num_GradSSDCSPG1,num_SubfSSDCSPG1] =...
                Err_Fun_Gap_Table_SS(SSDCSPGoutput(1:d,:),TimSSDCSPG, Num_Grad_ssdcspg,Num_subgrad_ssdcspg,X,lambda_spdca,Prox,Mygfun,'SSPDCSPG',eps);
            FunSSDCSPG       = FunSSDCSPG + funSSDCSPG1;
            GapSSDCSPG       = GapSSDCSPG + gapSSDCSPG1;
            Time_stopSSDCSPG = Time_stopSSDCSPG + time_stopSSDCSPG1;
            Time_iterSSDCSPG = Time_iterSSDCSPG + time_iterSSDCSPG1;
            Num_GradSSDCSPG  = Num_GradSSDCSPG  + num_GradSSDCSPG1;
            Num_SubfSSDCSPG  = Num_SubfSSDCSPG  + num_SubfSSDCSPG1;
        end
        fprintf(fid,'-----x0=%d------ Gfunindex = %s ==================\r\n',x0(1), name_pena);
        fprintf(fid,'VRSPA       & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f \r\n',  FunVRSPA/maxaver, GapVRSPA/maxaver, (Num_GradVRSPA/maxaver)/n,  Num_ProxVRSPA/maxaver,  Time_iterVRSPA/maxaver);
        fprintf(fid,'SPDCA       & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f \r\n', FunSPDCA/maxaver, GapSPDCA/maxaver, (Num_GradSPDCA/maxaver)/n,   Num_ProxSPDCA/maxaver, Time_iterSPDCA/maxaver);
        fprintf(fid,'SSPDCA      & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f \r\n', FunSSPDCA/maxaver, GapSSPDCA/maxaver, (Num_GradSSPDCA/maxaver)/n,   Num_SubfSSPDCA/maxaver, Time_iterSSPDCA/maxaver);
        fprintf(fid,'SSDCSPG     & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f \r\n', FunSSDCSPG/maxaver,GapSSDCSPG/maxaver, (Num_GradSSDCSPG/maxaver)/n,   Num_SubfSSDCSPG/maxaver, Time_iterSSDCSPG/maxaver);
    end
end