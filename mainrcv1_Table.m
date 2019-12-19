%======================================
clear all; clc
addpath Datasets
addpath Solvers_Table
addpath Subfunctions
%======================================
load('rcv1.mat')
%======================================
X     = sparse(X);
[n,d] = size(X);

ka    = 1e-4;       % kappa regularizer parameter
X0    = [rand(d,1) ones(d,1)];

maxaver = 5;
Max_grad_num = 30*n;
Max_prox_num = 1e4;
eps = 1e-3;

fid = fopen('C:\Users\Kaitu\Desktop\phd3_SDCA\SDCA\Results\Tables\Result_PCA_rcv1.txt','w');
fprintf(fid,'Algorithm     Fun,        Gap,       Num_Grad,     Num_Prox,     Time\r\n');
for xlength = 1:2
    x0    = X0(:,xlength);
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
        
        FunVRSPA=0; GapVRSPA=0; Time_stopVRSPA=0; Time_iterVRSPA=0; Num_GradVRSPA=0;Num_ProxVRSPA=0;
        FunSPDCA=0; GapSPDCA=0; Time_stopSPDCA=0; Time_iterSPDCA=0; Num_GradSPDCA=0; Num_ProxSPDCA=0;
        for aver = 1:maxaver
            % VRSPA SPDCA
            [VRSPAoutput, lambda_vrspa, Num_Grad_vrspa] = VRSPA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox);
            [SPDCAoutput,  lambda_spdca,Num_Grad_spdca] = SPDCA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox);
            
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
            
        end
        fprintf(fid,'-----x0=%d------ Gfunindex = %s ==================\r\n',x0(1), name_pena);
        fprintf(fid,'VRSPA       & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f  &%5.2f \r\n',   FunVRSPA/maxaver, GapVRSPA/maxaver, (Num_GradVRSPA/maxaver)/n,  Num_ProxVRSPA/maxaver,  Time_iterVRSPA/maxaver,  Time_stopVRSPA/maxaver);
        fprintf(fid,'SPDCA       & %5.4f   &%5.4e  &%5.4f   &%5.1f  &%5.2f  &%5.2f \r\n',   FunSPDCA/maxaver, GapSPDCA/maxaver, (Num_GradSPDCA/maxaver)/n,   Num_ProxSPDCA/maxaver, Time_iterSPDCA/maxaver,  Time_stopSPDCA/maxaver);
    end
end