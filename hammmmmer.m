% clear all 
% close all
% clc

%format short

% new algorithm for CACTI using TV denoising for each frame
% Xin Yuan; Bell Labs, Alcatel-Lucent
% xyuan@bell-labs.com
% initial date: 2015-7-02
%fname = '3Balls';
% fname = 'ColorData/3Balls';
% load(fname);
Row = 512;
Col = 512;

rotnum = 2;

y = Y(:,:,5:8);
ColT = 22;  % How many frames collapsed to 1 measurement
CodeFrame = 4;   % How many measurements to be used to reconstruct
T = ColT;
% 

y_R = zeros(Row/2, Col/2,CodeFrame);
y_B = zeros(Row/2, Col/2,CodeFrame);
y_G1 = zeros(Row/2, Col/2,CodeFrame);
y_G2 = zeros(Row/2, Col/2,CodeFrame);

Phi_R = zeros(Row/2, Col/2,ColT);
Phi_B = zeros(Row/2, Col/2,ColT);
Phi_G1 = zeros(Row/2, Col/2,ColT);
Phi_G2 = zeros(Row/2, Col/2,ColT);

% Xori_R = zeros(Row/2, Col/2,CodeFrame*ColT);
% Xori_B = zeros(Row/2, Col/2,CodeFrame*ColT);
% Xori_G1 = zeros(Row/2, Col/2,CodeFrame*ColT);
% Xori_G2 = zeros(Row/2, Col/2,CodeFrame*ColT);

for nR = 1:(Row/2)
    for nC = 1:(Col/2)
        y_R(nR,nC,:) = y((nR-1)*2+1,(nC-1)*2+1,:);
        y_B(nR,nC,:) = y((nR-1)*2+2,(nC-1)*2+2,:);
        y_G1(nR,nC,:) = y((nR-1)*2+1,(nC-1)*2+2,:);
        y_G2(nR,nC,:) = y((nR-1)*2+2,(nC-1)*2+1,:);
        Phi_R(nR,nC,:) = Phi((nR-1)*2+1,(nC-1)*2+1,:);
        Phi_B(nR,nC,:) = Phi((nR-1)*2+2,(nC-1)*2+2,:);
        Phi_G1(nR,nC,:) = Phi((nR-1)*2+1,(nC-1)*2+2,:);
        Phi_G2(nR,nC,:) = Phi((nR-1)*2+2,(nC-1)*2+1,:);
        
%         Xori_R(nR,nC,:) = X((nR-1)*2+1,(nC-1)*2+1,:);
%         Xori_B(nR,nC,:) = X((nR-1)*2+2,(nC-1)*2+2,:);
%         Xori_G1(nR,nC,:) = X((nR-1)*2+1,(nC-1)*2+2,:);
%         Xori_G2(nR,nC,:) = X((nR-1)*2+2,(nC-1)*2+1,:);
    end
end


UseF = CodeFrame*T;

%% set parameters
row = Row/2;
col = Col/2;
para.row = row;
para.col = col;
para.iter = 200;
para.lambda = 0.2;  
para.TVweight = 0.07;img

 % clear y
  %para.CSr = CSr;

    for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
        
        for rr=1:4
        %Recon_X = zeros(Row, Col, T);
            switch rr
                case 1
                      yuse = y_R(:,:,k);
                      Phi_use = Phi_R;
                      %para.ori_im = Xori_R(:,:,(k-1)*ColT+(1:ColT));
                case 2
                      yuse = y_G1(:,:,k);
                      Phi_use = Phi_G1;
                      %para.ori_im = Xori_G1(:,:,(k-1)*ColT+(1:ColT));
                case 3
                      yuse = y_G2(:,:,k);
                      Phi_use = Phi_G2;
                      %para.ori_im = Xori_G2(:,:,(k-1)*ColT+(1:ColT));
                case 4
                      yuse = y_B(:,:,k);
                      Phi_use = Phi_B;
                      %para.ori_im = Xori_B(:,:,(k-1)*ColT+(1:ColT));
            end
            
            A = @(z) A_xy(z, Phi_use);
        At = @(z) At_xy_nonorm(z, Phi_use);

         theta    =   TV_IST_CACTI( yuse,1, para, A,At);
         if(mod(k,2)==0)
         img_ist(:,:,(k-1)*ColT+(1:ColT),rr) = theta;
         else
            img_ist(:,:,(k-1)*ColT+(ColT:-1:1),rr) = theta;
         end
        end
    end


    
 [Img_recon_sensor_ist, X_recon_ist] = my_demosaic(Row,Col,CodeFrame,ColT,rotnum,img_ist);
 savename = [fname 'IST-TV_R' num2str(Row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
 CACTI_gene_moive(savename,CodeFrame,ColT,y,Img_recon_sensor_ist,X_recon_ist,rotnum);
 %[Img_recon_sensor_ist, X_recon_ist] = my_demosaic_reverse(Row,Col,CodeFrame,ColT,rotnum,img_ist);
%[PSNR, PSNR_mos] = my_genneratemoive(img_ist, Row, Col, CodeFrame, ColT,X,  y,X_ori_color,X_mos,fname,'IST-TV');
 
 %%
 
 para.iter = 200;
 para.acc= 0;
para.lambda = 1;
for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
         for rr=1:4
        %Recon_X = zeros(Row, Col, T);
            switch rr
                case 1
                      yuse = y_R(:,:,k);
                      Phi_use = Phi_R;
                      %para.ori_im = Xori_R(:,:,(k-1)*ColT+(1:ColT));
                case 2
                      yuse = y_G1(:,:,k);
                      Phi_use = Phi_G1;
                      %para.ori_im = Xori_G1(:,:,(k-1)*ColT+(1:ColT));
                case 3
                      yuse = y_G2(:,:,k);
                      Phi_use = Phi_G2;
                      %para.ori_im = Xori_G2(:,:,(k-1)*ColT+(1:ColT));
                case 4
                      yuse = y_B(:,:,k);
                      Phi_use = Phi_B;
                      %para.ori_im = Xori_B(:,:,(k-1)*ColT+(1:ColT));
            end
            Phi_sum = sum(Phi_use.^2,3);
            Phi_sum(Phi_sum==0)=1;
  
             para.Phi_sum = Phi_sum;
            A = @(z) A_xy(z, Phi_use);
          At = @(z) At_xy_nonorm(z, Phi_use);

         theta    =   TV_GAP_CACTI( yuse,1, para, A,At);
         if(mod(k,2)==0)
         img_gap(:,:,(k-1)*ColT+(1:ColT),rr) = theta;
         else
        img_gap(:,:,(k-1)*ColT+(ColT:-1:1),rr) = thetisa;     
         end
         end
  
 end

  

 [Img_recon_sensor_gap, X_recon_gap] = my_demosaic(Row,Col,CodeFrame,ColT,rotnum,img_gap);
  savename = [fname 'GAP-TV_R' num2str(Row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
 CACTI_gene_moive(savename,CodeFrame,ColT,y,Img_recon_sensor_gap,X_recon_gap,rotnum);
 %%
  para.acc= 1;
  para.iter =250;
for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
        for rr=1:4
         switch rr
                case 1
                      yuse = y_R(:,:,k);
                      Phi_use = Phi_R;
                    %  para.ori_im = Xori_R(:,:,(k-1)*ColT+(1:ColT));
                case 2
                      yuse = y_G1(:,:,k);
                      Phi_use = Phi_G1;
                      %para.ori_im = Xori_G1(:,:,(k-1)*ColT+(1:ColT));
                case 3
                      yuse = y_G2(:,:,k);
                      Phi_use = Phi_G2;
                     % para.ori_im = Xori_G2(:,:,(k-1)*ColT+(1:ColT));
                case 4
                      yuse = y_B(:,:,k);
                      Phi_use = Phi_B;
                     % para.ori_im = Xori_B(:,:,(k-1)*ColT+(1:ColT));
          end
            Phi_sum = sum(Phi_use.^2,3);
            Phi_sum(Phi_sum==0)=1;
  
             para.Phi_sum = Phi_sum;
            A = @(z) A_xy(z, Phi_use);
         At = @(z) At_xy_nonorm(z, Phi_use);

         theta    =   TV_GAP_CACTI( yuse,1, para, A,At);
         if(mod(k,2)==0)
          img_gap_acc(:,:,(k-1)*ColT+(1:ColT),rr) = theta;
         else
          img_gap_acc(:,:,(k-1)*ColT+(ColT:-1:1),rr) = theta;     
         end
        
  
        end
 end

  [Img_recon_sensor_gap_acc, X_recon_gap_acc] = my_demosaic(Row,Col,CodeFrame,ColT,rotnum,img_gap_acc);
   savename = [fname 'Acc-GAP-TV_R' num2str(Row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
 CACTI_gene_moive(savename,CodeFrame,ColT,y,Img_recon_sensor_gap_acc,X_recon_gap_acc,rotnum);
 % ADMM
% Phi_sum = sum(Phi.^2,3);
% %Phi_sum(Phi_sum==0)=1;
% % para.lambda = 1;
%  para.Phi_sum = Phi_sum;
 para.eta = 1e-3;
 para.iter = 200;
% para.acc= 0;
for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
        for rr=1:4
        switch rr
                case 1
                      yuse = y_R(:,:,k);
                      Phi_use = Phi_R;
                      %para.ori_im = Xori_R(:,:,(k-1)*ColT+(1:ColT));
                case 2
                      yuse = y_G1(:,:,k);
                      Phi_use = Phi_G1;
                     % para.ori_im = Xori_G1(:,:,(k-1)*ColT+(1:ColT));
                case 3
                      yuse = y_G2(:,:,k);
                      Phi_use = Phi_G2;
                     % para.ori_im = Xori_G2(:,:,(k-1)*ColT+(1:ColT));
                case 4
                      yuse = y_B(:,:,k);
                      Phi_use = Phi_B;
                     % para.ori_im = Xori_B(:,:,(k-1)*ColT+(1:ColT));
         end
            Phi_sum = sum(Phi_use.^2,3);
            %Phi_sum(Phi_sum==0)=1;
  
             para.Phi_sum = Phi_sum;
            A = @(z) A_xy(z, Phi_use);
         At = @(z) At_xy_nonorm(z, Phi_use);

        %theta    =   DCT_CACTI_weight( y_use, para, A,At);
        theta    =   TV_ADMM_CACTI( yuse,1, para, A,At);
        if(mod(k,2)==0)
        img_admm(:,:,(k-1)*ColT+(1:ColT),rr) = theta;
        else
        img_admm(:,:,(k-1)*ColT+(ColT:-1:1),rr) = theta;   
        end
        end
 end

    
   [Img_recon_sensor_admm, X_recon_admm] = my_demosaic(Row,Col,CodeFrame,ColT,rotnum,img_admm);
    savename = [fname 'ADMM-TV_R' num2str(Row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
 CACTI_gene_moive(savename,CodeFrame,ColT,y,Img_recon_sensor_admm,X_recon_admm,rotnum);
 % now we try gap_wavelet
 addpath(genpath('../CACTI_gap/'));
 stopc.iternum = 100;
stopc.err = 10^-5;
acc = 2;
ydim = Row/2*Col/2;


spbasis.space = 'wavelet';  % transform for space, 'wavelet' or 'dct'
spbasis.time  = 'dct';  % transform for spectrum, 'wavelet' or 'dct', dct is always used. If we use wavelet, T need to be the power of 2. we can use no, means no transformation
%spbasis.spectrum  = 'no';  % Here we use no, means no transfromation in spectrum

weighttype.space = 'tree';   % Here we can select:  'tree' or 'block'
weighttype.time = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2


weight_base.type = 'cos'; % here we can select 'exp' or 'cos'
if strcmp(weight_base.type,'exp')
    weight_base.space = 2;   % This weight is the base of exponential decay. should be large than 1 [1 2] is always used
    weight_base.time = 2;
end



% The block size for group
block.row = 2;
block.col = 2;
block.T = 6;

m_star = ceil(ydim/(block.row*block.col*block.T));

    for k=1:CodeFrame        
        fprintf('----- Reconstructing frame-block %d of %d\n',k,CodeFrame);
       for rr=1:4
        switch rr
                case 1
                      yuse = y_R(:,:,k);
                      Phi_use = Phi_R;
                     % para.ori_im = Xori_R(:,:,(k-1)*ColT+(1:ColT));
                case 2
                      yuse = y_G1(:,:,k);
                      Phi_use = Phi_G1;
                     % para.ori_im = Xori_G1(:,:,(k-1)*ColT+(1:ColT));
                case 3
                      yuse = y_G2(:,:,k);
                      Phi_use = Phi_G2;
                      %para.ori_im = Xori_G2(:,:,(k-1)*ColT+(1:ColT));
                case 4
                      yuse = y_B(:,:,k);
                      Phi_use = Phi_B;
                     % para.ori_im = Xori_B(:,:,(k-1)*ColT+(1:ColT));
        end
  
        
        theta_wL21 = GAP_3D_wL21_grayscale(yuse,Phi_use, Row/2,Col/2,ColT, block,spbasis, m_star,stopc,acc,weight_base,weighttype); 
        if(mod(k,2)==0)
        img_gap_w(:,:,(k-1)*ColT+(1:ColT),rr) = theta_wL21;
        else
            img_gap_w(:,:,(k-1)*ColT+(ColT:-1:1),rr) = theta_wL21;
        end
       end
      
    end
      [Img_recon_sensor_gap_w, X_recon_gap_w] = my_demosaic(Row,Col,CodeFrame,ColT,rotnum,img_gap_w);
      savename = [fname 'GAP-w_R' num2str(Row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
 CACTI_gene_moive(savename,CodeFrame,ColT,y,Img_recon_sensor_gap_w,X_recon_gap_w,rotnum);

    savename = ['3balls_row' num2str(row) '_T' num2str(ColT) '_F' num2str(CodeFrame)];
%    saveas(gcf,[ savename '.fig']);
% saveas(gcf,[ savename '.png']);

save([savename '.mat'],'X_recon_ist','X_recon_gap','X_recon_admm','X_recon_gap_acc','X_recon_gap_w');

