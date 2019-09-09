% clear all
% close all
% clc
% ====================== Load data ==============================
load('hammer_cacti_bayer2.mat');
meas = meas_bayer(:, :, 4);
res = single(zeros(256, 256, 22, 4));
phi = mask;

for i=1:4
    if i==1
       bb                   =        double(meas(1:2:end, 1:2:end))       ;  
       mask                 =        double(phi(1:2:end, 1:2:end, :))  ;
    end  
    if i==2
       bb                   =        double(meas(2:2:end, 1:2:end))      ;  
       mask                 =        double(phi(2:2:end, 1:2:end, :)) ;
    end
    if i==3
       bb                   =        double(meas(1:2:end, 2:2:end))      ;  
       mask                 =        double(phi(1:2:end, 2:2:end, :)) ;
    end
    if i==4
       bb                   =        double(meas(2:2:end, 2:2:end))     ;  
       mask                 =        double(phi(2:2:end, 2:2:end, :)) ;
    end

    %Original               =     mask_bayer(:,:,1:11)                 ;
    [n1,n2,n3]             =     size(mask)                            ;     
    A                      =     diag(sparse(double(mask(1:n1*n2))))   ;  
    S = [];
    for j=2:n3
       S=diag(sparse(double(mask(n1*n2*(j-1)+1:n1*n2*j))))             ;
       A=[A,S];
    end    
   
      alpha                  =        1                                ;
      maxItr                 =        200                              ; % maximum iteration
      rho                    =        0.005                            ;
      b                      =        bb(:)                            ; % available data
% ================ main process of cmpletion =======================
      size(A)
      size(b)

      X                      =    tensor_cpl_admm( A , b , rho , alpha , ...
                                [n1,n2,n3] , maxItr );
      X                      =        abs(reshape(X,[n1,n2,n3]))       ;
      res(:, :, :, i)          =       X                                         ;
     %X_dif                  =        X-Original                       ;
     %RSE                    =        norm(X_dif(:))/norm(Original(:)) ;
    
end

save hammer_res res;

final_rec = hammer_demosaic(res);

save hammer_res res final_rec;
% ======================== Result Plot =============================
%for n=1:n3
%    psnr_temp(n) = psnr(double(X(:,:,n)),double(Original(:,:,n)),max(max(max(double(Original(:,:,n))))));
%end

% 
% temp = max(max(max(double(X))));
% 
% for n=1:n3
%     psnr_temp(n) = psnr(double(X(:,:,n))/temp,double(Original(:,:,n))/temp);
%     ssim_(n) = ssim(double(X(:,:,n))/temp,double(Original(:,:,n)/temp));
% end

%figure(1);
%for i = 1:n3  
%    subplot(4,6,i);imagesc(X(:,:,i));
%    axis off; 
%     subplot(2,8,i+8);imagesc(Original(:,:,i));
%     axis off; colormap(gray); 
%end



   




