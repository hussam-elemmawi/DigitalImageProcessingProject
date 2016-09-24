 I = imread('cameraman.tif');
 y=imfinfo('cameraman.tif');
 display(y);
 I = im2double(I);%el dctmtx betala3 double 
 T = dctmtx(8);% e3mel dct of 8*8
 B = blkproc(I,[8 8],'P1*x*P2',T,T');%transform to ferq domain
 mask = [1 1 1 1 0 0 0 0; 1 1 1 0 0 0 0 0; 1 1 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];
 B2 = blkproc(B,[8 8],'P1.*x',mask);% ferq domain image multibly by mask
 I2 = blkproc(B2,[8 8],'P1*x*P2',T',T);%bakc to time domain
 imshow(I), figure, imshow(I2)
 %x=imfinfo('cameraman.tif');
% display(x);


% 
% function y = DCTImpl(x)
% N = length(x);
% if N == 1
%     y = x;
% else
%     x1 = [x(1:2:(N-1));
%         x(N:(-2):2)];
%     y = FFTImpl(x1);
%     rp = real(y);
%     ip = imag(y);
%     y = cos(pi*((0:(N-1))’)/(2*N)).*rp + sin(pi*((0:(N-1))’)/(2*N)).*ip;
%     y(2:N) = sqrt(2)*y(2:N);
% end

% 
%   function image_comp = dctII(image, b)
%     [h w] = size(image);
%     image = double(image) - 128;
%     block = zeros(b,b);
% 
%  image_t=zeros(size(image));
%  for k=1:b:h
%      for l=1:b:w
%         image_t(k:k+b-1,l:l+b-1)= image(k:k+b-1,l:l+b-1);
%         for u=1:b
%             for v=1:b
%                 if u == 0
%                     Cu = 1/sqrt(2);
%                 else
%                     Cu = 1;
%                 end
%                 if v == 0
%                     Cv = 1/sqrt(2);
%                 else
%                     Cv = 1;
%                 end
%                 Res_sum=0;
%                 for x=1:b;
%                     for y=1:b
%                         Res_sum = Res_sum + ((image_t(x,y))*cos(((2*x)+1)*u*pi/(2*b))*cos(((2*y)+1)*v*pi/(2*b)));  
%                     end
%                 end
%                 dct= (1/4)*Cu*Cv*Res_sum;
%                 block(u,v) = dct;
% 
%             end
%         end
%         image_comp(k:k+b-1,l:l+b-1)=block(u,v);
%      end
%  end
% end

%I = imread('cameraman.tif');

%I = im2double(I);
%T = dctmtx(8);
%dct = @(block_struct) T * block_struct.data * T';
%B = blockproc(I,[8 8],dct);
%mask = [1   1   1   1   0   0   0   0
 %       1   1   1   0   0   0   0   0
 %       1   1   0   0   0   0   0   0
  %      1   0   0   0   0   0   0   0
   %     0   0   0   0   0   0   0   0
    %    0   0   0   0   0   0   0   0
     %   0   0   0   0   0   0   0   0
      %  0   0   0   0   0   0   0   0];
%B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);
%invdct = @(block_struct) T' * block_struct.data * T;
%I2 = blockproc(B2,[8 8],invdct);
%imshow(I), figure, imshow(I2)