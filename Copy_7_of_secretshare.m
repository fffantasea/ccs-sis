%  本程序实现图像LENA的压缩传感（测量在小波域执行）
%  程序作者：沙威，香港大学电气电子工程学系，wsha@eee.hku.hk
%  算法采用正交匹配法，参考文献 Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit，IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
%  该程序没有经过任何优化

function Copy_7_of_secretshare

clc;clear

%%  CS测量

%  读文件
Z=imread('lena256.bmp');

    X=Z;
%     X= rgb2gray(X);
    X=double(X);

[a,b]=size(X);



%  测量矩阵生成(已知)
M=240;
k=0.3;    %chebyshv分形参数
distance=4; %采样间隔
n=M*a*distance;


x(1)=0.3121212;
for i=1:n-1
   if(x(i)>0 & x(i)<k)
        x(i+1)=x(i)/k;
    else
        x(i+1)=(1-x(i))/(1-k);
    end   %tent的映射关系
%     xx(i)=x(i)-0.5;
end
x_sample=downsample(x,distance);%该向量经过重新排列即可得到所需的混沌矩阵
p=1/a;
RR(:,:)=p*reshape(x_sample,M,a);   %压缩时使用的矩阵
RR1=65000*RR;

R=reshape(x_sample,M,a);   %压缩时使用的矩阵


hh4=99;
hh5=179;
% for k=1:hh5-hh4+1   
% h4=hh4+k-1
% h5=h4+1

h1=104;
h2=192;

%掩盖信息
  

 R1=R;
 R2=R;
 
 t(1)=0.57+0.3*rand;
     tt(1)=1;
    for j=1:256      
     tt(j)=0.57+0.3*rand;
    end
 
 
for i=h1:h2  %  列循环  

    R1(:,i)=-tt(i)*R(:,i);
%     R2(:,i)=-R(:,i);%攻击者使用R2恢复,方式1
    R2(:,i)=-rand*R(:,i);%攻击者使用R2恢复，方式2
%               R(:,i)=-1*R(:,i);
    end
 

 
    



%  小波变换矩阵生成
ww=DWT(a);

%  信号小波变换
%%X1=ww*X*ww';




X1=ww*X;




%  测量(已知)，在小波域测量
%%Y=R*X1;
R=R*ww'
R1=R1*ww'
Y(:,:)=R*X1
Y1(:,:)=R1*X1






Y3(:,hh4:hh5)=Y1(:,hh4:hh5);
Y3(:,1:(hh4-1))=Y(:,1:(hh4-1));
Y3(:,(hh5+1):a)=Y(:,(hh5+1):a);
Y5=Y3;
% Y3=Y3(1:204,:);
% R=R(1:204,:);
% R1=R1(1:204,:);




% Y4=Y3;
%%  CS恢复

%  OMP算法
% X2=zeros(a,b);  %  恢复矩阵
% for i=hh4:hh5   %  列循环       
%     rec=omp(Y3(:,i),R1,a);
%     X2(:,i)=rec;
% end


for i=hh4:hh5  %  列循环 
    rec=omp(Y3(:,i),R,a);
%     rec=omp(Y3(:,i),R2,a);
    X4(:,i)=rec;
end

for i=hh4:hh5   %  列循环       
    rec=omp(Y3(:,i),R1,a);
    X2(:,i)=rec;
end

    for i=1:(hh4-1)  %  列循环       
    rec=omp(Y3(:,i),R,a);
    X2(:,i)=rec;

end

for i=1:(hh4-1)  %  列循环       
    rec=omp(Y3(:,i),R,a);
    X4(:,i)=rec;

end

for i=(hh5+1):256  %  列循环       
    rec=omp(Y3(:,i),R,a);
    X2(:,i)=rec;

end

for i=(hh5+1):256  %  列循环       
    rec=omp(Y3(:,i),R,a);
    X4(:,i)=rec;

end


% X4=zeros(a,b);  %  恢复矩阵
% for i=1:256  %  列循环       
%     rec=omp(Y3(:,i),R,a);
%     X4(:,i)=rec;
% end

%%  CS结果显示

%  原始图像
figure(1);
imshow(uint8(Z));
% title('原始图像');

% %  变换图像
% figure(2);
% imshow(uint8(X2));
% title('恢复小波域的图像');

%  压缩传感恢复的图像

%%X3=ww'*sparse(X2)*ww;  %  小波反变换
X3=ww'*sparse(X2)
X3=full(X3);
W3(:,:)=uint8(X3);



%%X3=ww'*sparse(X2)*ww;  %  小波反变换
X5=ww'*sparse(X4)
X5=full(X5);
W5(:,:)=uint8(X5);


% end

Y4=0.0001*Y5+65000*RR; %没有其他使用，只为作图

 figure(7)
 imshow(Y4)


figure(3);
imshow(uint8(W3));
% title('恢复的空域图像-full');

figure(4);
imshow(uint8(W5));
% title('恢复的空域图像-semi');

% figure(5)
% title('混沌后图像直方图-gray');
% Y4gray=rgb2gray(Y4);
% imshow(uint8(Y4gray));
% % imhist(uint8(Y4gray),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])


figure(12)
% title('混沌后图像直方图-r');
imhist(uint8(Y4(:,:)),a);
set(gca,'FontName','Times New Roman','FontSize',10)
set(gcf,'unit','centimeters','position',[10 5 8 6])

% figure(13)
% title('混沌后图像直方图-g');
% imhist(uint8(Y4(:,:,2)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% 
% 
% figure(14)
% title('混沌后图像直方图-b');
% imhist(uint8(Y4(:,:,3)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])



% figure(6)
% title('原始图像直方图-gray');
% Zgray=rgb2gray(Z);
% imhist(uint8(Zgray),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');

figure(9)
% title('原始图像直方图-r');
imhist(uint8(Z(:,:)),a);
set(gca,'FontName','Times New Roman','FontSize',10)
set(gcf,'unit','centimeters','position',[10 5 8 6])
ylim([1,1450])
%title('histogram of original iamge');

% figure(10)
% title('原始图像直方图-g');
% imhist(uint8(Z(:,:,2)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');
% 
% figure(11)
% title('原始图像直方图-b');
% imhist(uint8(Z(:,:,3)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');
% 




figure(8);
W_5_semi=imcrop(W5,[hh4,h1,hh5-hh4,h2-h1]);
imshow(uint8(W_5_semi));
title('恢复的空域图像-仅semi');


figure(15);
Z_semi=imcrop(Z,[hh4,h1,hh5-hh4,h2-h1]);
imshow(uint8(Z_semi));
title('原始的空域图像-仅semi');


% 求信息熵
I=uint8(Y4);
[C,L]=size(I); %求图像的规格
Img_size=C*L; %图像像素点的总个数
G=256; %图像的灰度级
H_x=0;
nk=zeros(G,1);%产生一个G行1列的全零矩阵
for i=1:C
for j=1:L
Img_level=I(i,j)+1; %获取图像的灰度级
nk(Img_level)=nk(Img_level)+1; %统计每个灰度级像素的点数
end
end
for k=1:G  %循环
Ps(k)=nk(k)/Img_size; %计算每一个像素点的概率
if Ps(k)~=0; %如果像素点的概率不为零
H_x=-Ps(k)*log2(Ps(k))+H_x; %求熵值的公式
end
end
H_x_full=H_x
H_x_full  %显示熵值

%  误差(PSNR)
Z=double(Z);
W3=double(W3);
errorx=sum(sum(abs(W3-Z).^2));        %  MSE误差
psnr=10*log10(255*255/(errorx/a/b))   %  PSNR
% psnr_full=(psnr(:,:,1)+psnr(:,:,2)+psnr(:,:,3))/3
% fprintf('psnr(:,:,1)=%f',psnr(:,:,1))
% fprintf('psnr(:,:,2)=%f',psnr(:,:,2))
% fprintf('psnr(:,:,3)=%f',psnr(:,:,3))

%  OMP的函数
%  s-测量；T-观测矩阵；N-向量大小
function hat_y=omp(s,T,N)

Size=size(T);                                     %  观测矩阵大小
M=Size(1);                                        %  测量
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=s;                                            %  残差值

for times=1:M;                                    %  迭代次数（不会超过测量值）

    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    
    if (abs(aug_y(end))^2/norm(aug_y)<0.05)       %  自适应截断误差（***需要调整经验值）
        break;
    end
end

hat_y(pos_array)=aug_y;                           %  重构的向量


