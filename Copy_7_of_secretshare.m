%  ������ʵ��ͼ��LENA��ѹ�����У�������С����ִ�У�
%  �������ߣ�ɳ������۴�ѧ�������ӹ���ѧϵ��wsha@eee.hku.hk
%  �㷨��������ƥ�䷨���ο����� Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit��IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
%  �ó���û�о����κ��Ż�

function Copy_7_of_secretshare

clc;clear

%%  CS����

%  ���ļ�
Z=imread('lena256.bmp');

    X=Z;
%     X= rgb2gray(X);
    X=double(X);

[a,b]=size(X);



%  ������������(��֪)
M=240;
k=0.3;    %chebyshv���β���
distance=4; %�������
n=M*a*distance;


x(1)=0.3121212;
for i=1:n-1
   if(x(i)>0 & x(i)<k)
        x(i+1)=x(i)/k;
    else
        x(i+1)=(1-x(i))/(1-k);
    end   %tent��ӳ���ϵ
%     xx(i)=x(i)-0.5;
end
x_sample=downsample(x,distance);%�����������������м��ɵõ�����Ļ������
p=1/a;
RR(:,:)=p*reshape(x_sample,M,a);   %ѹ��ʱʹ�õľ���
RR1=65000*RR;

R=reshape(x_sample,M,a);   %ѹ��ʱʹ�õľ���


hh4=99;
hh5=179;
% for k=1:hh5-hh4+1   
% h4=hh4+k-1
% h5=h4+1

h1=104;
h2=192;

%�ڸ���Ϣ
  

 R1=R;
 R2=R;
 
 t(1)=0.57+0.3*rand;
     tt(1)=1;
    for j=1:256      
     tt(j)=0.57+0.3*rand;
    end
 
 
for i=h1:h2  %  ��ѭ��  

    R1(:,i)=-tt(i)*R(:,i);
%     R2(:,i)=-R(:,i);%������ʹ��R2�ָ�,��ʽ1
    R2(:,i)=-rand*R(:,i);%������ʹ��R2�ָ�����ʽ2
%               R(:,i)=-1*R(:,i);
    end
 

 
    



%  С���任��������
ww=DWT(a);

%  �ź�С���任
%%X1=ww*X*ww';




X1=ww*X;




%  ����(��֪)����С�������
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
%%  CS�ָ�

%  OMP�㷨
% X2=zeros(a,b);  %  �ָ�����
% for i=hh4:hh5   %  ��ѭ��       
%     rec=omp(Y3(:,i),R1,a);
%     X2(:,i)=rec;
% end


for i=hh4:hh5  %  ��ѭ�� 
    rec=omp(Y3(:,i),R,a);
%     rec=omp(Y3(:,i),R2,a);
    X4(:,i)=rec;
end

for i=hh4:hh5   %  ��ѭ��       
    rec=omp(Y3(:,i),R1,a);
    X2(:,i)=rec;
end

    for i=1:(hh4-1)  %  ��ѭ��       
    rec=omp(Y3(:,i),R,a);
    X2(:,i)=rec;

end

for i=1:(hh4-1)  %  ��ѭ��       
    rec=omp(Y3(:,i),R,a);
    X4(:,i)=rec;

end

for i=(hh5+1):256  %  ��ѭ��       
    rec=omp(Y3(:,i),R,a);
    X2(:,i)=rec;

end

for i=(hh5+1):256  %  ��ѭ��       
    rec=omp(Y3(:,i),R,a);
    X4(:,i)=rec;

end


% X4=zeros(a,b);  %  �ָ�����
% for i=1:256  %  ��ѭ��       
%     rec=omp(Y3(:,i),R,a);
%     X4(:,i)=rec;
% end

%%  CS�����ʾ

%  ԭʼͼ��
figure(1);
imshow(uint8(Z));
% title('ԭʼͼ��');

% %  �任ͼ��
% figure(2);
% imshow(uint8(X2));
% title('�ָ�С�����ͼ��');

%  ѹ�����лָ���ͼ��

%%X3=ww'*sparse(X2)*ww;  %  С�����任
X3=ww'*sparse(X2)
X3=full(X3);
W3(:,:)=uint8(X3);



%%X3=ww'*sparse(X2)*ww;  %  С�����任
X5=ww'*sparse(X4)
X5=full(X5);
W5(:,:)=uint8(X5);


% end

Y4=0.0001*Y5+65000*RR; %û������ʹ�ã�ֻΪ��ͼ

 figure(7)
 imshow(Y4)


figure(3);
imshow(uint8(W3));
% title('�ָ��Ŀ���ͼ��-full');

figure(4);
imshow(uint8(W5));
% title('�ָ��Ŀ���ͼ��-semi');

% figure(5)
% title('�����ͼ��ֱ��ͼ-gray');
% Y4gray=rgb2gray(Y4);
% imshow(uint8(Y4gray));
% % imhist(uint8(Y4gray),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])


figure(12)
% title('�����ͼ��ֱ��ͼ-r');
imhist(uint8(Y4(:,:)),a);
set(gca,'FontName','Times New Roman','FontSize',10)
set(gcf,'unit','centimeters','position',[10 5 8 6])

% figure(13)
% title('�����ͼ��ֱ��ͼ-g');
% imhist(uint8(Y4(:,:,2)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% 
% 
% figure(14)
% title('�����ͼ��ֱ��ͼ-b');
% imhist(uint8(Y4(:,:,3)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])



% figure(6)
% title('ԭʼͼ��ֱ��ͼ-gray');
% Zgray=rgb2gray(Z);
% imhist(uint8(Zgray),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');

figure(9)
% title('ԭʼͼ��ֱ��ͼ-r');
imhist(uint8(Z(:,:)),a);
set(gca,'FontName','Times New Roman','FontSize',10)
set(gcf,'unit','centimeters','position',[10 5 8 6])
ylim([1,1450])
%title('histogram of original iamge');

% figure(10)
% title('ԭʼͼ��ֱ��ͼ-g');
% imhist(uint8(Z(:,:,2)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');
% 
% figure(11)
% title('ԭʼͼ��ֱ��ͼ-b');
% imhist(uint8(Z(:,:,3)),a);
% set(gca,'FontName','Times New Roman','FontSize',10)
% set(gcf,'unit','centimeters','position',[10 5 8 6])
% ylim([1,1450])
% %title('histogram of original iamge');
% 




figure(8);
W_5_semi=imcrop(W5,[hh4,h1,hh5-hh4,h2-h1]);
imshow(uint8(W_5_semi));
title('�ָ��Ŀ���ͼ��-��semi');


figure(15);
Z_semi=imcrop(Z,[hh4,h1,hh5-hh4,h2-h1]);
imshow(uint8(Z_semi));
title('ԭʼ�Ŀ���ͼ��-��semi');


% ����Ϣ��
I=uint8(Y4);
[C,L]=size(I); %��ͼ��Ĺ��
Img_size=C*L; %ͼ�����ص���ܸ���
G=256; %ͼ��ĻҶȼ�
H_x=0;
nk=zeros(G,1);%����һ��G��1�е�ȫ�����
for i=1:C
for j=1:L
Img_level=I(i,j)+1; %��ȡͼ��ĻҶȼ�
nk(Img_level)=nk(Img_level)+1; %ͳ��ÿ���Ҷȼ����صĵ���
end
end
for k=1:G  %ѭ��
Ps(k)=nk(k)/Img_size; %����ÿһ�����ص�ĸ���
if Ps(k)~=0; %������ص�ĸ��ʲ�Ϊ��
H_x=-Ps(k)*log2(Ps(k))+H_x; %����ֵ�Ĺ�ʽ
end
end
H_x_full=H_x
H_x_full  %��ʾ��ֵ

%  ���(PSNR)
Z=double(Z);
W3=double(W3);
errorx=sum(sum(abs(W3-Z).^2));        %  MSE���
psnr=10*log10(255*255/(errorx/a/b))   %  PSNR
% psnr_full=(psnr(:,:,1)+psnr(:,:,2)+psnr(:,:,3))/3
% fprintf('psnr(:,:,1)=%f',psnr(:,:,1))
% fprintf('psnr(:,:,2)=%f',psnr(:,:,2))
% fprintf('psnr(:,:,3)=%f',psnr(:,:,3))

%  OMP�ĺ���
%  s-������T-�۲����N-������С
function hat_y=omp(s,T,N)

Size=size(T);                                     %  �۲�����С
M=Size(1);                                        %  ����
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=s;                                            %  �в�ֵ

for times=1:M;                                    %  �������������ᳬ������ֵ��

    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
    
    if (abs(aug_y(end))^2/norm(aug_y)<0.05)       %  ����Ӧ�ض���***��Ҫ��������ֵ��
        break;
    end
end

hat_y(pos_array)=aug_y;                           %  �ع�������


