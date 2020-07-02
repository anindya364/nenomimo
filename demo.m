clc
clear

%read the MRI image
B=imread('demo.jpg');
B=rgb2gray(B);
[row,col]=size(B);
imageSize=size(B);
nsize=row*col;
Y=reshape(B,nsize,1);
Y=double(Y);

%set the cluster number
K=3;

Nsamp=500;
z=randi(K,1,nsize);
zz=reshape(z,nsize,1);

%run the gibbsMSNBURRMM function
[miu,phi,al,p,z,churn] = gibbsMSNBURRMM(Y,K,10,Nsamp,z);

%%

ima=double(B);
%ima=imresize(ima,[512 512]);
copy=ima;
copy2=ima;% make a copy
ima=ima(:);         % vectorize ima
mi=min(ima);        % deal with negative 
ima=ima-mi+1;       % and zero values
m=max(ima);
s=length(ima);

%%
%3 cluster
mask = reshape(z,row,col);
for i=1:row
    for j=1:col
        if mask(i,j)==1
            final_img(i,j)=0;
        elseif mask(i,j)==2
            final_img(i,j)=255;
        else
            final_img(i,j)=128;
        end
    end
end
figure,imshow(final_img/255,[]); 

%% 4 CLUSTER
% mask = reshape(z,row,col);
% for i=1:row
%   for j=1:col
%       if mask(i,j)==1
%           final_img(i,j)=255;
%       elseif mask(i,j)==2
%           final_img(i,j)=200;
%       elseif mask(i,j)==3
%           final_img(i,j)=0;
%       else
%           final_img(i,j)=120;
%       end
%   end
% end
% figure,imshow(final_img/255,[]); 
%%
% 5 CLUSTER
% mask = reshape(z,row,col);
% for i=1:row
%   for j=1:col
%       if mask(i,j)==1
%           final_img(i,j)=255;
%       elseif mask(i,j)==2
%           final_img(i,j)=100;
%       elseif mask(i,j)==3
%           final_img(i,j)=0;
%       elseif mask(i,j)==4
%           final_img(i,j)=200;
%       else
%           final_img(i,j)=160;
%       end
%   end
% end
% figure,imshow(final_img/255,[]); 

%%
%subjective analysis
orgimg=copy;
im_bin=mask;
[r,c]=find(im_bin==3);
rc=[r,c];
for j=1:(numel(rc)/2)
    copy(r(j),c(j))=0;
end

figure, subplot(2,1,1);
imshow(copy,[]), title('NROI-Non region of interest');
nroi_image=copy; 
logimg=imsubtract(orgimg,copy);
roi_image=logimg;
subplot(2,1,2);
imshow(logimg,[]), title('ROI-region of interest');
 
%%
%OSC calculation for 3 cluster
data=Y;
klas=zz;
SI=silhouette(data,klas);
SI_cluster=[mean(SI(klas==1)) mean(SI(klas==2)) mean(SI(klas==3))];
SI_all=mean(SI_cluster);


