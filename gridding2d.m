%generate image
clear; clc;
addpath("D:\SemiPark\gridding\nufft_toolbox")
imsize = 256;
image = phantom(imsize);
image = fftshift(fft(fftshift(image,1),[],1),1);
image = fftshift(fft(fftshift(image,2),[],2),2);
image = ifftshift(ifft(ifftshift(image,1),[],1),1);
image = ifftshift(ifft(ifftshift(image,2),[],2),2);
nviews = 400;
nsamps =512;
%kspace tajectory
ang = linspace(0,2*pi,nviews)';
rad = linspace(-0.5,0.5,nsamps)';
ktraj = zeros(nsamps*nviews,2);
for i = 1:nviews
    m = repmat([cos(ang(i)) sin(ang(i))],[nsamps 1]);
    ktraj((i-1)*nsamps+1:i*nsamps,:) = m;
end
ktraj = repmat(rad,[nviews 2]).*ktraj;
om = ktraj*2*pi;
Nd = [imsize, imsize];
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
st = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
%image to kspace
kdata = nufft(image,st);
kdata_figure = reshape(kdata, nsamps,nviews);
%density compensation
cd = repmat(abs(linspace(-1,1,nsamps))',nviews,1);
%%
kdata = kdata.*cd;
ktraj_grid = ktraj*(nsamps-3)+(nsamps+1)/2;

ktraj_x = ktraj_grid(:,1);
ktraj_y = ktraj_grid(:,2);
clear ktraj_grid
fx = floor(ktraj_x);  cx = ceil(ktraj_x); 
x_index = repmat([fx-1 fx cx cx+1],[1 4]);
fy = floor(ktraj_y); cy = ceil(ktraj_y);
y_index = horzcat(repmat(fy-1,[1 4]),repmat(fy,[1 4]),repmat(cy,[1 4]),repmat(cy+1,[1 4]));
clear cx cy fx fy
% ktraj_repx = repmat(ktraj_x,[1 16]);
% ktraj_repy = repmat(ktraj_y,[1 16]);
dist_x = abs(x_index-ktraj_x);
dist_y = abs(y_index-ktraj_y);
distsq_x = dist_x.^2;
distsq_y = dist_y.^2;
distsq = distsq_x + distsq_y;
dist = sqrt(distsq);
clear distsq_x distsq_y distsq ktraj_x ktraj_y

kdata = repmat(kdata,[1 16]);
index_smth2 = find(dist <= 2);

x_index = x_index(index_smth2); y_index = y_index(index_smth2);
dist_x = dist_x(index_smth2); dist_y = dist_y(index_smth2);
dist = dist(index_smth2);
kdata = kdata(index_smth2);

%%
beta = 5.7567;
window = kaiser(40003,beta);
lookuptable = window(20002:40003);
dist_table_x = dist_x*10000+1; dist_table_y = dist_y*10000+1;
kerneldistance_x = (dist_table_x-floor(dist_table_x)).*lookuptable(ceil(dist_table_x))+(ceil(dist_table_x)-dist_table_x).*lookuptable(floor(dist_table_x));
kerneldistance_y = (dist_table_y-floor(dist_table_y)).*lookuptable(ceil(dist_table_y))+(ceil(dist_table_y)-dist_table_y).*lookuptable(floor(dist_table_y));
kerneldistance = kerneldistance_x.*kerneldistance_y;

final_value = kerneldistance.*kdata;
%%
kspace_data = zeros(nsamps,nsamps);
% kspace_data_padding = zeros(nsamps*2,nsamps*2);
% mask = zeros(nsamps,nsamps);
for i  = 1:size(index_smth2,1)
    kspace_data(x_index(i),y_index(i)) = kspace_data(x_index(i),y_index(i))+final_value(i);
end
kspace_data2 = accumarray(cat(2,x_index,y_index),final_value,[512 512],@sum);

IMG = kspace_data;
IMG = ifftshift(ifft(ifftshift(IMG,1),[],1),1);
IMG = ifftshift(ifft(ifftshift(IMG,2),[],2),2);
% IMG = IMG(129:384,129:384);

% f = linspace(-pi,pi,nsamps);
% z = sqrt((4*f/2).^2-beta^2)+eps;
% wf = abs( 4 .*sin(z)./besseli(0,beta)./z );
% wf = max(wf, 0.3);
% [W1,W2] = meshgrid(wf,wf);
% win_2d = W1.*W2;
% IMG = IMG./win_2d;

figure(91);imagesc(abs(IMG(129:384,129:384))); colormap gray
figure(2); imagesc(abs(image)); colormap gray

figure(3)
subplot(2,1,1)
imm = abs(IMG(nsamps/2,:));
plot(imm(129:384));
subplot(2,1,2)
plot(abs(image(128,:)));
