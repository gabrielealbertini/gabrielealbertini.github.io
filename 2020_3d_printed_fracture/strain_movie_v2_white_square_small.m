clear; %close all ; clc
addpath(fullfile('/home/ga288/sources/optio/solidmechanics/lefm'));
addpath(fullfile('/home/ga288/sources/optio/solidmechanics/linearelasticity'));
%addpath('../crc_FIGSM');
%addpath('../util/TLS');
%addpath('../util');
dictoolpath = '/home/ga288/sources/dic_tools/analysis_compare';
% this adds the correct folder to the matlab path
addpath(fullfile(dictoolpath));

%% fig setup
addpath(fullfile('/home/ga288/sources/matlab-utilities/altmany-export_fig'));
%axis ticklabels: Arial, size=6
fontsize = 12;
axesfontsize = fontsize;

% subplot letter a,b,c ... Times New Roman, bold, size=12

% axis titles - Cambria Math italic size=8 (but this might be
% windows only)

markersize = 12;
Ncolor = 21;

set(0,'DefaultTextFontsize',fontsize);
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',axesfontsize);
fontfct=1.3
set(0,'DefaultAxesTitleFontSizeMultiplier',1)%6/4)%2)
set(0,'DefaultAxesLabelFontSizeMultiplier',fontfct)%6/4)%5/4)

set(0,'DefaultLineMarkerSize',markersize);

%'defaultLineLineWidth'

%set(0,'defaultAxesLabelInterpreter','latex');
set(0,'defaultTextInterpreter','tex');%'latex');
set(0,'defaultAxesTickLabelInterpreter','tex');
set(0, 'DefaultLegendInterpreter', 'tex')%'latex')
set(0, 'DefaultAxesBox','on');
set(0, 'DefaultLegendBox','off');

set(0,'DefaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'DefaultPatchFaceColor','none');
set(0,'DefaultPatchMarkerSize',markersize);
set(0,'DefaultFigureColormap',jet);
% =======================================================

gray=[169,169,169]./256;%192,192,192

matcolororder =[  0.0000,    0.4470,    0.7410;
                  0.8500,    0.3250,    0.0980;
                  0.9290,    0.6940,    0.1250;
                  0.4940,    0.1840,    0.5560;
                  0.4660,    0.6740,    0.1880;
                  0.3010,    0.7450,    0.9330;
                  0.6350,    0.0780,    0.1840];

%drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:});

%------------------------------------------------

bname=split(pwd,'/');
bname=bname{end-1}
% setting some default plot options
%dpi = 200;
%fontsize = 12;
markersize = 12;
Ncolor = 21;

CorreliPath = '/home/ga288/sources/correli/';

% this adds the correct folders to the matlab path
addpath(fullfile(CorreliPath,'ml_lib_core'));
exprtpath='../';

% Read the image series
% ======================================================

load([exprtpath,'DIC_images.mat'],'images');
images=images.crack;
imfin = numel(images);


flip=false;%true;%false;
if flip 
    imgidx=fliplr(imgidx);
end

% Pixel size (meter)
%load('pixel_size.mat'); % PARAMETER TO ENTER
%pix2m = pixelsize;
load([exprtpath,'/pixel_size.mat']); % PARAMETER TO ENTER
pix2m = pixelsize;


% reference image
refim=1%2000
fxmin=100;

f = imread(fullfile(images{refim}));
f = f(:,100:end);
[n, m] = size(f);

% blur images to increase autocorelation length

blur = 0; % std of gaussian window

%f = image_blur(f,blur);

load_results_exist=true;%
sourcefile='single_pacman_reslts_crack_path_v0.mat'%'single_pacman_reslts.mat'
load(sourcefile,'results');

acquisition_fps=232.0e3;
dfp=232

%imgidx=results.imgidx;
imgidx = [%dfp:dfp:dfp*int64(14000/dfp),...
          dfp*int64(14100/dfp+1):dfp/2:dfp*int64(27000/dfp),...
          27050:  1:27090-6];%,...
if numel(imgidx)~=numel(unique(imgidx))
    error('imgidx must be unique')
end


%------------------------------
% load analysis file
load('../analysis/williams_compliance.mat','Young')

Young_will=2.5e9;
crack_0=0.01950;
delta_t= 1.0 /acquisition_fps;
results.time = double(results.imgidx).*delta_t;
results.crack_length=(results.a-results.a(1))*pix2m+crack_0
results.speed=compute_speed(results.crack_length,results.time,'back_diff');%'pchip');%
results.K1 = results.K1*1e6/Young_will*Young
results.K2 = results.K2*1e6/Young_will*Young
results.G = (results.K1.^2 + results.K2.^2)./Young;
stripe_start=0.03;
stripe_width=split(pwd,'Stripes');
stripe_width=split(stripe_width{2},'_');
stripe_width=str2num(stripe_width{2})*1e-3;
results.is_stripe =  @(x) (0.5*sign(mod(-x+0.0005+stripe_start,stripe_width*2)-stripe_width)+0.5).*...
    (0.5*sign(x-stripe_start)+0.5);
results.is_color = results.is_stripe(results.crack_length)'* matcolororder(2,:) + ...
 (1-results.is_stripe(results.crack_length))'* matcolororder(1,:);

%======================================
% exprt: adjust E
%======================================
load('/home/ga288/research/mode_i_heterogeneities_ponson/experiments_summary/crc_FIGSM/crc_mat.mat','mat_consolidated','exprtVC','exprtVWDW')

if true%false%for widx=1:length(will)
    for idx=1:length(exprtVC)
        if strcmp(exprtVC(idx).bname,bname)
            vcidx=idx;
            break
        end
    end
    
    for idx=1:length(exprtVWDW)
        if strcmp(exprtVWDW(idx).bname,bname)
            vwdwidx=idx;
            break
        end
    end

    l=results.crack_length; 
    l_up=min(l):1e-4:max(l);
    is_obst_up=results.is_stripe(l_up);
    a=[1];
    b=ones(1,20);
    b=b/length(b);
    is_obst_up=filtfilt(b,a,is_obst_up);
    
    is_obst=interp1(l_up,is_obst_up,l);
    
    results.G= results.G.*((1-is_obst).*exprtVC(vcidx).youngfct+...
                                 is_obst.*exprtVWDW(vwdwidx).youngfct*1.1);
end
% Crack definition
% ======================================================

load([exprtpath,'/crack_path.mat']);%crack_path_refim.mat'])
x=P(:,2)-2;   
y=P(:,1)-fxmin;
%v = flipud([y x]);%
v = [y x];

if false%true%
    figure('Position',[0 500 650 400])
    colormap('gray')
    imagesc(f)
    daspect([1 1 1]);
    colorbar
    hold on
    plot(y,x)
    
end    
crackpath = crack_makepath(v);

% initial crack length
a0 = min(y); %5.0e-3 /pix2m;

% initialize the crack
crackpos = crack_position(crackpath,a0);
crackang = crack_angle(crackpath,a0);


%%
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesTickLabelInterpreter','tex');
%fig1 = figure('Units','inches','Position',[0,0,3.375*1.1,2*1.1])
f=imread(images{1});
    
fxlim=[fxmin,900]
f=f(:,fxlim(1):fxlim(2));
[ny,nx]=size(f);

fig_fct=0.75
xpad=20
ypad=290
fig1=figure('Position',[0,0,(nx+xpad)*fig_fct,(ny*3+154+ypad)*fig_fct])

set(fig1,'Color','w')
axraw = subplot(11,1,1);%raw data
axgl =   subplot(11,1,2);%global dic
axw =   subplot(11,1,3);% williams
axwf =  subplot(11,1,4);
axglf =  subplot(11,1,5);
cE =    subplot(11,1,6);
axmm =  subplot(11,1,7);
axv_l  =  subplot(11,1,8);
axGvsl  =  subplot(11,1,9);
axv_lf =  subplot(11,1,10);
axGvslf =  subplot(11,1,11);

figxfct=(nx+xpad)/(nx+xpad);
figyfct=(ny*3+154)/(ny*3+154+ypad)
figyshift=(1-figyfct)/1.02;
ploty=ny/(ny*3+154);
timepos=[0.0,1.03,0]
cE.Position=   [0.07500*figxfct -0.2500*figyfct+figyshift 0.02500*figxfct 0.1500*figyfct];

axraw.Position=[0.000*figxfct 0.6700*figyfct+figyshift figxfct ploty];
mm5=axraw.Position(3)/(nx*pix2m)*0.005;
axmm.Position =[0.200*figxfct .96500*figyfct+figyshift mm5 0.003*figyfct];

axgl.Position=  axraw.Position;
axgl.Position(2)= 0.3500*figyfct+figyshift;
axw.Position=  axraw.Position;
axw.Position(2)=  0.0300*figyfct+figyshift;

axwf.Position = axw.Position;
axglf.Position = axgl.Position;

axv_l.Position =  [0.25*figxfct figyshift*0.3 figxfct*0.3, figyshift*0.6];
axGvsl.Position = [0.675*figxfct figyshift*0.3 figxfct*0.3, figyshift*0.6];
axv_lf.Position = axv_l.Position
axGvslf.Position = axGvsl.Position

cmap=colormap('gray');
jcmap=colormap('jet');

textcolor='k';%[1,1,1];%cmap(60,:);
rawlim=[0,3200]

xscale=1e3;
epsscale=1e3;


set(fig1,'CurrentAxes',axv_l)
%set(axv_l,'YColor','w')
%set(axv_l,'XColor','w')
%set(axv_l,'XTickLabel',[]);
xlabel('{\it{l}} (mm)')

hold  on
ylabel('{\it{v}} (m/s)')
scattxl=[0.015,0.06].*xscale;
xlim(scattxl)
vlim=[0,500];
ylim(vlim)

text(scattxl(1),vlim(2)*1.15,{'crack speed: {\it{v}}'}','fontsize',fontsize*1.1);%,'color','w
set(fig1,'CurrentAxes',axGvsl)

%set(axGvsl,'YColor','w')
%set(axGvsl,'XColor','w')
hold  on
ylabel('{{\bf{\Gamma}}} (J/m^2)')
xlabel('{\it{l}} (mm)')
xlim(scattxl)
Glim=[0,450];
ylim(Glim)
text(scattxl(1),Glim(2)*1.15,{'fracture energy: {{\bf{\Gamma}}}'}', ...
     'fontsize',fontsize*1.1);%,'color','w

set(fig1,'CurrentAxes',axv_lf)
hold on
set(axv_lf,'color','none');
set(axv_lf,'XTickLabel',[]);
set(axv_lf,'YTickLabel',[]);
xlim(scattxl)
ylim(vlim)

set(fig1,'CurrentAxes',axGvslf)
hold on
text(scattxl(2)*-1.45,Glim(2)*-0.4,'(Albertini et al., 2020)')%,'color','w')
set(axGvslf,'color','none');
set(axGvslf,'XTickLabel',[]);
set(axGvslf,'YTickLabel',[]);
xlim(scattxl)
ylim(Glim)
scatter([],[],markersize,matcolororder(1,:))
scatter([],[],markersize,matcolororder(2,:))
legend({'matrix','obstacle'})
%set(axGvslf,'color','none','visible','off');

set(fig1, 'CurrentAxes',axmm)
axis off
hold on
title(axmm, '5 mm','color',textcolor)
X=[0,1;
   0,1];
Y=[0,0;
   1,1];
pcolor(Y,X,[1,1;1,1]);
colormap(axmm,[0,0,0;0,0,0])

shading flat

set(fig1,'CurrentAxes',cE)
hold on
axis off
colormap(gca,'Jet')


clm=[0, 0.007];

crange=0:0.01:clm(2)*epsscale;
X=[crange;
   crange];
Y=[zeros(size(crange));
   ones(size(crange))];
pcolor(Y,X,X);
shading flat
a=[0,clm(2)/2,clm(2)].*epsscale;

for i=1:length(a)
    text(-0.5,a(i),num2str(a(i)),'color',textcolor,'interpreter','tex','HorizontalAlignment','right');
end
text(0.5,8/5*clm(2)*epsscale,'\epsilon_{\it{yy}} (10^{-3})','color',textcolor,...
     'interpreter','tex','HorizontalAlignment','center','fontsize',fontsize*fontfct)



fname_noloop = [bname,'_strain_images_v2_noloop_white_square_small.gif'];
fname = [bname,'_strain_images_v2_white_square_small.gif'];

wtime=1/12;

[~,i_start]=min(abs(results.imgidx-imgidx(1)));
for idx=imgidx(1:end-1)%i = 1:(length(imgidx)-7)
    %idx=imgidx(i);
    [~,i]=min(abs(results.imgidx-idx));
    didx=results.imgidx(i+1)-results.imgidx(i);
    
    fprintf('\nFrame %d --------------------\n',idx);
    defim = idx;

    crackpos = crack_position(crackpath,results.a(i));%-[100,0];
    xlm=[crackpos(1)-175,crackpos(1)+175];
    ylm=[0,ny];

    g = imread(images{defim});
    g = g(:,fxlim(1):fxlim(2));
    [ny,nx]=size(g);
    g(:,1:int64(xlm(1)))  =g(:,1:int64(xlm(1))).*0.4;
    g(:,int64(xlm(2)):nx)=g(:,int64(xlm(2)):nx).*0.4;

    % plot fig and box limits
    set(fig1,'CurrentAxes',axraw)
    cla(axraw)
    hold on
    axis off
    imagesc(3280-g,rawlim)
    colormap(gca,'gray')
    daspect([1,1,1])
    
    xlim([0,nx-1])
    ylim([0,ny])
    
    scatter(crackpos(1),crackpos(2),500,jcmap(60,:),'.')

    if false
        plot(xlm([1,1,2,2,1]),ylm([1,2,2,1,1]),'k') % box
        plot(y,x,'w') % crack path
    end

    title(strcat(num2str(round(double(idx)/acquisition_fps,6),'%6.6f s')),...%'frame #  ', string(i)),...
          'color',textcolor,'Units','normalized','interpreter','tex',...
          'Position',timepos,'HorizontalAlignment','left','VerticalAlignment','bottom')
   
    %--------------------------------------------
    % plot strain from global analysis
    set(fig1,'CurrentAxes',axgl)
    hold on
    axis off
    cla(axgl)
    
    colormap(axgl,'Jet')
    mesh_plot(results.meshw{i},...%Mesh_ini,...
              'Field',results.Eg{i}(:,2).*epsscale,...
              'Skin',false,'Wireframe',false, ...
              'Edgecolor','none','Pixelsize',pix2m*xscale);
    axgl.CLim=clm.*epsscale;
    daspect([1,1,1])

    xlim([0,nx-1]*pix2m*xscale)%xlim(xlm*pix2m*xscale)
    ylim(ylm*pix2m*xscale)    

    
    %--------------------------------------------
    % plot strain from williams integrated analysis
 
    % background image
    set(fig1,'CurrentAxes',axw)
    axis off
    cla(axw)
    colormap(axw,'Jet')
    p1=mesh_plot(results.meshg{i},...%D(1).Mesh,...
              'Field',results.Ew{i}(:,2).*epsscale,...
              'Skin',false,'Wireframe',false, ...
              'Edgecolor','none','Pixelsize',pix2m*xscale);
    axw.CLim=clm.*epsscale;
    daspect([1 1 1]);
        
    xlim([0,nx-1]*pix2m*xscale)%xlim(xlm*pix2m*xscale)
    ylim(ylm*pix2m*xscale) 
    %-------------------------------------
    % frame
    
    xfrm=[0:1:nx].*pix2m*xscale;
    yfrm=[0:1:ny].*pix2m*xscale;
    [Xfrm,Yfrm]=meshgrid(xfrm,yfrm);
    Cfrm=zeros(size(Xfrm));
    
    Cfrm(:,max(1,int64(xlm(1))):min(int64(xlm(2)),nx))=NaN;
    
    % foreground image   
    set(fig1,'CurrentAxes',axglf)
    axis off
    cla(axglf)
    pcolor(Xfrm,Yfrm,Cfrm)%,'alphadata',Cfrm>1e-2)
    shading flat
    colormap(axglf,[1,1,1; 1,1,1;1,1,1])
    set(axglf,'color','none','visible','off');
    %xlim([0,nx-1]*pix2m*xscale)%xlim(xlm*pix2m*xscale)
    ylim(ylm*pix2m*xscale) 
    %linkaxes([axgl axglf])
    
    set(fig1,'CurrentAxes',axwf)
    axis off
    cla(axwf)
    pcolor(Xfrm,Yfrm,Cfrm)%,'alphadata',Cfrm>1e-2)
    shading flat
    colormap(axwf,[1,1,1; 1,1,1;1,1,1])
    set(axwf,'color','none','visible','off');
    xlim([0,nx-1]*pix2m*xscale)%xlim(xlm*pix2m*xscale)
    ylim(ylm*pix2m*xscale) 
    %linkaxes([axw axwf])
    
    
    %-------------------------------------
    % plots scatter v(l) G(l)
    set(fig1,'CurrentAxes',axv_l)
    hold on

    scatter(results.crack_length(i_start:i)*xscale,results.speed(i_start:i),markersize, ...
            results.is_color(i_start:i,:))
    
    set(fig1,'CurrentAxes',axGvsl)
    hold on
    scatter(results.crack_length(i_start:i)*xscale,results.G(i_start:i),markersize, ...
            results.is_color(i_start:i,:))
    
    
    %-------------------------------------
    % save to .gif or gif
    drawnow

    frame = getframe(fig1);
    im=frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if idx==imgidx(1)
        imwrite(A,map,fname_noloop,'gif','LoopCount',1,'DelayTime',wtime);
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',wtime);
    else
        if didx==1
            wt=wtime*9;
        else
            wt=wtime;
        end
        imwrite(A,map,fname_noloop,'gif','WriteMode','append','DelayTime',wt);
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',wt);
    end
end


% last images
cla(axw)
cla(axgl)

for idx=imgidx(end):23:29000

    g = imread(images{idx});
    g = g(:,fxlim(1):fxlim(2));
    [ny,nx]=size(g);
    %g=g.*0.4;

    % plot fig and box limits
    set(fig1,'CurrentAxes',axraw)
    cla(axraw)
    hold on
    axis off
    imagesc(3280-g.*0.4,rawlim)
    colormap(gca,'gray')
    daspect([1,1,1])

    xlim([0,nx-1])
    ylim([0,ny])
    
    title(strcat(num2str(round(double(idx)/acquisition_fps,6),'%6.6f s')),...%'frame #  ', string(i)),...
          'color',textcolor,'Units','normalized','interpreter','tex',...
          'Position',timepos,'HorizontalAlignment','left','VerticalAlignment','bottom')
   
    drawnow
    wt=wtime;   
    
    frame = getframe(fig1);
    im=frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',wt);
    imwrite(A,map,fname_noloop,'gif','WriteMode','append','DelayTime',wt);
end