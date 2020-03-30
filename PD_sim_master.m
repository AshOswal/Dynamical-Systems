close all;
clear all;
p.D        = 1;
D = p.D;
p.dTsg     = D*6e-3;   % ms
p.dTgs     = D*6e-3;   % ms
p.dTgg     = D*4e-3;   % ms
p.dTcs     = D*5.5e-3;
p.dTsc     = D*20e-3;  %21.5e-3;
p.dTsgpi   = 6e-3;
p.dTgpegpi = 6e-3;
p.dTgpigpi = 4*1e-3;

p.Ts   = 12.8*1e-3;   % ms
p.Tg   = 20*1e-3;     % ms
p.Te   = 11e-3;       % ms
p.Ti   = 11e-3;       % ms
p.Tgpi = 14e-3;       % see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2585404/#!po=67.1053 for reference


p.Ctx  = 172.18; %27;       % spm/s
p.Str  = 8.46;   %  2;      % sp/s
p.Ms   = 300;               % sp/s
p.Bs   = 10;%17;     % sp/s
p.Mg   = 400;    % sp/s
p.Bg   = 20;%75;     % sp/s
p.Mgi  = 300;
p.Bgi  = 18;
p.Mi   = 205;    % normal range 200-330 sp/s
p.Bi   = 9.87;     % normal range 0-20 sp/s
p.Me   = 75;     % normal range 50-80 sp/s
p.Be   = 17.85;     % normal range 0-20 sp/s

% for the cortical populations
p.dTcs  = D*5.5e-3;   % ms - allowed range 1-10 ms
p.dTcc  = D*6e-3;     %5e-3%1e-3;%4.65e-3;  % 2;     % ms - allowed range 1-10 ms

p.Wsg  = 4.87;   %19 + p.K*(20-19);
p.Wgs  = 1.33;   %1.12 + p.K*(10.7-1.12);
p.Wgg  = 0.53;   %0.53;%6.6 + p.K*(12.3 - 6.6);
p.Wcs  = 0;      %10.98;%2.42 + p.K*(9.2 - 2.42);

% cortical parameters
p.Wsc     = 0;%0.4;%0.4;%0.4%8.93;
p.Wcc     = 4;
p.Wsgpi   = 3.7;
p.Wgpegpi = 1;
%p.axis    = [50 6];
p.S       = 0;
p.SN      = 0;
p.PhaseP  = 0;
p.Method  = 'D';
p.length  = 3;
p.Plot    = 1;

%%
Y = [];
try
    load('FF.mat');
catch
    YY = {};
    FF = {};
    i = 0:0.05:8;
    for n = 1:numel(i)
        p.Wcc = i(n);
        [Y,F] = simulate_PD(p);    
        if ~(p.Wcc==1 || p.Wcc==4)
           close gcf;
        end
        YY{n} = Y;
        FF{n} = F;
    end
    save('FF.mat','FF');
    save('YY.mat','YY')
end
pow   = collect(FF,'s');
powE = squeeze(pow(:,4,:));
powI = squeeze(pow(:,5,:));
figure;
subplot(1,2,1)
imagesc(0:0.05:8,FF{1}.f,powE',[0 6]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(0:0.05:8,FF{1}.f,powI',[0 6]);set(gca,'Ydir','normal','FontSize',12);

%% now keep Wcc fixed at 4 and systematically vary Te and Ti
p.Wcc = 4;
Y1 = [];
try
    load('FF1.mat');
catch
    YY1 = {};
    FF1 = {};
    i = [3:1:20].*1e-3;
    for n = 1:numel(i)
        p.Te = i(n);
        p.Ti = p.Te;
        [Y1,F1] = simulate_PD(p);
        close gcf;
        YY1{n} = Y1;
        FF1{n} = F1;
    end
    save('FF1.mat','FF1');
    save('YY1.mat','YY1')
end
pow   = collect(FF1,'s');
powE = squeeze(pow(:,4,:));
powI = squeeze(pow(:,5,:));
figure;
subplot(1,2,1)
imagesc([3:1:20].*1e-3,FF1{1}.f,powE',[0 12]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc([3:1:20].*1e-3,FF1{1}.f,powI',[0 12]);set(gca,'Ydir','normal','FontSize',12);

%% Finally assess the effect of delays on the frequency response
p.Wcs = 0;
p.Wsc = 0;
p.Te  = 11e-3;
p.Ti  = 11e-3;

try
    load('FF2.mat');
catch
    YY2 = {};
    FF2 = {};
    i = 1e-3:0.5*1e-3:10*1e-3;
    for n = 1:numel(i)
        p.dTcc = i(n);
        [Y2,F2] = simulate_PD(p);
        close gcf;
        YY2{n} = Y2;
        FF2{n} = F2;
    end
    save('FF2.mat','FF2');
    save('YY2.mat','YY2')
end
pow   = collect(FF2,'s');
powE = squeeze(pow(:,4,:));
powI = squeeze(pow(:,5,:));
figure;
subplot(1,2,1)
imagesc(1e-3:0.5*1e-3:10*1e-3,FF{1}.f,powE',[0 8]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(1e-3:0.5*1e-3:10*1e-3,FF1{1}.f,powI',[0 8]);set(gca,'Ydir','normal','FontSize',12);

%% now simulate full model with p.cs & p.sc set to full values

p.Te  = 11e-3;
p.Ti  = 11e-3;
p.Wcc = 4;
p.Wcs = 7;
p.Wsc = 0;
p.Wgs = 1.3;

p.PhaseP = 1;
p.Wsgpi   = 2.7;
simulate_PD(p);
s=findobj('type','legend');
delete(s);

p.Wsc = 0.5;
p.Wsgpi   = 2.7;
simulate_PD(p);
s=findobj('type','legend');
delete(s);

p.Wsc = 1;
p.Wsgpi   = 2.7;
simulate_PD(p);
s=findobj('type','legend');
delete(s);

p.Wsc = 2;
p.Wsgpi   = 2.7;
simulate_PD(p);
s=findobj('type','legend');
delete(s);


%%
p.Wcs = 7.98;
p.Wsc = 0.5;
p.Wsgpi = 2.6;
% simulate the effect of changing the hyperdirect pathway strength
Y3 = [];
try
    load('FF3.mat');
catch
    YY3 = {};
    FF3 = {};
    i = 0:0.5:16;
    for n = 1:numel(i)
        p.Wcs = i(n);
        [Y3,F3] = simulate_PD(p);
        if p.Wcs~=7
        close gcf;
        end
        YY3{n} = Y3;
        FF3{n} = F3;
    end
    save('FF3.mat','FF3');
    save('YY3.mat','YY3')
end
pow   = collect(FF3,'s');
powE = squeeze(pow(:,1,:));
powI = squeeze(pow(:,3,:));
figure;
subplot(1,2,1)
imagesc(0:0.5:16,FF3{1}.f,powE',[0 16]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(0:0.5:16,FF3{1}.f,powI',[0 16]);set(gca,'Ydir','normal','FontSize',12);

%% show the effect of varying cortical time constant on subcortical oscillations
p.Wcs = 7.98;
p.Wsgpi = 2.6;

try
    load('FF4.mat');
catch
    YY4 = {};
    FF4 = {};
    i = [3:1:20].*1e-3;
    for n = 1:numel(i)
        p.Te = i(n);
        p.Ti = i(n);
        [Y4,F4] = simulate_PD(p);
        close gcf;
        YY4{n} = Y4;
        FF4{n} = F4;
    end
    save('FF4.mat','FF4');
    save('YY4.mat','YY4')
end
pow   = collect(FF4,'s');
powE = squeeze(pow(:,1,:));
powI = squeeze(pow(:,3,:));
figure;
subplot(1,2,1)
imagesc([3:1:20].*1e-3,FF4{1}.f,powE',[0 10]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc([3:1:20].*1e-3,FF4{1}.f,powI',[0 10]);set(gca,'Ydir','normal','FontSize',12);

%% now finally the effect of cortical delays on shaping subcortical behaviour
p.Te = 11e-3;
p.Ti = 11e-3;
p.S  = 0;
try
    load('FF5.mat');
catch
    YY5 = {};
    FF5 = {};
    i = 1e-3:0.5*1e-3:8*1e-3;
    for n = 1:numel(i)
        p.dTcc = i(n);
        [Y5,F5] = simulate_PD(p);
        close gcf;
        YY5{n} = Y5;
        FF5{n} = F5;
    end
    save('FF5.mat','FF5');
    save('YY5.mat','YY5');
end
pow   = collect(FF5,'s');
powE = squeeze(pow(:,1,:));
powI = squeeze(pow(:,3,:));
figure;
subplot(1,2,1)
imagesc(1e-3:0.5*1e-3:8*1e-3,FF5{1}.f,powE',[0 10]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(1e-3:0.5*1e-3:8*1e-3,FF5{1}.f,powI',[0 10]);set(gca,'Ydir','normal','FontSize',12);

%%
p.dTcc  = D*6e-3;
p.Wsc   = 0;

try
    load('FF6.mat');
catch
    YY6 = {};
    FF6 = {};
    i = 1e-3:0.5*1e-3:20*1e-3;
    for n = 1:numel(i)
        p.dTgs = i(n);
        p.dTsg = p.dTgs;
        [Y6,F6] = simulate_PD(p);
        %close gcf;
        YY6{n} = Y6;
        FF6{n} = F6;
    end
    save('FF6.mat','FF6');
    save('YY6.mat','YY6');
end

pow   = collect(FF6,'s');
powE = squeeze(pow(:,1,:));
powI = squeeze(pow(:,3,:));
figure;
subplot(1,2,1)
imagesc(1e-3:0.5*1e-3:20*1e-3,FF6{1}.f,powE',[0 16]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(1e-3:0.5*1e-3:20*1e-3,FF6{1}.f,powI',[0 16]);set(gca,'Ydir','normal','FontSize',12);



