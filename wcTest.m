close all;
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
p.Wcs  = 4;      %10.98;%2.42 + p.K*(9.2 - 2.42);

% cortical parameters
p.Wsc     = 0;%0.4;%0.4;%0.4%8.93;
p.Wcc     = 4;
p.Wsgpi   = 3.7;
p.Wgpegpi = 1;
%p.axis    = [50 6];
p.S       = 0;
%%
Y = [];

YY = {};
FF = {};
i = 1:0.05:4.5;
for n = 1:numel(i)
    p.Wcc = i(n);
    [Y,F] = simulate_PD(p);
    YY{n} = Y;
    FF{n} = F;
end

pow   = collect(FF,'s');
powE = squeeze(pow(:,4,:));
powI = squeeze(pow(:,5,:));
figure;
subplot(1,2,1)
imagesc(3:0.05:4.5,FF{1}.f,powE',[0 12]);set(gca,'Ydir','normal','FontSize',12);
subplot(1,2,2)
imagesc(3:0.05:4.5,FF{1}.f,powI',[0 12]);set(gca,'Ydir','normal','FontSize',12);
