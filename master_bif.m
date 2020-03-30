clear all  
p.D        = 1;
D = p.D;
p.dTsg     = 6e-3;     % ms
p.dTgs     = 6e-3;     % ms
p.dTgg     = D*4e-3;   % ms
p.dTcs     = D*5.5e-3;
p.dTsc     = D*20e-3;
p.dTsgpi   = 6e-3;
p.dTgpegpi = 6e-3;
p.dTgpigpi = 4*1e-3;

p.Ts   = 12.8*1e-3;   % ms %12.8
p.Tg   = 20*1e-3;     % ms %20
p.Te   = 11e-3;%10e-3;  % ms
p.Ti   = 11e-3;%10e-3;  % ms
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
p.dTcc  = D*6e-3; %5e-3%1e-3;%4.65e-3;  % 2;     % ms - allowed range 1-10 ms

p.K = 1;
% disease
p.Wsg  = 4.87;%19 + p.K*(20-19);
p.Wgs  = 1.33;%1.12 + p.K*(10.7-1.12);
p.Wgg  = 0.53;%0.53;%6.6 + p.K*(12.3 - 6.6);
p.Wcs  = 8;%2.42 + p.K*(9.2 - 2.42);

% cortical parameters
p.Wsc     = 0;
p.Wcc     = 4;
p.Wsgpi   = 3.7;
p.Wgpegpi = 1;

% stochastic component - GWN
p.S     = 0;
sinput  = 0:2.5:800;
try
    load('bif_data.mat');
    bif = bif_data.max_f_s - bif_data.min_f_s > 1;
    min_f_s = bif_data.min_f_s;
    max_f_s = bif_data.max_f_s;
    min_f_g = bif_data.min_f_g;
    max_f_g = bif_data.max_f_g;
    freq_s  = bif_data.freq_s;
    freq_g  = bif_data.freq_g;
    freq    = bif_data.freq;
catch
    
    for i = 1:numel(sinput)
        p.sinput = sinput(i);
        [Y,F,S] = simulate_PD_bif(p);
        % maximum and minimum firing rates for STN/GP
        min_f_s(i) = min(Y.y(1,Y.x>1.5));
        max_f_s(i) = max(Y.y(1,Y.x>1.5));
        
        min_f_g(i) = min(Y.y(2,Y.x>1.5));
        max_f_g(i) = max(Y.y(2,Y.x>1.5));
        
        [pow,ind]      =  max(F.s,[],2);
        freq_s(i)    =  F.f(ind(1));
        freq_g(i)    =  F.f(ind(2));
        
    end
    
    bif_data = [];
    bif_data.min_f_s = min_f_s;
    bif_data.max_f_s = max_f_s;
    bif_data.min_f_g = min_f_g;
    bif_data.max_f_g = max_f_g;
    bif_data.freq_s  = freq_s;
    bif_data.freq_g  = freq_g;
    save('bif_data.mat','bif_data');
    bif = max_f_s - min_f_s > 1;
end

x   = sinput(~bif);
%%
figure;
subplot(2,1,1)
plot(x,min_f_s(~bif),'k.');hold on;
scatter(sinput(bif),max_f_s(bif),8,freq_s(bif),'filled');hold on;
scatter(sinput(bif),min_f_s(bif),8,freq_s(bif),'filled');hold on;
colorbar
xlabel('Cortical input (spk/s)');
ylabel('STN firing rate (spk/s)');
set(gca,'FontSize',12);
subplot(2,1,2)
plot(x,min_f_g(~bif),'k.');hold on;
scatter(sinput(bif),max_f_g(bif),8,freq_g(bif),'filled');hold on;
scatter(sinput(bif),min_f_g(bif),8,freq_g(bif),'filled');hold on;
colorbar
xlabel('Cortical input (spk/s)');
ylabel('GPe firing rate (spk/s)');
set(gca,'FontSize',12);
