function [] = RunDir(dirstring)

try 
    matlabpool 
catch
    matlabpool close 
    matlabpool
end



addpath('/home/mleite/Documents/MATLAB/my_toolbox')
addpath(genpath('/home/mleite/Documents/MATLAB/spm12'))
addpath /home/mleite/Documents/MATLAB/FokkerPlanck_Special 

TestParam5
pathstring = ['/home/mleite/Documents/MATLAB/FokkerPlanck_Special/parallel/' dirstring];

cd(pathstring)

mkdir IterVar
FmaxP = -Inf;
FmaxF = -Inf;
F0   = zeros(20,1);
k=1;
while k
    
    disp('*************')
    disp(['*** ' dirstring '***     - Iteration:' num2str(k)])
    disp('*************')
    
    try        
        [Ep,Cp,Eh,F,dFdp,dFdpp] = spm_nlsi_GN_LC(M,0,Y);
    end
    
    load( [pathstring '/IterVar/EpMax'] )
    M.P = Ep;
    
end





