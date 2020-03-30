%% Genetic algorithm

%% Initialize data

TestParam8CA1spinit
mkdir pop
mkdir jobwait
mkdir jobdone
mkdir jobrequest
mkdir jobrequestdenied
mkdir jobprogress
mkdir generations

delete jobwait/*
delete jobdone/*
delete jobrequest/*
delete jobrequestdenied/*
delete jobprogress/*

%% Initialize population

fprintf('Initializing Population...\n')
ngenerations = 4;
npop = 10; % number of individuals
nmates = 1;b = 0;
while b <= npop
    nmates = nmates+1;
    b = nchoosek(nmates,2) + nmates;
end
nmates = nmates-1;

population  = {};
populationi = {};
populationf = zeros(npop,1);

% Decrease variance on connectivity matrices to avoid explosions
Cp = M.pC;
Cp.ang= exp(-32);
Cp.A = spm_unvec(exp(log(spm_vec(Cp.A)) - 4),Cp.A);
Cp.C = spm_unvec(exp(log(spm_vec(Cp.C)) - 4),Cp.C);
Cp =  .25*diag(spm_vec(Cp));
MixM = chol(Cp);
np   = length(MixM);


for jobid = 1:npop
    
    M.P = spm_unvec( spm_vec(M.pE) + MixM*rand(np,1), M.pE);
%     
%     M.P.A(1).M(1,2) = log(normY(4)/normCE);
%     M.P.A(2).M(2,1) = log(normY(1)/normCI);
%     M.P.A(3).M(2,2) = log(normY(2)/normCE);
%     M.P.A(4).M(1,1) = log(normY(3)/normCI);
%     
    filename= ['pop/M' num2str(jobid)];
    save(filename,'M','Y','jobid')
    
end



%% Optimize individuals

fprintf('Optimizing first generation...\n')
movefile pop/* jobwait/

waitlist     = dir('jobwait');
progresslist = dir('jobprogress');
jobprogress  = [];

while length(waitlist)>2
    
    % check for requests
    requestlist = dir('jobrequest');
    
    if length(requestlist)>2
        
        load(['jobrequest/' requestlist(3).name])
        
        if sum(jobprogress == jobid)==0
            jobprogress = [jobprogress jobid];
            movefile(['jobrequest/' requestlist(3).name],...
                ['jobprogress/' requestlist(3).name])
            delete(['jobwait/' requestlist(3).name])
            fprintf('Job %i accpted for machine %i \n', jobid, machineID)
        else
            movefile(['jobrequest/' requestlist(3).name],...
                ['jobrequestdenied/' requestlist(3).name])
            fprintf('Job %i denied for machine %i\n', jobid, machineID)
        end
        
        % eliminage repeated requests
        for k = 4:length(requestlist)
            load(['jobrequest/' requestlist(k).name])
            if sum(jobprogress == jobid)>0
                movefile(['jobrequest/' requestlist(k).name],...
                         ['jobrequestdenied/' requestlist(k).name])
                fprintf('Job %i denied for machine %i (b) \n', jobid, machineID)
            end
        end
        
    end
    
    waitlist     = dir('jobwait');
    pause(2)

end

%wait for jobs to finish
fprintf('Waiting for jobs to finish\n')
donelist = dir('jobdone');
while length(donelist)<npop+2
    fprintf('.')
    pause(20)
    donelist = dir('jobdone');
end
    
%% While not convergence
mkdir generations/0
copyfile jobdone/* generations/0/
%%
% Mate the most apt individuals
for iter = 1:ngenerations
    % wait for jobs to complete
    
    mkdir(['generations/' num2str(iter)])
    copyfile('jobdone/*', ['generations/' num2str(iter) '/'])
    
    fprintf('\nGenerating offspring %i ...\n', iter)
    donelist = dir('jobdone');
    
    for k = 1:npop
        load(['jobdone/' donelist(k+2).name])
        populationi{k} = M;
        populationf(k)   = M.F;
    end
    
    % sort population according to fitness
    [populationf,IX] = sort(populationf,'descend');
    for k = 1:npop
        population{k} = populationi{IX(k)};
    end
    
    % Generate mating pairs
    matingpairs = nchoosek(1:nmates,2);
    
    for k = 1:length(matingpairs)
        M1 = population{matingpairs(k,1)};
        M2 = population{matingpairs(k,2)};
        M  = mate(M1,M2);
%         M3 = M;
%         M3.Ep = M3.pE;
%         M3.Cp = diag(spm_vec(M3.pC));
%         M = mate(M,M3);
        M.converged = 0;

        population{k+nmates} = M;
    end
    
    % Optimize offspring
    
    fprintf('Optimizing offspring %i ...\n', iter)
    
    delete pop/*
    delete jobwait/*
    delete jobdone/*
    delete jobrequest/*
    delete jobrequestdenied/*
    delete jobprogress/*
    
    for jobid = 1:npop
        M = population{jobid};
        if  M.converged
            filename= ['jobdone/M' num2str(jobid)];            
        else
            filename= ['jobwait/M' num2str(jobid)];
        end
        save(filename,'M','Y','jobid')
    end
           
    waitlist     = dir('jobwait');
    progresslist = dir('jobprogress');
    jobprogress  = [];
    fprintf('Distributing jobs\n')
    while length(waitlist)>2
        
        % check for requests
        requestlist = dir('jobrequest');
        
        if length(requestlist)>2
            try
                load(['jobrequest/' requestlist(3).name])
            end
            if sum(jobprogress == jobid) == 0
                jobprogress = [jobprogress jobid];
                movefile(['jobrequest/' requestlist(3).name],...
                    ['jobprogress/' requestlist(3).name])
                delete(['jobwait/' requestlist(3).name])
                fprintf('Job %i accpted for machine %i \n', jobid, machineID)
            else
                movefile(['jobrequest/' requestlist(3).name],...
                    ['jobrequestdenied/' requestlist(3).name])
                fprintf('Job %i denied for machine %i \n', jobid, machineID)
            end
            
        end
        
        waitlist     = dir('jobwait');
        progresslist = dir('jobprogress');
        pause(2)
        
    end
    
    %wait for jobs to finish    
    fprintf('Waiting for jobs to finish\n')
    donelist = dir('jobdone');
    while length(donelist)<npop+2
        fprintf('.')
        pause(20)
        donelist = dir('jobdone');
    end
end
    
    fprintf('\nFINISHED\n')
    
    
    
    
    
    
    
    
    


