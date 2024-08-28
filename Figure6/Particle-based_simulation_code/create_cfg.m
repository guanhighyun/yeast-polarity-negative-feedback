 % Create the directory that stores the simulations
directory = 'Simulations';
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop=4000;
% Store the coordinates of molecules every 10 seconds. Units in 0.1 ms.
samplingrate=100000;
% Range of molecule numbers. Modify if needed.
n_BemGEF = 500; 
n_Cdc42 = 5000;
n_GAP = 0;
P1 = [0.018]; % value of lambda3*dt. Change as needed.
k2 = [0.01,0.015,0.02,0.025,0.03,0.035,0.045,0.055,0.065,0.1,...
    6.85,6.95,7.05,7.15,7.25,7.35,7.45,7.55,7.65,7.75,7.85,7.95,8,8.1,8.3]; % value of k2. Change as needed.
% Number of realizations
random_seeds = 1;

% Create the script file for SLURM to run the jobs.
fid=fopen(sprintf('%s/run.sh',directory),'w');
fprintf(fid,'#!/bin/bash\n\n');
for j = random_seeds
    for u = 1:numel(n_Cdc42)
        for v = 1:numel(n_BemGEF)
            for w = 1:numel(n_GAP)
                for m = 1:numel(P1)
                    for n = 1:numel(k2)
                        curr_seed  = j;
                        curr_fileprefix = sprintf('Cdc42_%g-Bem_%g-GAP_%g-P1_%g-k2_%g-%d',n_Cdc42(u),n_BemGEF(v),n_GAP(w),P1(m),k2(n),j);
                        smoldyn_cfg(curr_fileprefix,directory,tstop,samplingrate,n_Cdc42(u),n_BemGEF(v),n_GAP(w),P1(m),k2(n),curr_seed);
                        cfg_name = sprintf('%s.cfg',curr_fileprefix);
                        fprintf(fid,'sbatch -p general -N 1 -J Smoldyn -t 168:00:00 --mem=8g --wrap="/nas/longleaf/home/kaiyun/Smoldyn/cmake/smoldyn %s"\n',cfg_name);
                    end
                end
            end
        end
    end
end
fclose(fid);