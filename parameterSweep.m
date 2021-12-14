function [results,summary] = parameterSweep(De, p_off, t_total, n_each)

n_runs = length(De) * length(p_off) * n_each;
results(1:n_runs) = struct('De',[],'p_off',[],'ParticleDist',zeros(1,10));
n_conditions = length(De) * length(p_off);
summary(1:n_conditions) = struct('De',[],'p_off',[],'avgParticleDist',zeros(1,10),'pctAnterior',[]);
runcounter = 1;
condcounter = 1;
for d = De
    for p = p_off
        for n = 1:n_each
            fprintf('Now running simulation %d of %d\n',runcounter, n_runs)
            Tracks = MainMatrix_EOD_Dec8th2021(1000,0.05,4,d,7.7,p,t_total,100);
            distr_x = Tracks(:,end, 1);
            binsize = 6;
            for b = 1:10
              loweredge = (b-1)*binsize-30;
              upperedge = b*binsize-30;
              bin_index = [distr_x] > loweredge & [distr_x] < upperedge;
              results(runcounter).ParticleDist(b) = sum(bin_index) / 1000;
            end
            results(runcounter).De = d;
            results(runcounter).p_off = p;
            runcounter = runcounter + 1;
        end
        summary(condcounter).particleDist = mean(cell2mat({results((condcounter-1)*n_each+1:condcounter*n_each).ParticleDist}'),1);
        summary(condcounter).pctAnterior = sum(summary(condcounter).particleDist(1:6));
        summary(condcounter).De = d;
        summary(condcounter).p_off = p;
        condcounter = condcounter+1;
    end
end