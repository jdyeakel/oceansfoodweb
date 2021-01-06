# Simulate food webs from xmax
using Distributions
using LinearAlgebra
using DataFrames
using CSV




#import likelihood-maximizing alpha/beta/gamma parameters
#Benguela
xmax = [-1.1776880101309444, 0.412575554787146, -0.023876195960886182];
alpha = xmax[1];
beta = xmax[2];
gamma = xmax[3];

nw_bs = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/bodysizedata.csv",header=true);
nw_fmatrix = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/foraging_matrix.csv",header=true);

nw_minmass = nw_bs[!,Symbol("Mass min")];
nw_maxmass = nw_bs[!,Symbol("Mass max")];
nw_mass = vec(mean([nw_minmass nw_maxmass],dims=2));
trophic = nw_bs[!,Symbol("trophic")];

occur = findall(!iszero,nw_bs[!,Symbol("1700")]);
ns = length(occur);
amatrix = Array{Int64}(undef,ns,ns) .* 0;

#Loop across consumers
for i=1:ns
    if trophic[occur[i]] > 0
        # Loop across prey
        # The while loop means that all consumers will have at least one prey
        while sum(amatrix[:,i]) < 1
            for j=1:ns
                #consumer mass
                mi = nw_mass[occur[i]];
                #prey mass
                mj =nw_mass[occur[j]];

                #The probability of an interaction
                pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));

                #draw
                linkdraw = rand();
                if linkdraw < pij
                    amatrix[j,i] = 1;
                else
                    amatrix[j,i] = 0;
                end
            end
        end
    end
end



