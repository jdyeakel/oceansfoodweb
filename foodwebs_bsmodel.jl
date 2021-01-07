# Simulate food webs from xmax
using Distributions
using LinearAlgebra
using DataFrames
using CSV
using RCall



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
nc = length(findall(x->x>1,trophic));
amatrix = Array{Int64}(undef,ns,ns) .* 0;

#Loop across consumers
for i=1:ns
    if trophic[occur[i]] > 1
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

#sort by body mass
bm_sortperm = sortperm(nw_mass[occur]);
amatrix_sort = amatrix[bm_sortperm,bm_sortperm];

#Calculate trophic level
R"library(NetIndices)"
R"ctl <- TrophInd($(amatrix_sort))[,1]"
@rget ctl;
ctlcolors = Int64.(round.(ctl .* 10));
ctlcolors = ctlcolors .- minimum(ctlcolors) .+ 1;

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_rough.pdf";
R"""
library(igraph)
library(plot.matrix)
library(RColorBrewer)
pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max($ctlcolors))
agraph = graph_from_adjacency_matrix($amatrix_sort)
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl
pdf($namespace,width=12,height=6)
par(mfrow=c(1,2))
plot($amatrix_sort,axis.col=3,main='',key=NULL)
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white')
dev.off()
"""
