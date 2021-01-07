using Distributed
using SharedArrays
using Distributions
using CSV
using RCall
using Optim
using Plots
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/logitlikelihood.jl")
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nichemodelweb.jl")
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/lfunc.jl")
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/fc.jl")
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nullweb.jl")
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/bsfoodweb.jl")

#NOTE: 10/28/20 - Build in rectangular matrix; calculate %correct and %incorrect for both links and nonlinks

#Load food web interaction matrix
# Aload = ;
# massload = ;

Adata_bg = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BenguelaSea/bengula_A_matrix.csv",header=true);
numsp = size(Adata_bg)[1];
massdata_bg = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BenguelaSea/bengula_masses.csv",header=true);

## BENGUELA
#Use species 5-29 (species 1-4 have missing diet data)
Aprime_bg = Matrix(Adata_bg[1:numsp,2:numsp+1]);
# nolinks = findall(iszero,sum(Aprime_bg,dims=1));
# A_bg = Matrix(Adata[5:29,6:30]);
A_bg = copy(Aprime_bg);
# massvec_bg = massdata[:mass][5:numsp];
massvec_bg = copy(vec(massdata_bg[!,:mass_kg]));


#######################
# FITTING
#######################



## BENGUELA
x0 = [0.0,0.0,0.0];
results_bg = optimize(x->lfunc(x,A_bg,massvec_bg),x0,NelderMead());
results_bg.minimizer
xmax_bg = results_bg.minimizer;
fcorr_bg, Apredict_bg = fc(xmax_bg,A_bg,massvec_bg)

x0 = [0.0,0.0,0.0];
Anull_bg = nullweb(A_bg);
resultsnull_bg = optimize(x->lfunc(x,Anull_bg,massvec_bg),x0,NelderMead());
resultsnull_bg.minimizer
xmax = resultsnull_bg.minimizer;
fcorrnull_bg, Apredictnull_bg = fc(xmax,Anull_bg,massvec_bg)

sortsp_bg = sortperm(massvec_bg);
Asort_bg = A_bg[sortsp_bg,sortsp_bg];

# use xmax to build a random food web

nw_bs = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/bodysizedata.csv",header=true);
nw_fmatrix = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/foraging_matrix.csv",header=true);
nw_minmass = nw_bs[!,Symbol("Mass min")];
nw_maxmass = nw_bs[!,Symbol("Mass max")];
nw_mass = vec(mean([nw_minmass nw_maxmass],dims=2));
trophic = nw_bs[!,Symbol("trophic")];
date = "1700";
occur = findall(!iszero,nw_bs[!,Symbol(date)]);
bgmatrix_1700 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);

date = "2000";
occur = findall(!iszero,nw_bs[!,Symbol(date)]);
bgmatrix_2000 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);


preymasses = collect(1:10:1000);
predmasses = collect(1:10:1000);
lm = length(preymasses);
pijarray = Array{Float64}(undef,lm,lm);
alpha = xmax_bg[1];
beta = xmax_bg[2];
gamma = xmax_bg[3];
for i=1:lm
    mi = predmasses[i];
    for j=1:lm
        mj = preymasses[j];
        pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));
        pijarray[i,j] = pij;
    end
end





# p1=heatmap(Array{Int64}(Asort_bg),zlims=(0,1));
# p2=heatmap(Array{Int64}(Asort_pc),zlims=(0,1));
# p3=heatmap(Array{Int64}(Asort_bts),zlims=(0,1));
# p4=heatmap(Array{Int64}(Anull_bts),zlims=(0,1));
# plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
# savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/fig_cons.pdf")


# sortsp = sortperm(massvec_bts);
# Asort_bts = A_bts[sortsp,sortsp];
# Anull_bts = Anull_bts[sortsp,sortsp];
# p1=heatmap(Array{Int64}(Asort_bts),zlims=(0,1));
# p2=heatmap(Array{Int64}(Anull_bts),zlims=(0,1));
# p3=heatmap(Array{Float64}(Apredict_bts[sortsp,sortsp]),zlims=(0,1));
# p4=heatmap(Array{Float64}(Apredictnull_bts[sortsp,sortsp]),zlims=(0,1));
# plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
# savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/Anull_bts_cons.pdf")


# sortsp = sortperm(massvec_bg);
# Asort_bg = A_bg[sortsp,sortsp];
# Anull_bg = Anull_bg[sortsp,sortsp];
# p1=heatmap(Array{Int64}(Asort_bg),zlims=(0,1));
# p2=heatmap(Array{Int64}(Anull_bg),zlims=(0,1));
# p3=heatmap(Array{Float64}(Apredict_bg[sortsp,sortsp]),zlims=(0,1));
# p4=heatmap(Array{Float64}(Apredictnull_bg[sortsp,sortsp]),zlims=(0,1));
# plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
# savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/Anull_bg_cons.pdf")

# sortsp = sortperm(massvec_pc);
# Asort_pc = A_pc[sortsp,sortsp];
# Anull_pc = Anull_pc[sortsp,sortsp];
# p1=heatmap(Array{Int64}(Asort_pc),zlims=(0,1));
# p2=heatmap(Array{Int64}(Anull_pc),zlims=(0,1));
# p3=heatmap(Array{Float64}(Apredict_pc[sortsp,sortsp]),zlims=(0,1));
# p4=heatmap(Array{Float64}(Apredictnull_pc[sortsp,sortsp]),zlims=(0,1));
# plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))



#Calculate trophic level
R"library(NetIndices)"
R"ctl_1700 <- TrophInd($(bgmatrix_1700))[,1]"
@rget ctl_1700;
ctlcolors_1700 = Int64.(round.(ctl_1700 .* 10));
ctlcolors_1700 = ctlcolors_1700 .- minimum(ctlcolors_1700) .+ 1;

R"ctl_2000 <- TrophInd($(bgmatrix_2000))[,1]"
@rget ctl_2000;
ctlcolors_2000 = Int64.(round.(ctl_2000 .* 10));
ctlcolors_2000 = ctlcolors_2000 .- minimum(ctlcolors_2000) .+ 1;

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panel.pdf";
R"""
library(igraph)
library(plot.matrix)
library(RColorBrewer)
pdf($namespace,width=14,height=12)
par(mfrow=c(2,2))
pijpal = brewer.pal(9,'YlGnBu')

plot($(pijarray),main='',border=NA,col=pijpal,xlab='Pred mass (kg)',ylab='Prey mass (kg)')

plot(1-$(Apredict_bg[sortsp_bg,sortsp_bg]),axis.col=3,main='')

# plot($(bgmatrix_1700),axis.col=3,main='',key=NULL)
pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max(c($ctlcolors_1700,$ctlcolors_2000)))
agraph = graph_from_adjacency_matrix($(bgmatrix_1700))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_1700
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_1700],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white')

# plot($(bgmatrix_2000),axis.col=3,main='',key=NULL)
# pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max($ctlcolors_2000))
agraph = graph_from_adjacency_matrix($(bgmatrix_2000))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_2000
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_2000],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white')
dev.off()
"""
