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
#Sort by body size
sortsp_bg = sortperm(massvec_bg);
Asort_bg = A_bg[sortsp_bg,sortsp_bg];

#Import Potters Cove and Barrents Sea adjacency matrices
Adata_pc = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/adjmatrix.csv",header=true);
massdata_pc = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/masslist.csv",header=true);
## POTTERS COVE
Aprime_pc = Matrix(Adata_pc[1:88,2:89])
#Only keep Role >1 (do not have length estimates for Role = 1 ~ primary producers)
keep = findall(x->x>1,vec(massdata_pc[!,:role]));
A_pc = copy(Aprime_pc[keep,keep]);
massvec_pc = massdata_pc[!,:length_cm][keep];
#delete completely disconnected species
keep2 = findall(!iszero,vec(sum(A_pc,dims=1)).+vec(sum(A_pc,dims=2)));
A_pc = A_pc[keep2,keep2];
massvec_pc = massvec_pc[keep2];
#Sorting
sortsp = sortperm(massvec_pc);
Asort_pc = A_pc[sortsp,sortsp];


Adata_bts = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BarentsSea/Barents_adjmatrix.csv",header=true);
massdata_bts = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BarentsSea/Barents_masslist.csv",header=true);
## BARENTS SEA
# NOTE: NEED TO REMOVE DETRITUS, ETC
Aprime_bts = Matrix(Adata_bts[1:159,2:160]);
role_bts = massdata_bts[!,:role];
#Exclude primary producers (role = 1), baleen whales (role = 4) and birds (role=3)
keep = findall(x->x>1 && x<3,role_bts);
A_bts = Aprime_bts[keep,keep];
massvec_bts = massdata_bts[!,:length_cm][keep];
namesvec_bts = Adata_bts[!,:species][keep];
keep2 = findall(x->x>1,massvec_bts);
A_bts = A_bts[keep2,keep2];
massvec_bts = massvec_bts[keep2];
namesvec_bts = namesvec_bts[keep2];
#Sorting
sortsp = sortperm(massvec_bts);
Asort_bts = A_bts[sortsp,sortsp];

#######################
# FITTING
#######################



## BENGUELA
x0 = [0.0,0.0,0.0];
results_bg = optimize(x->lfunc(x,A_bg,massvec_bg),x0,NelderMead());
results_bg.minimizer
xmax_bg = results_bg.minimizer;
fcorr_bg, Apredict_bg = fc(xmax_bg,A_bg,massvec_bg)
## POTTERS COVE
x0 = [0.0,0.0,0.0];
results_pc = optimize(x->lfunc(x,A_pc,massvec_pc),x0,NelderMead());
results_pc.minimizer
xmax_pc = results_pc.minimizer;
fcorr_pc, Apredict_pc = fc(xmax_pc,A_pc,massvec_pc)
## BARENTS SEA
x0 = [0.0,0.0,0.0];
results_bts = optimize(x->lfunc(x,A_bts,massvec_bts),x0,NelderMead());
results_bts.minimizer
xmax_bts = results_bts.minimizer;
fcorr_bts, Apredict_bts = fc(xmax_bts,A_bts,massvec_bts);

#Sort and isolate consumers only from Apredict
Apredict_bg_sort = Apredict_bg[sortsp_bg,sortsp_bg];
cons = findall(x->x>0,vec(sum(Apredict_bg_sort,dims=1)));
Apredict_bg_sort_cons = Apredict_bg_sort[cons,cons];

# Apredict_bg_sort = Apredict_bg[sortsp_bg,sortsp_bg];
cons = findall(x->x>0,vec(sum(Apredict_pc,dims=1)));
Apredict_pc_sort_cons = Apredict_pc[cons,cons];

# Apredict_bg_sort = Apredict_bg[sortsp_bg,sortsp_bg];
cons = findall(x->x>0,vec(sum(Apredict_bts,dims=1)));
Apredict_bts_sort_cons = Apredict_bts[cons,cons];



# use xmax to build a random food web

nw_bs = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/bodysizedata.csv",header=true);
nw_fmatrix = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/NWAtlantic/foraging_matrix.csv",header=true);
nw_minmass = nw_bs[!,Symbol("Mass min")];
nw_maxmass = nw_bs[!,Symbol("Mass max")];
nw_mass = vec(mean([nw_minmass nw_maxmass],dims=2));
trophic = nw_bs[!,Symbol("trophic")];
date = "1700";
occur = findall(!iszero,nw_bs[!,Symbol(date)]);
id1700 = collect(1:length(nw_mass))[occur];
bgmatrix_1700 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);

date = "2000";
occur = findall(!iszero,nw_bs[!,Symbol(date)]);
id2000 = collect(1:length(nw_mass))[occur];
bgmatrix_2000 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);


preymasses = [10^i for i=collect(1:0.1:4)];
predmasses = [10^i for i=collect(1:0.1:4)];
lpy = length(preymasses);
lpd = length(predmasses);
pijarray = Array{Float64}(undef,lpd,lpy);
alpha = xmax_bg[1];
beta = xmax_bg[2];
gamma = xmax_bg[3];
for i=1:lpd
    mi = predmasses[i];
    for j=1:lpy
        mj = preymasses[j];
        pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));
        pijarray[i,j] = pij/(1+pij);
    end
end

# Calculate some network metrics
R"""
library(UNODF)
library(igraph)
library(NetIndices)
library(fields)
library(plot.matrix)
library(RColorBrewer)
"""
reps = 5000;
rich = Array{Float64}(undef,reps,2);
conn = Array{Float64}(undef,reps,2);
nest = Array{Float64}(undef,reps,2);
mod = Array{Float64}(undef,reps,2);
nind2 = Array{Float64}(undef,reps,2);
nind3 = Array{Float64}(undef,reps,2);
maxtl = Array{Float64}(undef,reps,2);


nind2_1700_sp = Array{Float64}(undef,reps,23);
nind3_1700_sp = Array{Float64}(undef,reps,23);

nind2_2000_sp = Array{Float64}(undef,reps,17);
nind3_2000_sp = Array{Float64}(undef,reps,17);


for r=1:reps
    date = "1700";
    occur = findall(!iszero,nw_bs[!,Symbol(date)]);
    id1700 = collect(1:length(nw_mass))[occur];
    #This is a directed web
    bgmatrix_1700 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);
        # #What happens if web is undirected?
        # R"""
        # gdir = graph_from_adjacency_matrix($bgmatrix_1700);
        # gundir = as.undirected(gdir);
        # bgmatrix_1700 = as.matrix(as_adjacency_matrix(gundir));
        # """
        # @rget bgmatrix_1700;
    rich_1700 = size(bgmatrix_1700)[1];
    R"""
    unodfvalue = suppressWarnings(unodf($bgmatrix_1700,selfloop=TRUE))
    nestvalue_c = unodfvalue$UNODFc
    nestvalue_r = unodfvalue$UNODFr
    """
    @rget nestvalue_c;
    @rget nestvalue_r;
    nest_1700 = mean([nestvalue_c,nestvalue_r]);

    R"""
    g_ud <- as.undirected(graph.adjacency($bgmatrix_1700));
    wtc <- cluster_walktrap(g_ud);
    mod_1700 <- modularity(g_ud,membership(wtc));
    """
    @rget mod_1700;

    #number of indirect paths
    # 3
    aprime2_1700 = bgmatrix_1700^2;
    nind2_1700 = sum(aprime2_1700);
    aprime3_1700 = bgmatrix_1700^3;
    nind3_1700 = sum(aprime3_1700);

    #TL
    R"ctl_1700 <- TrophInd($(bgmatrix_1700))[,1]"
    @rget ctl_1700;
    

    
    rich[r,1] = rich_1700;
    conn[r,1] = sum(bgmatrix_1700)/size(bgmatrix_1700)[1]^2;
    nest[r,1] = nest_1700;
    mod[r,1] = mod_1700;
    nind2[r,1] = nind2_1700;
    nind3[r,1] = nind3_1700;
    maxtl[r,1] = maximum(ctl_1700);

    date = "2000";
    occur = findall(!iszero,nw_bs[!,Symbol(date)]);
    id2000 = collect(1:length(nw_mass))[occur];
    bgmatrix_2000 = bsfoodweb(xmax_bg,nw_mass,trophic,occur);
        # #What happens if web is undirected?
        # R"""
        # gdir = graph_from_adjacency_matrix($bgmatrix_2000);
        # gundir = as.undirected(gdir);
        # bgmatrix_2000 = as.matrix(as_adjacency_matrix(gundir));
        # """
        # @rget bgmatrix_2000;
    rich_2000 = size(bgmatrix_2000)[1];
    R"""
    unodfvalue = suppressWarnings(unodf($bgmatrix_2000,selfloop=TRUE))
    nestvalue_c = unodfvalue$UNODFc
    nestvalue_r = unodfvalue$UNODFr
    """
    @rget nestvalue_c;
    @rget nestvalue_r;
    nest_2000 = mean([nestvalue_c,nestvalue_r]);

    R"""
    g_ud <- as.undirected(graph.adjacency($bgmatrix_2000));
    wtc <- cluster_walktrap(g_ud);
    mod_2000 <- modularity(g_ud,membership(wtc));
    """
    @rget mod_2000;

    #number of indirect paths
    # 3
    aprime2_2000 = bgmatrix_2000^2;
    aprime3_2000 = bgmatrix_2000^3;
    nind2_2000 = sum(aprime2_2000);
    nind3_2000 = sum(aprime3_2000);

    #TL
    R"ctl_2000 <- TrophInd($(bgmatrix_2000))[,1]"
    @rget ctl_2000;
    

    rich[r,2] = rich_2000;
    conn[r,2] = sum(bgmatrix_2000)/size(bgmatrix_2000)[1]^2;
    nest[r,2] = nest_2000;
    mod[r,2] = mod_2000;
    nind2[r,2] = nind2_2000;
    nind3[r,2] = nind3_2000;
    maxtl[r,2] = maximum(ctl_2000);

    # Contribution of each species
    #1700
    splist = collect(1:rich_1700);
    for i=1:rich_1700
        spadj = setdiff(splist,i);
        a_adj = bgmatrix_1700[spadj,spadj];
        aprime2_adj = a_adj^2;
        nind2_1700_sp[r,i] = nind2_1700 / sum(aprime2_adj);
        aprime3_adj = a_adj^3;
        nind3_1700_sp[r,i] = nind3_1700 / sum(aprime3_adj);
    end
    splist = collect(1:rich_2000);
    for i=1:rich_2000
        spadj = setdiff(splist,i);
        a_adj = bgmatrix_1700[spadj,spadj];
        aprime2_adj = a_adj^2;
        nind2_2000_sp[r,i] = nind2_2000 / sum(aprime2_adj);
        aprime3_adj = a_adj^3;
        nind3_2000_sp[r,i] = nind3_2000 / sum(aprime3_adj);
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
library(fields)
library(plot.matrix)
library(RColorBrewer)
pdf($namespace,width=14,height=12)
par(mfrow=c(2,2))
# layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
# par(list(oma = c(3, 3, 0.5, 0.5), mar = c(3, 3, 3, 1)))


pal = rev(brewer.pal(9,'YlOrRd'))
pijarray = $pijarray
rownames(pijarray) = $predmasses
colnames(pijarray) = $preymasses
# plot((pijarray),main='',border=NA,col=pal,xlab='Pred mass (kg)',ylab='Prey mass (kg)')
image.plot(x=$predmasses,y=$preymasses,z=$pijarray,col=pal,legend.args=list(text='     Pr(link)',
side=3,
line=.5), xlab='Pred mass (kg)', ylab = 'Prey mass (kg)',cex.axis=1.5, cex.lab=1.5
)

pal = rev(brewer.pal(9,'YlGnBu'))
plot($(Apredict_bg_sort_cons),axis.col=3,main='',col=pal,fmt.key="%.2f",xlab='',ylab='',cex.axis=1.5, cex.lab=1.5)

# plot($(bgmatrix_1700),axis.col=3,main='',key=NULL)
pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max(c($ctlcolors_1700,$ctlcolors_2000)))
agraph = graph_from_adjacency_matrix($(bgmatrix_1700))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_1700
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_1700],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white',vertex.label=$id1700)
# legend(0,10,lv,col=pal[lv],pch=16)

# plot($(bgmatrix_2000),axis.col=3,main='',key=NULL)
# pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max($ctlcolors_2000))
agraph = graph_from_adjacency_matrix($(bgmatrix_2000))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_2000
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_2000],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white',vertex.label=$id2000)
dev.off()
"""



#inset figure
namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panelinset.pdf";
R"""
pdf($namespace,width=5,height=4)
# par(mfrow=c(2,2))
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(list(oma = c(3, 2, 0.5, 0.5), mar = c(1, 3, 0, 1)))
boxplot($rich,names=c('',''),boxwex=0.25)
boxplot($conn,names=c('',''),boxwex=0.25)
boxplot($nest,names=c('1700','2000'),boxwex=0.25)
boxplot($mod,names=c('1700','2000'),boxwex=0.25)
dev.off()
"""


#inset figure
namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panelinset_trim.pdf";
R"""
pdf($namespace,width=5,height=2)
# par(mfrow=c(2,2))
layout(matrix(c(1,2),1,2,byrow=TRUE), widths=c(1,1), heights=c(0.5,0.5))
par(list(oma = c(1, 0.25, 0.5, 0.25), mar = c(1, 2, 0, 1)))
boxplot($rich,names=c('1700','2000'),boxwex=0.25)
boxplot($conn,names=c('1700','2000'),boxwex=0.25)
dev.off()
"""

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panelinset_trim_indirect.pdf";
R"""
pdf($namespace,width=2,height=3.5)
# par(mfrow=c(2,2))
pal = brewer.pal(3,'Set1')
layout(matrix(c(1,2),2,1,byrow=TRUE), widths=c(1,1), heights=c(0.5,0.5))
par(list(oma = c(1, 0.25, 0.5, 0.25), mar = c(1, 2, 0, 1)))
boxplot($nind2,names=c('',''),boxwex=0.25,col=pal[2])
boxplot($nind3,names=c('1700','2000'),boxwex=0.25,col=pal[3])
dev.off()
"""


namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panelinset_trim_indirectnorm.pdf";
R"""
pdf($namespace,width=2,height=3.5)
# par(mfrow=c(2,2))
pal = brewer.pal(3,'Set1')
layout(matrix(c(1,2),2,1,byrow=TRUE), widths=c(1,1), heights=c(0.5,0.5))
par(list(oma = c(1, 0.25, 0.5, 0.25), mar = c(1, 2, 0, 1)))
boxplot($nind2 / $rich,names=c('',''),boxwex=0.25,col=pal[2])
boxplot($nind3 / $rich,names=c('1700','2000'),boxwex=0.25,col=pal[3])
dev.off()
"""



#legend
namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panel_legend.pdf";
R"""
pdf($namespace,width=5,height=6)
plot(seq(1,20),seq(1,20))
lv = seq(0,max(c($ctlcolors_1700,$ctlcolors_2000)),length.out=4); lv[1] = 1;
lvnum = c(1,2,3,4)
legend(15,12,rev(lvnum),col=rev(pal[lv]),pch=16,cex=1.5)
dev.off()
"""



# indirect paths species by species
date = "1700";
occur1700 = findall(!iszero,nw_bs[!,Symbol(date)]);

date = "2000";
occur2000 = findall(!iszero,nw_bs[!,Symbol(date)]);

occurboth = intersect(occur1700,occur2000)
namesboth = nw_bs[!,:Species][occurboth];

alist = [findall(x->x==occurboth[i],occur1700) for i=1:17];
apos = collect(Iterators.flatten(alist));

# scatterplot(mspindirect2_1700[apos],mspindirect2_2000)
# scatterplot(mspindirect3_1700[apos],mspindirect3_2000)

mspindirect2_1700 = vec(mean(nind2_1700_sp,dims=1));
mspindirect2_2000 = vec(mean(nind2_2000_sp,dims=1));
mspindirect3_1700 = vec(mean(nind3_1700_sp,dims=1));
mspindirect3_2000 = vec(mean(nind3_2000_sp,dims=1));

fullind2_1700 = nind2_1700_sp[:,apos];
fullind3_1700 = nind3_1700_sp[:,apos];
fullind2_2000 = nind2_2000_sp;
fullind3_2000 = nind3_2000_sp;

namesboth[findall(x->x>1.5,mspindirect3_2000)]
apos[findall(x->x>1.5,mspindirect3_2000)]


namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_4panelinset_trim_indirect_species.pdf";
R"""
pdf($namespace,width=3.5,height=3.5)
# par(mfrow=c(2,2))
layout(matrix(c(1),1,1,byrow=TRUE), widths=c(1), heights=c(0.5))
par(list(oma = c(3, 2, 0.5, 0.25), mar = c(1, 2, 0, 1)))
pal = brewer.pal(5,'Set1')
plot($(mspindirect2_1700[apos]),$(mspindirect2_2000),pch=21,bg=pal[2],col='black',xlim=c(1,1.8),ylim=c(1,1.8),cex=1.25,xlab='',ylab='')
# points($(vec(fullind2_1700)),$(vec(fullind2_2000)),pch='.',col=paste(pal[2],'10',sep=''))
# points($(vec(fullind3_1700)),$(vec(fullind3_2000)),pch='.',col=paste(pal[3],'10',sep=''))
lines(seq(0.5,2.5),seq(0.5,2.5),col='gray')
points($(mspindirect2_1700[apos]),$(mspindirect2_2000),pch=21,bg=pal[2],col='black',cex=1.25)
points($(mspindirect3_1700[apos]),$(mspindirect3_2000),pch=21,bg=pal[3],col='black',cex=1.25)
points($(mspindirect2_1700[apos]),$(mspindirect2_2000),pch=16,col=pal[2],cex=1.25)
points($(mspindirect3_1700[apos]),$(mspindirect3_2000),pch=16,col=pal[3],cex=1.25)
# text($(mspindirect3_1700[apos]),$(mspindirect3_2000),$namesboth,cex=0.5)
# lines(seq(0,2,length.out=10),rep(1,10),col=pal[4])
# lines(rep(1,10),seq(0,2,length.out=10),col=pal[4])
mtext(side=1,expression(paste("Ratio indirect interactions (1700)")),line=2.5)
mtext(side=2,expression(paste("Ratio indirect interactions (2000)")),line=2)
dev.off()
"""





# layoutmatrix = reshape([1,1,3,3,6,6,1,1,3,3,6,6,1,1,4,4,7,7,2,2,4,4,7,7,2,2,5,5,8,9,2,2,5,5,8,9],(6,6))
namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/foodweb_9panel.pdf";
R"""
library(igraph)
library(fields)
library(plot.matrix)
library(RColorBrewer)
pdf($namespace,width=20,height=20)
# par(mfrow=c(2,2))
# layout(($layoutmatrix), widths=rep(1,6*6), heights=rep(1,6*6))
# par(list(oma = c(3, 3, 3, 3), mar = c(3, 3, 3, 3)))
par(list(oma = c(3, 3, 0, 0), mar = c(2, 2, 1, 0)))


pal = rev(brewer.pal(9,'YlOrRd'))
pijarray = $pijarray
rownames(pijarray) = $predmasses
colnames(pijarray) = $preymasses
margin = 0.05
# PLOT 1
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0*(1+0.01), 0.33*(1-0.01), .66*(1+0.01), 1*(1-0.01))))
image(t(apply(t($Asort_bts), 1, rev)),col=c('white','black'),legend=FALSE,xaxt='n',yaxt='n')
# axis(1,labels=FALSE,tick=FALSE)
# axis(2,labels=FALSE,tick=FALSE)
# axis(3,labels=FALSE,tick=FALSE)
# axis(4,labels=FALSE,tick=FALSE)
# axis.col=NULL, axis.row=NULL,
# PLOT 2
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(.33*(1+0.01), 0.66*(1-0.01), .66*(1+0.01), 1*(1-0.01))))
# plot($Asort_pc,col=c('white','black'),key=NULL, axis.col=NULL, axis.row=NULL,border=NA,main='',xlab='',ylab='')
image(t(apply(t($Asort_pc), 1, rev)),col=c('white','black'),legend=FALSE,xaxt='n',yaxt='n')

# PLOT 3
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(.66*(1+0.0), 1*(1-0.01), .66*(1+0.01), 1*(1-0.01))))
# plot($Asort_bg,col=c('white','black'),key=NULL, axis.col=NULL, axis.row=NULL,border=NA,main='',xlab='',ylab='')
image(t(apply(t($Asort_bg), 1, rev)),col=c('white','black'),legend=FALSE,xaxt='n',yaxt='n')

# PLOT 4
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0*(1+margin), 0.5*(1-margin), .33*(1+margin), 0.66*(1-0))))
image.plot(x=$predmasses,y=$preymasses,z=$pijarray,col=pal,legend.args=list(text='     Pr(link)',
side=3,
line=.5), xlab='Pred mass (kg)', ylab = 'Prey mass (kg)',cex.axis=1.5, cex.lab=1.5
)
pal = rev(brewer.pal(9,'YlGnBu'))

# PLOT 5
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0.55*(1+margin), 1*(1-margin), .33*(1+margin), 0.66*(1-0))))
plot($(Apredict_bg_sort_cons),axis.col=1,main='',col=pal,fmt.key="%.2f",xlab='',ylab='',cex.axis=1.5, cex.lab=1.5)

# Plot 6
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0*(1+margin), 0.33*(1-margin), 0*(1+margin), 0.33*(1-margin))))
# plot($(bgmatrix_1700),axis.col=3,main='',key=NULL)
pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max(c($ctlcolors_1700,$ctlcolors_2000)))
agraph = graph_from_adjacency_matrix($(bgmatrix_1700))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_1700
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_1700],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white',vertex.label=$id1700)
# legend(0,10,lv,col=pal[lv],pch=16)

# Plot 7
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(.33*(1+margin), 0.66*(1-margin), 0*(1+margin), 0.33*(1-margin))))
# plot($(bgmatrix_2000),axis.col=3,main='',key=NULL)
# pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))(max($ctlcolors_2000))
agraph = graph_from_adjacency_matrix($(bgmatrix_2000))
lay <- layout.fruchterman.reingold(agraph)
# lay <- layout_on_sphere(agraph)
lay[,2] <- ctl_2000
plot.igraph(agraph,layout = lay,vertex.color=pal[$ctlcolors_2000],vertex.size=18,edge.arrow.size=0.5,vertex.label.color='white',vertex.label=$id2000)

pal = brewer.pal(3,'Set1')
# PLOT 8
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0.66*(1+margin), 1*(1-margin), 0.16*(1+margin), 0.33*(1-margin))))
boxplot($nind2,names=c('',''),boxwex=0.25,col=pal[2])

# PLOT 9
# c(x1, x2, y1, y2)
par(list(new=TRUE, plt=c(0.66*(1+margin), 1*(1-margin), 0*(1+margin), 0.16*(1-margin))))
boxplot($nind3,names=c('1700','2000'),boxwex=0.25,col=pal[3])

dev.off()
"""
