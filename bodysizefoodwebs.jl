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

#NOTE: 10/28/20 - Build in rectangular matrix; calculate %correct and %incorrect for both links and nonlinks

#Load food web interaction matrix
# Aload = ;
# massload = ;

Adata_bg = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BenguelaSea/bengula_A_matrix.csv",header=true);
numsp = size(Adata_bg)[1];
massdata_bg = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BenguelaSea/bengula_masses.csv",header=true);

Adata_pc = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/adjmatrix.csv",header=true);
massdata_pc = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/masslist.csv",header=true);

Adata_bts = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BarentsSea/Barents_adjmatrix.csv",header=true);
massdata_bts = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/BarentsSea/Barents_masslist.csv",header=true);

## BENGUELA
#Use species 5-29 (species 1-4 have missing diet data)
Aprime_bg = Matrix(Adata_bg[1:numsp,2:numsp+1]);
# nolinks = findall(iszero,sum(Aprime_bg,dims=1));
# A_bg = Matrix(Adata[5:29,6:30]);
A_bg = copy(Aprime_bg);
# massvec_bg = massdata[:mass][5:numsp];
massvec_bg = copy(vec(massdata_bg[!,:mass_kg]));

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

## BARENTS SEA
Aprime_bts = Matrix(Adata_bts[1:159,2:160]);
role_bts = massdata_bts[!,:role];
#Exclude primary producers (role = 1), baleen whales (role = 4) and birds (role=3)
keep = findall(x->x>1 && x<3,role_bts);
A_bts = Aprime_bts[keep,keep];
massvec_bts = massdata_bts[!,:length_cm][keep];
keep2 = findall(x->x<5,massvec_bts);
A_bts = A_bts[keep2,keep2];
massvec_bts = massvec_bts[keep2];


#######################
# FITTING
#######################



## BENGUELA
x0 = [0.0,0.0,0.0];
results_bg = optimize(x->lfunc(x,A_bg,massvec_bg),x0,NelderMead());
results_bg.minimizer
xmax = results_bg.minimizer;
fcorr_bg, Apredict_bg = fc(xmax,A_bg,massvec_bg)

x0 = [0.0,0.0,0.0];
Anull_bg = nullweb(A_bg);
resultsnull_bg = optimize(x->lfunc(x,Anull_bg,massvec_bg),x0,NelderMead());
resultsnull_bg.minimizer
xmax = results_bg.minimizer;
fcorrnull_bg, Apredictnull_bg = fc(xmax,Anull_bg,massvec_bg)

## POTTERS COVE
x0 = [0.0,0.0,0.0];
results_pc = optimize(x->lfunc(x,A_pc,massvec_pc),x0,NelderMead());
results_pc.minimizer
xmax = results_pc.minimizer;
fc(xmax,A_pc,massvec_pc)

## BARENTS SEA
x0 = [0.0,0.0,0.0];
results_bts = optimize(x->lfunc(x,A_bts,massvec_bts),x0,NelderMead());
results_bts.minimizer
xmax = results_bts.minimizer;
fcorr_bts, Apredict_bts = fc(xmax,A_bts,massvec_bts);

Anull_bts = nullweb(A_bts);
resultsnull_bts = optimize(x->lfunc(x,Anull_bts,massvec_bts),x0,NelderMead());
resultsnull_bts.minimizer
#NOTE: these maxL params result in flat probabilities... so guessing the problem is in the fc function
xmax = resultsnull_bts.minimizer;
# xmax = results_bts.minimizer;
fcorrnull_bts, Apredictnull_bts = fc(xmax,Anull_bts,massvec_bts);


# x0 = [0.0,0.0,0.0];
# results_pc = optimize(x->lfunc(x,A_pc,massvec_pc),x0,SimulatedAnnealing(),Optim.Options(iterations=50000));
# results_pc.minimizer


# sortsp = sortperm(vec(sum(A_pc,dims=1)));
sortsp = sortperm(massvec_pc);
Asort_pc = A_pc[sortsp,sortsp];
# sortsp = sortperm(vec(sum(A_bg,dims=1)));
sortsp = sortperm(massvec_bg);
Asort_bg = A_bg[sortsp,sortsp];

sortsp = sortperm(massvec_bts);
Asort_bts = A_bts[sortsp,sortsp];
Anull_bts = nullweb(Asort_bts);

p1=heatmap(Array{Int64}(Asort_bg),zlims=(0,1));
p2=heatmap(Array{Int64}(Asort_pc),zlims=(0,1));
p3=heatmap(Array{Int64}(Asort_bts),zlims=(0,1));
p4=heatmap(Array{Int64}(Anull_bts),zlims=(0,1));
plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/fig.pdf")


sortsp = sortperm(massvec_bts);
Asort_bts = A_bts[sortsp,sortsp];
Anull_bts = Anull_bts[sortsp,sortsp];
p1=heatmap(Array{Int64}(Asort_bts),zlims=(0,1));
p2=heatmap(Array{Int64}(Anull_bts),zlims=(0,1));
p3=heatmap(Array{Float64}(Apredict_bts[sortsp,sortsp]),zlims=(0,1));
p4=heatmap(Array{Float64}(Apredictnull_bts[sortsp,sortsp]),zlims=(0,1));
plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/Anull_bts.pdf")


sortsp = sortperm(massvec_bg);
Asort_bg = A_bg[sortsp,sortsp];
Anull_bg = Anull_bg[sortsp,sortsp];
p1=heatmap(Array{Int64}(Asort_bg),zlims=(0,1));
p2=heatmap(Array{Int64}(Anull_bg),zlims=(0,1));
p3=heatmap(Array{Float64}(Apredict_bg[sortsp,sortsp]),zlims=(0,1));
p4=heatmap(Array{Float64}(Apredictnull_bg[sortsp,sortsp]),zlims=(0,1));
plot(p1,p2,p3,p4,layout=(2,2),size=(1500,1000))
savefig("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/Anull_bg.pdf")
