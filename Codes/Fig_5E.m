%% Y212 (mCherry codon optimized) vs Y197 (non mCherry codon optimized)

%% load dta

clear all
cd('~/Google Drive File Stream/Mi unidad/CareyLab/Feedback Detector/RBP library/RBP Feedback detector Manuscript/New Work/Data_GitHub/')
load 'mCherry_codon_bias.mat'
DataSet_Old = DS;

%% Plot Y197 vs Y212

close all
clc

%define DS:
idx = find(DataSet_Old.NCellsAfterFiltering >= 300 & DataSet_Old.TotalEvents <= 70000);
DS = DataSet_Old(idx,:);
DS.StrainGood = cellfun(@(x) x(1:3),DS.Strain,'UniformOutput',false);

Orange = [1,0.4,0];
Gray = [.7,.7,.7];
Black = [0,0,0];
LightPurple = [168,112,219]./255;
DarkPurple = [55,33,75]./255;
Strains = {'Y19','COm'};
colors = [DarkPurple;LightPurple];

fh = figure('units','centimeters','position',[5 5 6.5 6.5]);
hold on; 
for S = 1:2
    idx212 = find(strcmp(DS.StrainGood,'COm'));
    Y212_GFP = log2(DS.mean_FITC_A(idx212));
    Y212_RFP = log2(DS.mean_PE_TexasRed_A(idx212));
    
    idxS = find(strcmp(DS.StrainGood,Strains(S)));
    GFP = log2(DS.mean_FITC_A(idxS));
    RFP = log2(DS.mean_PE_TexasRed_A(idxS));
    [~,o] = sort(GFP);
    GFP = GFP(o);
    RFP = RFP(o);
        
    [ GFPfit,RFPfit_noF,RFPfit_F ] = Fit_RfG_model( GFP, RFP, Y212_GFP , Y212_RFP );  
    
    plot(GFP,RFP,'.','color',colors(S,:),'markerfacecolor',colors(S,:),'linewidth',2,'markersize',12)
    plot(GFPfit,RFPfit_noF,'-','color',Gray,'linewidth',2)
    plot(GFPfit,RFPfit_F,'-','color',colors(S,:),'linewidth',2)
    
    FoldChangeInExpression = max(RFPfit_F) - max(RFPfit_noF)
    
    %xlim([7.9,11.2])
    %ylim([8,10.5])                
      
    xlabel('log_2 GFP')
    ylabel('log_2 mCherry')
end





