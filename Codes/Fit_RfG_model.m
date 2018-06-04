function [ GFPfit,RFPfit_noF,RFPfit_F,r2_basal,r2_F,SD,BE,F] = Fit_RfG_model(GFP,RFP,Y212_GFP,Y212_RFP,n_tries)
%This function takes GFP, RFP, GFP_Y212 and RFP_Y212 and outputs GFPfit,
%RFPfit_F, RFPfit_noF, FeedbackStrength, r2 of the basal model, r2 of the F
%model, Synthesis/Degradation, mCherry Basal Expression and Ftl

Cg_range = [0,1e6];
BE_r_range = [0,1];
BE_g_range = [0,1];
Tc_r_range = [0,1e7];
Tl_r_range = [0,1e7];
D_r_range = [0,1e6];
Ftl_range = [-100,10000]; 

Y212_lb  = [min(Cg_range),min(BE_r_range),min(BE_g_range),min(Tc_r_range),min(Tl_r_range),min(D_r_range),0];
Y212_ub  = [max(Cg_range),max(BE_r_range),max(BE_g_range),max(Tc_r_range),max(Tl_r_range),max(D_r_range),0];

Best_Cg = 0;
Best_BE_r = 0;
Best_BE_g = 0;
Best_Tc_r = 0;
Best_Tl_r = 0;
Best_D_r = 0;
Best_Ftl = 0;

GFPdata = linspace(min(Y212_GFP),max(Y212_GFP),10000)';

Nfits = 0;
Y212_R2 = 0; 
while Nfits <= n_tries & Y212_R2 < 0.999       
    b0_Y212 = random('uniform',Y212_lb,Y212_ub); 
    try    
        RFPeq = @(b,Y212_GFP) log2(((((((2.^(Y212_GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76) + (1 - b(3).*(1 + (2.^(Y212_GFP - log2(b(1)*b(3))) - 1))).*b(2))./((1 - b(3).*(1 + (2.^(Y212_GFP - log2(b(1)*b(3))) - 1))) + ((2.^(Y212_GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76))) ./ (1 - ((((2.^(Y212_GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76) + (1 - b(3).*(1 + (2.^(Y212_GFP - log2(b(1)*b(3))) - 1))).*b(2))./((1 - b(3).*(1 + (2.^(Y212_GFP - log2(b(1)*b(3))) - 1))) + ((2.^(Y212_GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76))).*((b(7)*b(4))./b(6)))).*(b(4)*b(5))) ./ b(6));
        [b_RFP, lsq_resnorm,lsq_residuals, lsq_exitflag, lsq_output] = lsqcurvefit(RFPeq, b0_Y212,Y212_GFP,Y212_RFP, Y212_lb, Y212_ub);  %b_RFP cointains the desired 7 parameters  
        Cg = real(b_RFP(1));
        BE_r = real(b_RFP(2));
        BE_g = real(b_RFP(3));
        Tc_r = real(b_RFP(4));
        Tl_r = real(b_RFP(5));
        D_r = real(b_RFP(6));
        Ftl = real(b_RFP(7));  
        RFP_fit_R2_Y212 = log2(((((((2.^(Y212_GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(Y212_GFP - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(Y212_GFP - log2(Cg*BE_g)) - 1))) + ((2.^(Y212_GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76))) ./ (1 - ((((2.^(Y212_GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(Y212_GFP - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(Y212_GFP - log2(Cg*BE_g)) - 1))) + ((2.^(Y212_GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76))).*((Ftl*Tc_r)./D_r))).*(Tc_r*Tl_r)) ./ D_r);       
        [r2,rmse] = rsquare(Y212_RFP,RFP_fit_R2_Y212);
        if r2 > Y212_R2 
            Y212_R2 = r2;          
            Cg_212 = Cg;
            BE_r_212 = BE_r;
            BE_g_212 = BE_g;
            Tc_r_212 = 	Tc_r;
            Tl_r_212 = Tl_r;
            D_r_212 = D_r;
            Ftl_212 = Ftl;
        end 
        Nfits = Nfits + 1;
     end
end

GFPdata = linspace(min(GFP),max(GFP),100000)'; 
for J = 1:2 
    Cg = Cg_212;
    BE_r = BE_r_212;
    BE_g = BE_g_212;
    Tc_r = Tc_r_212;
    Tl_r = Tl_r_212;
    D_r = D_r_212;
    Ftl = 0;    
    if J == 1; BE_r = min(BE_r_range); Tl_r = min(Tl_r_range);end  
    if J == 2; BE_r = min(BE_r_range); Tl_r = min(Tl_r_range);Ftl = min(Ftl_range); end                        
    lb = [Cg,BE_r,BE_g,Tc_r,Tl_r,D_r,Ftl];        
    if J == 1; BE_r = max(BE_r_range); Tl_r = max(Tl_r_range);end    
    if J == 2; BE_r = max(BE_r_range); Tl_r = max(Tl_r_range);Ftl = max(Ftl_range); end                    
    ub = [Cg,BE_r,BE_g,Tc_r,Tl_r,D_r,Ftl];    
    Nfits = 0;
    TriedFits = 0;
    if J == 1
        R2 = 0; 
    end
    if J == 2
        R2 = R2withoutF; 
    end                
    while Nfits <= n_tries & R2 < 0.999 & TriedFits <= 100 
        if J == 1 
            b0 = random('uniform',lb,ub); 
        end
        if J == 2  
            b0 = Best_parms_without_F;
        end
        try
            RFPeq = @(b,GFP) log2(((((((2.^(GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76) + (1 - b(3).*(1 + (2.^(GFP - log2(b(1)*b(3))) - 1))).*b(2))./((1 - b(3).*(1 + (2.^(GFP - log2(b(1)*b(3))) - 1))) + ((2.^(GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76))) ./ (1 - ((((2.^(GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76) + (1 - b(3).*(1 + (2.^(GFP - log2(b(1)*b(3))) - 1))).*b(2))./((1 - b(3).*(1 + (2.^(GFP - log2(b(1)*b(3))) - 1))) + ((2.^(GFP - log2(b(1)*b(3))) - 1).*b(3).*0.76))).*((b(7)*b(4))./b(6)))).*(b(4)*b(5))) ./ b(6));
            [b_RFP, lsq_resnorm,lsq_residuals, lsq_exitflag, lsq_output] = lsqcurvefit(RFPeq, b0,GFP,RFP, lb, ub);
            Cg = real(b_RFP(1));
            BE_r = real(b_RFP(2));
            BE_g = real(b_RFP(3));
            Tc_r = real(b_RFP(4));
            Tl_r = real(b_RFP(5));
            D_r = real(b_RFP(6));
            Ftl = real(b_RFP(7));  
            RFP_fit_R2 = log2(((((((2.^(GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFP - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFP - log2(Cg*BE_g)) - 1))) + ((2.^(GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76))) ./ (1 - ((((2.^(GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFP - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFP - log2(Cg*BE_g)) - 1))) + ((2.^(GFP - log2(Cg*BE_g)) - 1).*BE_g.*0.76))).*((Ftl*Tc_r)./D_r))).*(Tc_r*Tl_r)) ./ D_r);
            [r2,rmse] = rsquare(RFP,RFP_fit_R2);
            if r2 > R2                 
                R2 = r2; 
                BestFtl = Ftl;  
                Best_Cg = Cg;
                Best_BE_r = BE_r;
                Best_BE_g = BE_g;
                Best_Tc_r = Tc_r;
                Best_Tl_r = Tl_r;
                Best_D_r = D_r;
                Best_Ftl = Ftl; 
            end 
            Nfits = Nfits + 1;
        end
        TriedFits = TriedFits + 1;
    end
    Cg = Best_Cg;
    BE_r = Best_BE_r;
    BE_g = Best_BE_g;
    Tc_r = Best_Tc_r;
    Tl_r = Best_Tl_r;
    D_r = Best_D_r;
    Ftl = Best_Ftl;         
    if J == 1 
        Best_parms_without_F = [Cg,BE_r,BE_g,Tc_r,Tl_r,D_r,0];
        R2withoutF = R2;
        RFP_fit_without_F_ranging = log2(((((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))) ./ (1 - ((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))).*((Ftl*Tc_r)./D_r))).*(Tc_r*Tl_r)) ./ D_r);            
    end
    if J == 2  
        R2withF = R2;
        RFP_fit_with_F = log2(((((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))) ./ (1 - ((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))).*((Ftl*Tc_r)./D_r))).*(Tc_r*Tl_r)) ./ D_r);                   
        Ftl = 0;
        RFP_fit_without_F = log2(((((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))) ./ (1 - ((((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76) + (1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))).*BE_r)./((1 - BE_g.*(1 + (2.^(GFPdata - log2(Cg*BE_g)) - 1))) + ((2.^(GFPdata - log2(Cg*BE_g)) - 1).*BE_g.*0.76))).*((Ftl*Tc_r)./D_r))).*(Tc_r*Tl_r)) ./ D_r);            
        FoldChangeInExpression = max(RFP_fit_with_F) - max(RFP_fit_without_F);
    end
end

GFPfit = GFPdata;
RFPfit_noF = RFP_fit_without_F;
RFPfit_F = RFP_fit_with_F;
r2_basal = R2withoutF;
r2_F = R2withF;
SD = Tl_r;
BE = BE_r;
F = Best_Ftl;

end
