%%
close all;

finalFnFitValLog = dlmread('finalFnFitValLog.csv');

dValueMatrix = dlmread('sep_id_dCntr_dMN_dVal_ttaRad_ttaDeg.csv');


%%
dvalues = dValueMatrix(:,7);
logDval = log10(dvalues);

RijUnique = unique(finalFnFitValLog(:,1));

oddLogDval = logDval(1:2:64);
evenLogDval = logDval(2:2:64);

% figure,plot(RijUnique,evenLogDval);
% hold on;
% plot(RijUnique,oddLogDval, 'r');
% hold off;

%%
Di = finalFnFitValLog(:,4);
%Si = 1./(1+exp(-b*(Di-c)));
Si = heaviside(Di+5);

Dj = finalFnFitValLog(:,6);
%Sj = 1./(1+exp(-b*(Dj-c)));
Sj = heaviside(Dj+5);

Ttai = finalFnFitValLog(:,5);
Ttaj = finalFnFitValLog(:,7);

Cos_Ttai_Ttaj = cos(Ttai-Ttaj);

finalFnVal = finalFnFitValLog(:,1:3);

finalFnVal(:,1) = (finalFnVal(:,1)+110)*10^-7;
finalFnVal(:,4) = Si;
finalFnVal(:,5) = Sj;
finalFnVal(:,6) = Cos_Ttai_Ttaj;
finalFnVal(:,7) = finalFnFitValLog(:,8);

%%

SiSjCos = Si.*Sj.*Cos_Ttai_Ttaj;

Rij = finalFnVal(:,1);

EijAll = finalFnVal(:,7);

%%

%figure, scatter3(Rij, SiSjCos, EijAll, 'b*');

%figure, plot3(Rij ,SiSjCos ,EijAll, 'b*');

%%
DiDjCosIdx = find(SiSjCos > 9.989559367286183e-01);
f_DiDjCos = SiSjCos(DiDjCosIdx);
f_Eij = EijAll(DiDjCosIdx);
f_Rij = Rij(DiDjCosIdx);

%%
rij_didjCos_Eij = [Rij SiSjCos EijAll];
%dlmwrite('rij_didjCos_Eij.csv',rij_didjCos_Eij);

%%
LLG_sep_indvE_interE = dlmread('LLG_sep_indvE_interE.csv');

LLG_sep_indvE_interE(:,2) = LLG_sep_indvE_interE(:,3)*1e-7;

%%
gamma_a = -3.735e-06;
rho_b= -2.452e+05;

aeRij = gamma_a*exp(rho_b*Rij);
aeRij_DiDjCos = aeRij.*SiSjCos;
aeRij_Cos = aeRij.*Cos_Ttai_Ttaj;

%%
beta_a = 2.5*1e-7;
omega_d = 6.4*1e-8;

ESi = (beta_a.*Si + omega_d);
ESj = (beta_a.*Sj + omega_d);
ETot = aeRij_DiDjCos + ESi + ESj;

%%
% sep_sepIncm_id1_id2_Si_Sj_CosTi_Tj_ESi_ESj_Ec_Et
superFinalVal =  finalFnFitValLog(:,1);
superFinalVal(:,2:7) =  finalFnVal(:,1:6);
superFinalVal(:,8) = ESi;
superFinalVal(:,9) = ESj;
superFinalVal(:,10) = aeRij_DiDjCos;
superFinalVal(:,11) = ETot;


%%
sep_fittedE_llgE_fittedTotalE = LLG_sep_indvE_interE(:,1);

dot1Id = 1;
dot2Id = 2;

for i=1:size(sep_fittedE_llgE_fittedTotalE)
    
    sep = sep_fittedE_llgE_fittedTotalE(i,1);
    
    fittedE_idx = find(finalFnFitValLog(:,1) == sep & finalFnFitValLog(:,2) == dot1Id & finalFnFitValLog(:,3) == dot2Id);
    llgE_idx = find(LLG_sep_indvE_interE(:,1)==sep);
    
    sep_fittedE_llgE_fittedTotalE(i,2) = aeRij_DiDjCos(fittedE_idx,1);
    sep_fittedE_llgE_fittedTotalE(i,3) = LLG_sep_indvE_interE(llgE_idx, 3);
    sep_fittedE_llgE_fittedTotalE(i,4) = superFinalVal(fittedE_idx, 11);
    sep_fittedE_llgE_fittedTotalE(i,5) = Rij(fittedE_idx, 1);
    sep_fittedE_llgE_fittedTotalE(i,6) = aeRij_Cos(fittedE_idx, 1);
    
    dot1Id = dot1Id + 2;
    dot2Id = dot2Id + 2;
    
end

%%
%% PLOTS
fittedE = sep_fittedE_llgE_fittedTotalE(2:end,2);
llgE = sep_fittedE_llgE_fittedTotalE(2:end,3)*1e-1;

center2centerFactor = 110;
trueRij = (center2centerFactor+sep_fittedE_llgE_fittedTotalE(2:end,1))*1e-9;

figure, plot(trueRij,abs(fittedE*-1),'b*');
hold on
plot(trueRij, abs(fittedE*-1),'b');

plot(trueRij, llgE,'r*');
plot(trueRij, llgE,'r');

legend('abs-Fitted Coupling E','','LLG E','')
grid on
hold off;

%%
%% Delete me

Rij_1 = sep_fittedE_llgE_fittedTotalE(2:end,5);
aerijcos_1 = sep_fittedE_llgE_fittedTotalE(2:end,6);

gamma_a = -3.735e-05;
% gamma from manual calc
%gamma_a = -8.735e-05;

rho_b= -2.452e+05;

aeRij = gamma_a*exp(rho_b*Rij);
aeRij_DiDjCos = aeRij.*SiSjCos;
aeRij_Cos = aeRij.*Cos_Ttai_Ttaj;
