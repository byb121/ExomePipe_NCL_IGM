detLogLike<-function(VarFreqs,NoVarFreqs,Sens,Spez)
{
  # 04/05/2011
  # Edited 15/06/2011
  #2*s-2*s*f-s^2+s^2*f+f-u*f+u*f^2+s*u*f-s*u*f^2-u^2*f^2
  ObsVar<-2*Sens-2*Sens*VarFreqs-Sens^2+Sens^2*VarFreqs+VarFreqs-Spez*VarFreqs+Spez*VarFreqs^2+
          Sens*Spez*VarFreqs-Sens*Spez*VarFreqs^2-Spez^2*VarFreqs^2
  #  1-f-2*s+2*s*f+s^2-s^2*f+u*f-u*f^2-s*u*f+s*u*f^2+u^2*f^2     
  ObsNoVar<-1-(2*Sens-2*Sens*NoVarFreqs-Sens^2+Sens^2*NoVarFreqs+NoVarFreqs-Spez*NoVarFreqs    
              +Spez*NoVarFreqs^2+Sens*Spez*NoVarFreqs-Sens*Spez*NoVarFreqs^2-
              Spez^2*NoVarFreqs^2)
  sum(log(ObsVar))+sum(log(ObsNoVar))
}

estimateShell_III<-function(Observed,Freqs,Selected.SNPs,
                            Start.Spez=0.5,Start.Sens=0.5,Eps=1e-7)
{
  # 04/05/2011,derived from estimateShell_II (last edit 06/09/2010). 
  # Edited: 14/03/2013
  # Input:
  #   Observed: Vector (NbrOfSNPs): 
  #                           <1: NA (should be marked already in selected);
  #                            1: WT homoz; 
  #                            2: Het; 
  #                            3: Variant homoz.
  #   Freqs: Vector (NbrOfSNPs) with variant freqency.
  #   Selected.SNPs: (NbrSelectedSNPs) Indices (in vectors above) 
  #       of  SNPs to use
  detLikeInt<-function(x){-detLogLike(VarFreqs=CleanFreqs[[1]],NoVarFreqs=CleanFreqs[[2]],Sens=x[2],Spez=x[1])}
  Total.SNP.Nbr<-length(Observed)
  Selected<-rep(F,Total.SNP.Nbr)

  if (any(is.na(Selected.SNPs))) Selected[1:Total.SNP.Nbr]<-T
    else Selected[Selected.SNPs]<-T
  VarSNPs<-Observed==2 | Observed==3
  OrigSNPs<-Observed==1
  CleanFreqs<-list(VarFreqs=Freqs[Selected & VarSNPs],OrigFreqs=Freqs[OrigSNPs & Selected])
  x<-c(Start.Spez,Start.Sens)
  Res<-optim(par=x, detLikeInt, gr = NULL,method = "L-BFGS-B",lower=c(Eps,Eps),upper=c(1-Eps,1-Eps))
  list(Sens=Res$par[2],Spez=Res$par[1],LogLike=Res)
}

analyseData0320<-function(Freqs,Genotypes,MinFreq=1e-5,Eps_Freq=1e-5,
                          UpCompatible=NA,DownCompatible=NA,
                          SegmentSize=1,NBoot=1000)
{
  # Derived 20/03/2013 from analyseData0505(last edited 15/06/2011)
  # Edited          
  # Genotypes: Vector (NbrOfSNPs): 
  #         <1: Not covered;
  #         1: WT homoz; 
  #         2: Het; 
  #         3: Variant homoz. 
  # MinFreq: Minimum Variant Frequency used for calculation
  # Eps_Freq: Frequency considered undistuingishable from 0 
  #        (i.e. to be set to MinFreq)
  # UpCompatible (DownCompatible): Vector containing for each variant the 
  #        location (in the Genotypes and Freqs vectors) of the 
  #        nearest variant (going either towards higher(Up) or 
  #        lower(Down) indices in the order used to sort 
  #        Genotypes and Freqs), that is either on another 
  #        chromosome or distant enough(see MinIndDist in 
  #        prepareAnalysis.CI)
  # SegmentSize:  Parameter used in choosing the variants analysed at iteration
  #       (default 1) 
  # NBoot: Number of iterations to be analysed. In each iteration a subset of 
  #       the originial set of variants is used.  These subsets are chosen so 
  #       that any two variants of one subset are either on different chromosomes 
  #       or separated by at least MinIndDist bp.
  #
  Genotypes<-Genotypes[,1,drop=F]
  N.in.Target<-nrow(Genotypes)       
  if (all(is.na(UpCompatible))) UpCompatible<-c(1:N.in.Target)
  if (all(is.na(DownCompatible))) DownCompatible<-c(1:N.in.Target)
  Freqs[Freqs<Eps_Freq]<-MinFreq
  Freqs[Freqs>(1-Eps_Freq)]<-(1-MinFreq)
  Intermediate<-matrix(ncol=2,nrow=NBoot)
  Results<-matrix(nrow=2,ncol=3)
  dimnames(Results)<-list(c("Sensitivity","Specificity"),
                          c("Median","0.05","0.95"))
  Quantls<-c(0.5,0.05,0.95)
  for (i in seq(NBoot)){
       Indx<-rep(F,N.in.Target)
       FirstChoosen<-sample(N.in.Target,1)
       Choosen<-FirstChoosen
       repeat {
          if (Choosen<1) break;
          Indx[Choosen]<-T
          Boundary<-min(Choosen+SegmentSize,N.in.Target )
          Choosen<-sample(UpCompatible[Choosen:Boundary],1)
       }
       Choosen<-FirstChoosen
       repeat {
          if (Choosen<1) break;
          Indx[Choosen]<-T
          Boundary<-max(1,Choosen-SegmentSize)
          Choosen<-sample(DownCompatible[Boundary:Choosen],1)
       }
       Indx<-(Genotypes[,1]>0) & Indx
       RunningGenos<-Genotypes[Indx,1]
       RunningFreqs<-Freqs[Indx]
     OptiRes<- estimateShell_III(Observed=RunningGenos,Freqs=RunningFreqs,
                                   Selected.SNPs=NA,Start.Spez=0.5,Start.Sens=0.5)
     Intermediate[i,]<-c(OptiRes[[1]], OptiRes[[2]])
  }
  Results["Sensitivity",]<-quantile(Intermediate[,1],Quantls)
  Results["Specificity",]<-quantile(Intermediate[,2],Quantls)
  Results
}
prepareAnalysis.CI<-function(AllData,MinFreq=1e-5,Eps_Freq=1e-5,
                            MinIndDist=500000, NBoot=1000,
                            ExcludeXY=T,ExcludeImpossibles=F)
{
  # derived 20/03/2013 from prepareAnalysis (last edited 15/06/2011)
  # Edited 20/03/2013  
  # MinFreq: Minimum Variant Frequency used for calculation
  # Eps_Freq: Frequency considered undistuingishable from 0 
  #       (i.e. to be set to MinFreq)
  # MinIndDist: Minimal distance in bp between two variants that is considered
  #       sufficient to ensure that both variants are not in 
  #       linkage disequilibrium  (set to 0.5 Mbp)
  # NBoot: Number of iterations to be analysed. In each iteration a subset of 
  #       the originial set of variants is used.  These subsets are chosen so 
  #       that any two variants of one subset are either on different chromosomes 
  #       or separated by at least MinIndDist bp.(set to 1000).
  # ExcludeXY: Excludes entries on the X and Y chromosomes and on the 
  #       mitochondrial genome
  # ExcludeImpossible: Excludes entries which
  #       show code 1 in all samples but frequency is 0
  #       show code 3 in all samples but frequency is 1   
  # Genotypes coded as: 
  #       <1: Not covered;
  #       1: WT homoz; 
  #       2: Het; 
  #       3: Variant homoz. 
  IndxXYM<-rep(F,nrow(AllData))
  InvariantImpWT<-rep(F,nrow(AllData)) 
  InvariantImpVar<-rep(F,nrow(AllData)) 
  if (ExcludeXY) IndxXYM<-((AllData$Chromosome=="chrX")|  
                           (AllData$Chromosome=="chrY")|  
                           (AllData$Chromosome=="chrM"))
  if  (ExcludeImpossibles) {
    InvariantImpWT<-apply(AllData[,-c(1:5)],1,function(x){all(x==1)})&(AllData$Frequency< Eps_Freq )
    InvariantImpVar<-apply(AllData[,-c(1:5)],1,function(x){all(x==3)})&(AllData$Frequency>(1- Eps_Freq))
  }
  Indx<-(!IndxXYM) & (!InvariantImpWT) & (!InvariantImpVar)
  AllData<-AllData[Indx,]
  NSites<-nrow(AllData)
  Chrs<-sort(unique(AllData$Chromosome))
  MinOffset<-max(AllData$Position)+MinIndDist
  ChrOffSet<-c(0:(length(Chrs)-1))*MinOffset
  ModPos<-AllData$Position
  for (i in seq(length(Chrs))){
    ChrIndx<-AllData$Chromosome==Chrs[i]
    ModPos[ChrIndx]<-ModPos[ChrIndx]+ChrOffSet[i]
  }
  OrderIndx<-order(ModPos)
  AllData<-AllData[OrderIndx,]
  ModPos<-ModPos[OrderIndx]
  UpCompatible<-rep(0,NSites)
  DownCompatible<-rep(0,NSites)
  for (i in seq(NSites)){
    SufficientAppart<-(ModPos>=(ModPos[i]+MinIndDist))
    if ((any(SufficientAppart))) 
      UpCompatible[i]<-min(which(SufficientAppart))
    SufficientAppartDown<-(ModPos<=(ModPos[i]-MinIndDist))
    if ((any(SufficientAppartDown))) 
      DownCompatible[i]<-max(which(SufficientAppartDown))
  } 
  Genotypes<-AllData[,-c(1:5),drop=F]
  Freqs<-AllData$Frequency
  analyseData0320(Freqs,Genotypes,MinFreq=MinFreq,Eps_Freq=Eps_Freq,
                  UpCompatible=UpCompatible,DownCompatible=DownCompatible,
                  NBoot=NBoot)
}
                               
# 14/03/2013
library(stringr)
args<- commandArgs(trailingOnly = TRUE)
inputfile<-args[1]   
outputfile<-str_replace(inputfile,"Input.txt","SensSpec.txt")
Table<-read.table(inputfile,head=T)
Results<-prepareAnalysis.CI(AllData=Table)
write.table(Results[-3,], file=outputfile, append=F, quote=F, sep="\t", row.names=T, col.names=F)
