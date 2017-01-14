library(plotrix)
library(lmtest)
library(qpcR)
library(nlstools)

require(qpcR) 
require(rpanel)
require(lattice)
require(latticeExtra)
require(grid)
require(lattice)
require(RColorBrewer)
library(ggplot2)

source("lib/nlsResiduals.r")
source("lib/ridiculas_functions.r")
source("lib/msFunc.r")


# ---------------------------------
# Definição dos modelos encontrados
# 

# Compressão
MOD.fc <- function(eva, fib){
  return(21.08*(1-(8.809*eva^2+681.187*fib^2-70.469*eva*fib)))
}

# Tração
MOD.ft <- function(eva, fib){
  return(3.23*(1-(2.313*eva+11.553*fib-91.437*eva*fib)))
}

# Deformacao
MOD.ec <- function(eva, fib){
  return(36.71*(1-(1.444*eva+201.752*fib^2)))
}

# Maximização da Fibra (Compressão)
MOD.max_fib <- function(eva){
  return ((1485.49/28718.8)*eva)
}
# ---------------------------------


### DADOS
dados <- read.table(file="dados.txt", header=T, sep=" ")

## Filtra apenas o EVA fino
dados <- dados[dados$TIPO!='EG15'&dados$TIPO!='1EG15',]


paste("R^2 -> ",R2(MOD.fc(dados$EVA, dados$FIBRA), dados$COMPRESSAO))

paste("R^2 -> ",R2(MOD.ft(dados$EVA, dados$FIBRA), dados$TRACAO))

paste("R^2 -> ",R2(MOD.ec(dados$EVA, dados$FIBRA), dados$MODULO))

## PONTOS A SEREM PLOTADOS

iEva=seq(0, 0.25, by=0.005)

######################## COMPRESSAO ################################
paste("R^2 -> ",R2(MOD.fc(dados$EVA, dados$FIBRA), dados$COMPRESSAO))

df     <- data.frame()

for (i in seq(1,length(iEva))) {
  df <- rbind(df, data.frame(tipo='CF', 
                             eva=iEva[i]*100, 
                             fc= MOD.fc(iEva[i], MOD.max_fib(iEva[i])))) 
  
  
  df <- rbind(df, data.frame(tipo='SF', 
                             eva=iEva[i]*100, 
                             fc= MOD.fc(iEva[i], 0))) 
}


p1 <- ggplot(data=df, aes(x=eva, y=fc, group=tipo, colour=tipo)) +
  geom_line()+
  scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
  theme_minimal()+
  labs(title="Maximiza\u{E7}\u{E3}o da Resist\u{00EA}ncia com a Utiliza\u{E7}\u{E3}o da Fibra", 
       x="EVA (%)", y = "Resistência à Compress\u{E3}o (MPa)", colour="Tipo")+
  theme(axis.text =element_text(size=12),
        axis.title=element_text(size=16),
        plot.title=element_text(size=20,face="bold"))

print(p1)  

######################## TRACAO ################################
paste("R^2 -> ",R2(MOD.ft(dados$EVA, dados$FIBRA), dados$TRACAO))

df     <- data.frame()

for (i in seq(1,length(iEva))) {
  df <- rbind(df, data.frame(tipo='CF', 
                             eva=iEva[i]*100, 
                             fc= MOD.ft(iEva[i], MOD.max_fib(iEva[i])))) 
  
  
  df <- rbind(df, data.frame(tipo='SF', 
                             eva=iEva[i]*100, 
                             fc= MOD.ft(iEva[i], 0))) 
}


p2 <- ggplot(data=df, aes(x=eva, y=fc, group=tipo, colour=tipo)) +
  geom_line()+
  scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
  theme_minimal()+
  labs(title="Maximiza\u{E7}\u{E3}o da Resist\u{00EA}ncia com a Utiliza\u{E7}\u{E3}o da Fibra", 
       x="EVA (%)", y = "Resistência à Trac\u{E3}o (MPa)", colour="Tipo")+
  theme(axis.text =element_text(size=12),
        axis.title=element_text(size=16),
        plot.title=element_text(size=20,face="bold"))

print(p2) 


######################## Deformação ################################
paste("R^2 -> ",R2(MOD.ft(dados$EVA, dados$FIBRA), dados$TRACAO))

df     <- data.frame()

for (i in seq(1,length(iEva))) {
  df <- rbind(df, data.frame(tipo='CF', 
                             eva=iEva[i]*100, 
                             fc= MOD.ec(iEva[i], MOD.max_fib(iEva[i])))) 
  
  
  df <- rbind(df, data.frame(tipo='SF', 
                             eva=iEva[i]*100, 
                             fc= MOD.ec(iEva[i], 0))) 
}


p3 <- ggplot(data=df, aes(x=eva, y=fc, group=tipo, colour=tipo)) +
  geom_line()+
  scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
  theme_minimal()+
  labs(title="Maximiza\u{E7}\u{E3}o da Resist\u{00EA}ncia com a Utiliza\u{E7}\u{E3}o da Fibra", 
       x="EVA (%)", y = "Modulo de Deformacao (GPa)", colour="Tipo")+
  theme(axis.text =element_text(size=12),
        axis.title=element_text(size=16),
        plot.title=element_text(size=20,face="bold"))

print(p3) 
#ggsave(file=filename)


multiplot(p1,p2,p3, cols=2)



