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
source("lib/func_porosidade.r")

# Função principal para a geração dos gráficos
run_all <- function(modelo, dados){
  residuos <- fcc-fcc_r
  
  print(summary(modelo))
  print(paste("R^2 -> ",R2(fcc_r, fcc)))
  
}

### DADOS
dados <- read.table(file="dados.txt", header=T, sep=" ")

## Filtra apenas o EVA fino
dados <- dados[dados$TIPO!='EG15'&dados$TIPO!='1EG15',]

## Dados de compressão do CR (Concreto de Referencia)
fcc0  <- dados[dados$EVA==0&dados$FIBRA==0,]$COMPRESSAO

## Dados de Compressao
y   <- dados$COMPRESSAO

## Transforma dados porosidade
dados$POROSIDADE <- dados$POROSIDADE/100

######################### Rergessao  #########################

### MODELO - Balshin (powers)
powers     <- lm(log(y)~log(1-dados$POROSIDADE))
powers['k1']  <- exp(coef(powers)[1]) #Intercept
powers['k2']  <- coef(powers)[2]   
  
### MODELO - Ryshkevitch
ryshkevitch     <- lm(log(y)~dados$POROSIDADE)
ryshkevitch['k1']  <- exp(coef(ryshkevitch)[1]) #Intercept
ryshkevitch['k2']  <- coef(ryshkevitch)[2]   

### MODELO - Hasselmann
hasselmann     <- lm(y~dados$POROSIDADE)
hasselmann['k1']  <- coef(hasselmann)[1] #Intercept
hasselmann['k2']  <- coef(hasselmann)[2]   


### MODELO - Schiller
schiller     <- lm(y~log(dados$POROSIDADE))
schiller['k1']  <- coef(schiller)[2]*-1
schiller['k2']  <- exp(coef(schiller)[1]/schiller$k1)


### MODELO - PROPOSTO
proposto     <- lm(y~dados$POROSIDADE*dados$EVA*dados$FIBRA)
proposto['k1']  <- coef(proposto)[1] #Intercept
proposto['k2']  <- coef(proposto)[2]   


der    <- deriv3(~(k1*(1-POROSIDADE)^k2)*(1-(k3*EVA)),
                 c("k1", "k2", "k3"), 
                 function(POROSIDADE, FIBRA, EVA, k1, k2, k3){ NULL })

modelo <- nls(y~der(POROSIDADE, FIBRA, EVA, k1, k2, k3), data=dados, start=list(k1=1, k2=1, k3=1))

summary(modelo)
P.plot(powers, ryshkevitch, hasselmann, schiller, proposto, dados)

