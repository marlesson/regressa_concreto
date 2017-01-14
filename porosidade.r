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
dados <- read.table(file="dados_porosidade.txt", header=T, sep=" ")

## Filtra apenas o EVA fino
dados <- dados[dados$TIPO!='EG15'&dados$TIPO!='1EG15',]

## Dados de compressão do CR (Concreto de Referencia)
#fcc0  <- dados[dados$EVA==0&dados$FIBRA==0,]$COMPRESSAO

## Dados de Porosidade
y   <- dados$POROSIDADE/100


# Porosidade 
porosidade <- function(modelo, eva, fib, ms, mi){
  k1 <- coef(modelo)['k1']
  k2 <- coef(modelo)['k2']
  k3 <- coef(modelo)['k3']
  
  return((fib*k1+eva*k2+k3)/(fib*k1+eva*k2+k3+ms-mi))
}

######################### Rergessao  ############


### MODELO ABSORÇÃO

der    <- deriv3(~((FIBRA*k1+EVA*k2+k3)/(FIBRA*k1+EVA*k2+k3+MS-MI)), 
                 c("k1", "k2", "k3"), 
                 function(FIBRA, k1, EVA, k2, k3,MS,MI){ NULL })

modelo <- nls(y~der(FIBRA, k1, EVA, k2, k3,MS,MI), data=dados, start=list(k1=1, k2=1, k3=1))

y_r    <- predict(modelo,dados)
          
summary(modelo)
paste("R^2 -> ",R2(y_r, y))

porosidade(modelo, 0, 0, dados$MS, dados$MI)
# GRAFICOS

PrintResiduoNormalizado(pathImagens("porosidade", "residuo_normalizado.png"), 
                        modelo, residuals(modelo), y_r)

# Funcao
printModeloAjuste(pathImagens("porosidade", "modelo_ajuste.png"), 
                  modelo, y_r, y, ylim=c(0.10, 0.25), xlim=c(0.10, 0.25))

