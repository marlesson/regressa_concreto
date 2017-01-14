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

library("lattice") 
library(latticeExtra) # for mergedTrellisLegendGrob()

source("lib/nlsResiduals.r")
source("lib/ridiculas_functions.r")
source("lib/msFunc.r")


# Função principal para a geração dos gráficos
run_all <- function(modelo, dados){
  
  # vendo a superfície ajustada
  args(panel.3d.contour)
  body(panel.3d.contour)
  
  fcc_r    <- M.predict(modelo, dados$EVA, dados$FIBRA)
  fcc      <- dados$TRACAO
  
  residuos <- fcc-fcc_r
  
  reducao_res <- (M.predict(modelo, 0.25, 0)-M.predict(modelo, 0, 0))/M.predict(modelo, 0, 0)
  print(summary(modelo))
  print(paste("R^2 -> ", Q.R2(fcc_r, fcc)))
  print(paste("SSE -> ", Q.SSE(fcc_r, fcc)))
  press <- Q.PRESS(modelo, 'TRACAO', fcc0, dados)
  print(paste("PRESS -> ", press))
  print(paste("R^2 PRESS -> ", Q.R2Pred(fcc, press)))
  
  ## TEstes de hipotese
  
  #Shapiro-wilk
  print(shapiro.test(residuals(modelo)))
  
  
  #Durbin-Watson test
  print(dwtest(modelo))
  
  #Breusch-Pagan test
  print(bptest(modelo))

  
  ## PLOT  ##
  
  # Residuo
  PrintResiduoNormalizado(pathImagens("tracao", "residuo_normalizado.png"), 
                          modelo, residuos, fcc_r)
  
  PrintResiduoVarReg(pathImagens("tracao", "residuo.png"),
                     "Resistência à Tração", residuos, dados$TRACAO)
  
  
  PrintResiduoVarReg(pathImagens("tracao", "residuo_eva.png"),
                     "EVA (%)", residuos, dados$EVA)
  
  PrintResiduoVarReg(pathImagens("tracao", "residuo_fibra.png"), 
                     "Fibras (%)", residuos, dados$FIBRA)
  
  # Funcao
  printModeloAjuste(pathImagens("tracao", "modelo_ajuste.png"), 
                    modelo, fcc_r, fcc, ylim=c(1,4), xlim=c(1,4))
  
  # Graficos 3d
  lim_xyz <- c(0.26, 0.021, 0, 4)
  lbl_xyz <- c("EVA (%)", "Fibra (%)", "Resistência à Tra\u{E7}\u{E3}o (MPa)")
  
  print3DFunc(pathImagens("tracao", "3d_1.png"), 
              modelo, 
              list(z=-50, x=-70, y=0), 
              lbl_xyz,  
              lim_xyz)
  
  print3DFunc(pathImagens("tracao", "3d_2.png"), 
              modelo, 
              list(z=-120, x=-70, y=0), 
              lbl_xyz,  
              lim_xyz)
  
  # Graficos 3d EFICIENCIA
  print3DEfici(pathImagens("tracao", "3d_ef.png"), 
               modelo, 
               list(z=140, x=-70, y=-00), 
               c("EVA (%)", "Fibra (%)", "Eficiência (%)"),  
               c(0.26, 0.021, 0, 60))
  
  # Curvas de Nível
  printCurvasNivelFunc(pathImagens("tracao", "curva_nivel.png"), modelo, lbl_xyz,lim_xyz)
  
  # Max resistencia Barra
  maxFib <- optimizeFibModelo(modelo, iEva=seq(0.05, 0.25, by=0.05), TRUE)
  printMaxResistenciaFibra(pathImagens("tracao", "barra_max.png"), maxFib, lbl_xyz[3])
  
  # Max resistencia Evolucao
  maxFib <- optimizeFibModelo(modelo, iEva=seq(0,0.25, by=0.005), TRUE)
  printMaxResistenciaFibraLine(pathImagens("tracao", "eva_max.png"), maxFib, lbl_xyz[3])
  
  # Resistencia por eva, com 3 curvas 0,1,2 % fibra
  printResistenciaFibraLinha(pathImagens("tracao", "eva_0_1_2.png"), modelo, lbl_xyz[3])
  
  # Resistencia por fibra, com 3 curvas 5, 15, 25 % fibra
  printResistenciaEVALinha(pathImagens("tracao", "fibra_5_15_25.png"), modelo, lbl_xyz[3])
  
}


### DADOS
dados <- read.table(file="dados.txt", header=T, sep=" ")

## Filtra apenas o EVA fino
dados <- dados[dados$TIPO!='EG15'&dados$TIPO!='1EG15',]

## Dados de compressão do CR (Concreto de Referencia)
fcc0  <- dados[dados$EVA==0&dados$FIBRA==0,]$TRACAO

## Dados de Compressao 
fcc   <- dados$TRACAO

## Calcula o y da regressao
y     <- fcc-fcc0


######################### Rergessao  #########################

### MODELO 1 COMPLETO
modelo      <-lm(y~0+EVA+FIBRA+EVA:FIBRA+I(EVA^2)+I(FIBRA^2),data=dados) 
modelo$fcc0 <- fcc0

summary(modelo)
anova(modelo)

fcc_r    <- M.predict(modelo, dados$EVA, dados$FIBRA)
run_all(modelo, dados)


### MODELO 1 RESUMO
modelo      <-lm(y~0+EVA+FIBRA+EVA:FIBRA,data=dados) 
modelo$fcc0 <- fcc0

summary(modelo)
anova(modelo)

fcc_r    <- M.predict(modelo, dados$EVA, dados$FIBRA)
run_all(modelo, dados)
