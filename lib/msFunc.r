# Funções de auxilio para a regressão e impressão dos gráficos
#
# Marlesson Santana
#


# Path de imagens para salvar
pathImagens <- function(tipo, filename){
  return(paste("Graficos/",tipo, "/", filename, sep = ""))
}

# Definicao do modelo, 
M.predict <- function(modelo, eva, fibra){
  return(array(modelo$fcc0+predict(modelo, newdata = data.frame(EVA=eva, FIBRA=fibra))))
}

# Calcula o R2
Q.R2 <- function(y_pred, y, ...){
  
  sqt  <- sum((y     - median(y))^2)
  
  sqr  <- Q.SSE(y_pred, y)
  
  return( 1 - (sqr / sqt) )
}

# Calcula o R2
Q.R2Pred <- function(y, press){
  
  sqt  <- sum((y - median(y))^2)
  
  return( 1 - (press / sqt) )
}

# Calcula o soma dos quadrados do resuduo
Q.SSE <- function(y_pred, y, ...){
  return( sum((y_pred - y)^2) )
}

# Calcula o soma dos quadrados do resuduo
Q.MSE <- function(modelo, y_pred, y, ...){
  
  if(is.na(coef(modelo)["(Intercept)"])){
    coefs <- length(coef(modelo))
  }else{
    coefs <- length(coef(modelo)) - 1
  }
  
  return(Q.SSE(y_pred, y)/(length(y)-coefs-1))
}

# Calcula o PRESS
Q.PRESS <- function(modelo, y_str, fcc0, dados){
  press <- 0

  for(i in 1:length(dados$TIPO)){
    newDados       <- dados[-i, ]
    
    ## Calcula o y da regressao
    newDados['y']  <- newDados[y_str] - fcc0
      
    newModelo      <- update(modelo, data=newDados) 
    newModelo$fcc0 <- fcc0
    
    newPred        <-  M.predict(newModelo, dados[i, 'EVA'], dados[i, 'FIBRA'])
    realValue      <-  dados[i, y_str]

    press <- press + (newPred- realValue)^2
  }
  
  return(press)
}


# Plot Residuo Normalizado
PrintResiduoNormalizado <- function(filename, modelo, residuos, fcc_r, ...){
  sqe           <- sum(residuos^2)
  qme           <- sqe/(length(resid(modelo))-length(coef(modelo))-1)
  
  nres          <- residuos/sqrt(qme)
  
  df    <- data.frame(res=nres, fcc=fcc_r)
  
  p <- ggplot(df, aes(x = fcc, y = res)) +  ylim(-2.5, 2.5) +
    geom_point(size=3.5)+
    geom_abline(slope = 0, intercept = 2, colour = "darkred", linetype = 2, size=0.6) + 
    geom_abline(slope = 0, intercept = 0, colour = "darkred", linetype = 2, size=0.7) + 
    geom_abline(slope = 0, intercept = -2, colour = "darkred", linetype = 2, size=0.7) + 
    labs(x="Valores Ajustados", y = "Res\u{ED}duo Normalizado")+
    theme_minimal()+
    theme(axis.text =element_text(size=15),
          axis.title=element_text(size=18),
          plot.title=element_text(size=20, face="bold"))
  
  #print(p)
  
  ggsave(file=filename)
}


# Plot Residuo EVA
PrintResiduoVarReg <- function(filename, x_lbl, residuos, var, ...){
  
  df    <- data.frame(residuos=residuos, var=var)
  
  p <- ggplot(df, aes(x = var, y = residuos))+
    geom_point(size=3.5)+
    geom_abline(slope = 0, intercept = 0, colour = "darkred", linetype = 2) + 
    labs(x=x_lbl, y = "Res\u{ED}duos")+
    theme_minimal()+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))
  
  print(p)
}


# Print Grafico da função 
printModeloAjuste  <- function(filename, modelo, fcc_r, fcc, ylim, xlim, ...){

  #-----------------------------------------------------------------------------
  # ajuste Pontos
  pred    <- data.frame(y=fcc,x=fcc_r)
  m0      <- lm(y~x, data = pred)
  
  #-----------------------------------------------------------------------------
  # fazendo a predição intervalar num grid regular mais fino
  seq_pred <- seq(min(fcc_r), max(fcc_r), length=10)
  pred2    <- data.frame(x=seq_pred)
  pred$x2  <- seq_pred
  pred$y2  <- predict(m0, newdata=pred2)
  
  #str(pred)
  r2      <- format(c(coef(m0), summary(m0)$r.squared)[3]*100, digits=4)

  #-----------------------------------------------------------------------------
  #Correla\u{E7}\u{00E3}o do Modelo
  ggplot(pred, aes(x = x, y = y)) +  ylim(ylim) + xlim(xlim) +
    geom_point(size=3)+
    geom_line(aes(y=y2, x=x2), colour = "darkred", linetype = 2)+
    geom_text(data = NULL, x =xlim[2], y = ylim[2], 
              label = paste("R\u{00B2}", r2), size=6)+
    labs(title="Valores Ajustados X Valores Observados", x="Valores Ajustados (MPa)", y = "Valores Observados (MPa)")+
    theme_minimal()+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))
  
  ggsave(file=filename)
  
}


# Print grafico 3D Eficiencia
print3DEfici <- function(filename, modelo, screen, lbl, lim){
  
  # Limite 
  max_eva  <- lim[1]
  max_fib  <- lim[2]
  
  pred       <- expand.grid(bloco="I", EVA=seq(0,max_eva,l=35), FIBRA=seq(0,max_fib,l=35))
  pred$FCC   <- ((M.predict(modelo,pred$EVA, pred$FIBRA)-M.predict(modelo,pred$EVA,0*pred$FIBRA))/M.predict(modelo,pred$EVA,0*pred$FIBRA))*100
  pred$EVA   <- pred$EVA*100
  pred$FIBRA <- pred$FIBRA*100
  pred <- pred[pred$FCC > 0,]
  #-----------------------------------------------------------------------------
  # vendo a superfície ajustada
  args(panel.3d.contour)
  body(panel.3d.contour)
  
  colr <- brewer.pal(11, "RdYlGn")
  colr <- colorRampPalette(colr, space="rgb")
  
  xlab <- lbl[1]
  ylab <- lbl[2]
  zlab <- lbl[3]
  
  png(file=filename, width=400, height=400)
  fig <- wireframe(FCC~EVA*FIBRA, data=pred,
                   scales  =list(arrows=FALSE, axis=list(text=list(cex=1.5))), 
                   main    =list(cex=1.5),
                   zlab    =list(zlab, rot=90, cex=1.2),
                   xlab    =list(xlab, rot=24, cex=1.2), 
                   ylab    =list(ylab, rot=-37, cex=1.2),
                   cex.lab =0.2,
                   zlim    =c(lim[3],lim[4]), 
                   col     ="gray50", 
                   col.contour=1,
                   panel.3d.wireframe="panel.3d.contour", 
                   type=c("on","bottom"),
                   col.regions=colr(100),  
                   drape=TRUE,
                   screen=screen)
  print(fig)
  #
  dev.off()
}


# Print grafico 3D
print3DFunc <- function(filename, modelo, screen, lbl, lim){

  # Limite 
  max_eva  <- lim[1]
  max_fib  <- lim[2]
  
  pred       <- expand.grid(bloco="I", EVA=seq(0,max_eva,l=35), FIBRA=seq(0,max_fib,l=35))
  pred$FCC   <- M.predict(modelo,pred$EVA,pred$FIBRA)
  pred$EVA   <- pred$EVA*100
  pred$FIBRA <- pred$FIBRA*100
  
  #-----------------------------------------------------------------------------
  # vendo a superfície ajustada
  args(panel.3d.contour)
  body(panel.3d.contour)
  
  colr <- brewer.pal(11, "RdYlGn")
  colr <- colorRampPalette(colr, space="rgb")
  
  xlab <- lbl[1]
  ylab <- lbl[2]
  zlab <- lbl[3]
  
  #png(file=filename, width=850, height=850)
  fig <- wireframe(FCC~EVA*FIBRA, data=pred,
            scales  =list(arrows=FALSE, axis=list(text=list(cex=2.5))), 
            main    =list(cex=2),
            zlab    =list(zlab, rot=90, cex=2),
            xlab    =list(xlab, rot=24, cex=2), 
            ylab    =list(ylab, rot=-37, cex=2),
            cex.lab =1.2,
            zlim    =c(lim[3],lim[4]), 
            col     ="gray50", 
            col.contour=1,
            panel.3d.wireframe="panel.3d.contour", 
            type=c("on","bottom"),
            col.regions=colr(100),  
            drape=TRUE,
            screen=screen)
  print(fig)
  #
  #dev.off()
}


# Cuervas de Nivel G3d
printCurvasNivelFunc <- function(filename, modelo, lbl, lim){
  
  # Limite 
  max_eva  <- lim[1]
  max_fib  <- lim[2]
  
  colr <- brewer.pal(11, "RdYlGn")
  colr <- colorRampPalette(colr, space="rgb")
  
  xlab <- lbl[1]
  ylab <- lbl[2]
  zlab <- lbl[3]
  
  
  pred     <- expand.grid(bloco="I", EVA=seq(0,max_eva,l=200), FIBRA=seq(0,max_fib,l=200))
  pred$FCC <- M.predict(modelo, pred$EVA, pred$FIBRA)
  pred$EVA   <- pred$EVA*100
  pred$FIBRA <- pred$FIBRA*100
  
  #png(file=filename, width=950, height=850)
  p <- levelplot(FCC~EVA*FIBRA, data=pred, col.regions=colr(100),
                 xlab = list(label=xlab, cex=2),
                 ylab = list(label=ylab, cex=2),
                 panel= function(..., at, contour=FALSE, labels=NULL){
                   panel.levelplot(..., at=at, contour=contour, labels=labels)
                   panel.contourplot(..., at=at, contour=TRUE,
                                     labels=list(labels=format(at, digits=3), cex=2))
                 },
                 par.settings=list(
                   layout.widths=list(right.padding=4, cex=2)))
  
  p$legend$right <- list(fun=mergedTrellisLegendGrob(p$legend$right,
                                                     list(fun =textGrob, 
                                                          args=list(zlab, rot=-90, x=2)),
                                                     vertical=FALSE))
  print(p)
  #dev.off()
}

# MAximizando modelo para utilização da fibra
optimizeFibModelo <- function(modelo, iEva, maximum){
  
  # Intervalo de utilização da fibra
  iFib <- c(-0.0001, 0.0201) 
  
  
  fcc_sem_fibra  <- c()
  fcc_com_fibra  <- c()
  fib            <- c() 

  i <- 1
  for(eva in iEva){
    fibMax <- optimize(M.predict, iFib, tol = 0.0001, maximum=maximum, modelo = modelo,  eva = eva)
    
    fcc_sem_fibra[i] <- M.predict(modelo,eva,0)
    if(maximum){
      fcc_com_fibra[i] <- M.predict(modelo,eva,fibMax$maximum)
      fib[i]           <- fibMax$maximum    
    }else{
      fcc_com_fibra[i] <- M.predict(modelo,eva,fibMax$minimum)
      fib[i]           <- fibMax$minimum    
    }
    
    i <- i+1
  }
  
  return(list(eva =iEva, fib=fib,fcc_sem_fibra=fcc_sem_fibra,fcc_com_fibra=fcc_com_fibra))
}

# Plot de Barra do máximo de resistencia com a utlização
# da fibra
printMaxResistenciaFibra <- function(filename, maxFib, lbl_y){
  #tipo eva    fib           len      fcc
  #1   SF   0         2.108000e+01 21.08000
  #2   CF   0        -6.107621e-05 21.07994
  #3   SF   5         2.061578e+01 20.61578
  #4   CF   5 0.26 %  9.604657e-02 20.71182
  #5   SF  10         1.922310e+01 19.22310
  #6   CF  10 0.52 %  3.841863e-01 19.60729
  
  df       <- data.frame()

  # Melhora em porcentagem
  pMelhora <- round(((maxFib$fcc_com_fibra-maxFib$fcc_sem_fibra)/maxFib$fcc_sem_fibra)*100,2)
  
  
  for (i in seq(1,length(maxFib$eva))) {
    if(maxFib$fib[i] > 1e-3){
      fib     <- ""
      melhora <- paste(round(maxFib$fib[i]*100, 2),"% Fibra\n","+",pMelhora[i], "%", sep = "")
    }else{
      fib     <- ""
      melhora <- "Sem\nMelhora"
    }
    
     df <- rbind(df, data.frame(tipo='SF', 
                               eva=maxFib$eva[i]*100, 
                               fib=fib,
                               cMelhora='',
                               len = maxFib$fcc_sem_fibra[i],
                               fcc = maxFib$fcc_sem_fibra[i]))
    
  
    
    df <- rbind(df, data.frame(tipo='CF', 
                               eva=maxFib$eva[i]*100, 
                               fib='',
                               cMelhora=melhora,
                               len= maxFib$fcc_com_fibra[i]-maxFib$fcc_sem_fibra[i],
                               fcc= maxFib$fcc_com_fibra[i])) 
    
  }
  
 
  p <- ggplot(data=df, aes(x=eva, y=len, fill=tipo)) +
    geom_bar(stat="identity")+
    geom_label(aes(y=fcc+0.2,  label=cMelhora), size=4.2, nudge_x = 0.5)+
    geom_text(aes(y=len+0.05, label=fib), size=4.2, nudge_x = 0.5)+
    scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
    theme_minimal()+
    labs(x="EVA (%)", y = lbl_y, fill="Tipo")+
    theme(axis.text =element_text(size=16),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))
  
  print(p)
  #ggsave(file=filename)
}

# Retorna um vetor de indices para os pontos de interceptação entre
# as fuas retas
interceptLine <- function(vet1, vet2, error){
  
  i = c()
  
  for(k in seq(1, length(vet1))){
    sub <- abs(vet1[k]-vet2[k])
    
    i[k] <- (sub < error) 
  }
  
  return(i)
}

# Plot de linha do máximo de resistencia com a utlização
# da fibra com linha sem utilização da fibra
printMaxResistenciaFibraLine <- function(filename, maxFib, lbl_y){
  df     <- data.frame()
  
  for (i in seq(1,length(maxFib$eva))) {
    df <- rbind(df, data.frame(tipo='CF', 
                               eva=maxFib$eva[i]*100, 
                               fcc= maxFib$fcc_com_fibra[i])) 
    
    
    df <- rbind(df, data.frame(tipo='SF', 
                           eva=maxFib$eva[i]*100, 
                           fcc = maxFib$fcc_sem_fibra[i]))
  }
  

  p <- ggplot(data=df, aes(x=eva, y=fcc, group=tipo, colour=tipo)) +
    geom_line()+
    scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
    theme_minimal()+
    labs(title="Maximiza\u{E7}\u{E3}o da Resist\u{00EA}ncia com a Utiliza\u{E7}\u{E3}o da Fibra", 
         x="EVA (%)", y = lbl_y, colour="Tipo")+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))

  #Pontos de divisão
  intercept_lines <- interceptLine(maxFib$fcc_sem_fibra, maxFib$fcc_com_fibra, error=1e-02)
  
  # Se tiver intercecção com as linhas
  if(length(which(intercept_lines)) > 0){
    i_max   <- max(which(intercept_lines))
    
    x1 <- maxFib$eva[i_max]*100
    y1 <- maxFib$fcc_com_fibra[i_max]
    points_inter <- data.frame(tipo='CF',eva = x1, fcc = y1)    

    # Plot point
    p <- p + geom_point(data=points_inter, colour="red") +
      geom_text(data=head(points_inter, 1), 
                label=paste("(",round(x1, 2),", ", round(y1,2),")", sep = ""), 
                colour = "black", size=4, hjust=1)    
  }
  
  
  #print(p)  
  
  ggsave(file=filename)
}

# Plot de linha com valores de fibra 0, 1 e 2 
printResistenciaFibraLinha <- function(filename, modelo, lbl_y){
  iEva     <- seq(0, 0.25, by=0.0025)

  fcc_0    <- M.predict(modelo,iEva, rep(0,    length(iEva)))
  fcc_1    <- M.predict(modelo,iEva, rep(0.01, length(iEva)))
  fcc_2    <- M.predict(modelo,iEva, rep(0.02, length(iEva)))
  
  df     <- data.frame()
  
  for (i in seq(1,length(iEva))) {
    df <- rbind(df, data.frame(tipo='SF', 
                               eva = iEva[i]*100, 
                               fcc = fcc_0[i]))
    
    df <- rbind(df, data.frame(tipo='CF1', 
                               eva=iEva[i]*100, 
                               fcc=fcc_1[i]))
    
    
    df <- rbind(df, data.frame(tipo='CF2', 
                               eva=iEva[i]*100, 
                               fcc=fcc_2[i])) 
  }
  
  
  p <- ggplot(data=df, aes(x=eva, y=fcc, group=tipo, colour=tipo)) +
    geom_line()+
    scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
    theme_minimal()+
    labs(title="Resist\u{00EA}ncia X EVA", x="EVA (%)", 
         y = lbl_y, colour="Tipo")+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))
  
  
  # Verifica se tem ponto de intercepcao
  intercept_lines <- interceptLine(fcc_0, fcc_1, error=1e-01)
  
  # Se tiver intercecção com as linhas fcc_0, fcc_1
  if(length(which(intercept_lines)) > 0){
    i_max   <- median(which(intercept_lines))
    
    x1 <- iEva[i_max]*100
    y1 <- fcc_0[i_max]
    points_inter <- data.frame(tipo='CF2',eva = x1, fcc = y1)    
    
    # Plot point
    p <- p + geom_point(data=points_inter, colour="red") +
      geom_text(data=head(points_inter, 1), 
                label=paste("(",round(x1, 2),", ", round(y1,2),")", sep = ""), 
                colour = "black", size=4, hjust=1)    
  }

  intercept_lines <- interceptLine(fcc_0, fcc_2, error=1e-01)
  
  # Se tiver intercecção com as linhas fcc_0, fcc_1
  if(length(which(intercept_lines)) > 0){
    i_max   <- median(which(intercept_lines))
    
    x1 <- iEva[i_max]*100
    y1 <- fcc_0[i_max]
    points_inter <- data.frame(tipo='CF2',eva = x1, fcc = y1)    
    
    # Plot point
    p <- p + geom_point(data=points_inter, colour="red") +
      geom_text(data=head(points_inter, 1), 
                label=paste("(",round(x1, 2),", ", round(y1,2),")", sep = ""), 
                colour = "black", size=4, hjust=1)    
  }
  
  #print(p)
  
  ggsave(file=filename)
}


# Plot de linha com valores de eva 5, 15, 25
printResistenciaEVALinha <- function(filename, modelo, lbl_y){
  iFib     <- seq(0, 0.02, by=0.00025)
  
  fcc_0    <- M.predict(modelo, rep(0.05, length(iFib)), iFib)
  fcc_1    <- M.predict(modelo, rep(0.15, length(iFib)), iFib)
  fcc_2    <- M.predict(modelo, rep(0.25, length(iFib)), iFib)
  
  df     <- data.frame()
  
  for (i in seq(1,length(iFib))) {
    df <- rbind(df, data.frame(tipo='EVA5', 
                               fib = iFib[i]*100, 
                               fcc = fcc_0[i]))
    
    df <- rbind(df, data.frame(tipo='EVA15', 
                               fib=iFib[i]*100, 
                               fcc=fcc_1[i]))
    
    
    df <- rbind(df, data.frame(tipo='EVA25', 
                               fib=iFib[i]*100, 
                               fcc=fcc_2[i])) 
  }
  
  #Pontos de intercecção
  y  <- min(abs(fcc_0-fcc_1))
  i  <- which(abs(fcc_0-fcc_1) == y)
  x1 <- iFib[i]*100
  y1 <- fcc_0[i]
  
  y <- min(abs(fcc_2-fcc_0))
  i  <- which(abs(fcc_2-fcc_0) == y)
  x2 <- iFib[i]*100
  y2 <- fcc_0[i]
  
  p <- ggplot(data=df, aes(x=fib, y=fcc, group=tipo, colour=tipo)) +
    geom_line()+
    scale_fill_brewer(palette="Paired", labels=c("EVA", "EVA+Fibra"))+
    theme_minimal()+
    labs(title="Resist\u{00EA}ncia X Fibra", x="Fibra (%)", 
         y = lbl_y, colour="Tipo")+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))

  
  # print(p)
  
  ggsave(file=filename)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}