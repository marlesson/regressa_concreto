

P.powers <- function(modelo, p){
  return(array(modelo$k1*(1-p)^modelo$k2))
}

P.ryshkevitch <- function(modelo, p){
  return(array(modelo$k1*exp(modelo$k2*p)))
}

P.hasselmann <- function(modelo, p){
  return(array(modelo$k1+(modelo$k2*p)))
}

P.schiller <- function(modelo, p){
  return(array(modelo$k1*log(modelo$k2/p)))
}

P.proposto <- function(modelo, p){
  return(predict(modelo))  
}

P.plot   <- function(powers, ryshkevitch, hasselmann, schiller, proposto, dados){

  x   <- seq(0.10, 0.25, length.out=100)
  
  #Powers
  df  <- data.frame(
            func=rep('Powers', length(dados$POROSIDADE)), 
            x=x, y=P.powers(powers, x))
  
  #Ryshkevitch
  df  <- rbind(df, data.frame(
            func=rep('Ryshkevitch', length(dados$POROSIDADE)), 
            x=x, y=P.ryshkevitch(ryshkevitch, x)))
  
  #hasselmann
  df  <- rbind(df, data.frame(
            func=rep('Hasselmann', length(dados$POROSIDADE)), 
            x=x, y=P.hasselmann(hasselmann, x)))

  
  #schiller
  df  <- rbind(df, data.frame(
            func=rep('Schiller', length(dados$POROSIDADE)), 
            x=x, y=P.schiller(schiller, x)))
  
  
  #proposto
  #df  <- rbind(df, data.frame(
  #  func=rep('Proposto', length(dados$POROSIDADE)), 
  #  x=x, y=P.proposto(proposto, x)))
  
  
  # Pontos de coleta
  data_points <- data.frame(func='CF', x = dados$POROSIDADE, y = dados$COMPRESSAO)    
  
  p <- ggplot(data=df, aes(x=x*100, y=y, group=func, colour=func)) +
    geom_line()+
    scale_fill_brewer(palette="Paired", labels=c("Powers", "Ryshkevitch", "Hasselmann", "Schiller"))+
    geom_point(data=data_points, colour="black")+
    labs(title="", 
         x="Porosidade (%)", y = "Resistência à Compress\u{E3}o (MPa)", colour="Equa\u{00E7}\u{00E3}o")+
    theme(axis.text =element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=20,face="bold"))+
    theme_minimal()
  
  print(p)  
  
  #ggsave(file=filename)
  
}