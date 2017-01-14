#==========================================================================================
# Funções desenvolvidas para matérias do blog "Ridículas <- dicas curtas sobre R"
# www.ridiculas.wordpress.com, por Walmes Zeviani
# envie sugestões para walmes<at>ufpr.br
#==========================================================================================

#------------------------------------------------------------------------------------------
# painel para gráficos com contornos

panel.3d.contour <- function(x, y, z, rot.mat, distance, type="on",
                             nlevels=20, zlim.scaled, col.contour=1, ...){
  clines <- contourLines(x, y, matrix(z, nrow=length(x), byrow=TRUE),
                         nlevels=nlevels)
  if(any(type%in%c("bottom"))){
    for(ll in clines){
      n <- ltransform3dto3d(rbind(ll$x, ll$y, zlim.scaled[1]),
                            rot.mat, distance)
      panel.lines(n[1,], n[2,], col=col.contour, lty=1, lwd=1)
    }}
  panel.3dwire(x, y, z, rot.mat, distance, zlim.scaled=zlim.scaled, ...)
  if(any(type%in%c("on"))){
    for(ll in clines){
      n <- ltransform3dto3d(rbind(ll$x, ll$y, ll$level), rot.mat, distance)
      panel.lines(n[1,], n[2,], col=col.contour, lty=1, lwd=1)
    }}
  if(any(type%in%c("top"))){
    for(ll in clines){
      n <- ltransform3dto3d(rbind(ll$x, ll$y, zlim.scaled[2]),
                            rot.mat, distance)
      panel.lines(n[1,], n[2,], col=col.contour, lty=1, lwd=1)
    }}
}

# exemplo de uso
# da <- expand.grid(x=seq(0,1,l=20), y=seq(0,1,l=20))
# da$z <- with(da, 10+1*x+0.5*y+1.2*x*y-0.9*x^2-0.05*y^2)
# require(lattice)
# wireframe(z~x+y, da, zlim=c(8,12), drape=TRUE, col="white", col.contour="red",
#           panel.3d.wireframe="panel.3d.contour", type=c("on","bottom"),
#           screen=list(z=40, x=-75))

#------------------------------------------------------------------------------------------
# painéis para fazer curvas de regressão com bandas de confiança

prepanel.ciH <- function(x, y, ly, uy, subscripts, ...){
  x <- as.numeric(x)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  list(ylim=range(uy, ly, finite=TRUE), xlim=range(x))
}

panel.ciH <- function(x, y, ly, uy, subscripts, ...){
  y <- as.numeric(y)
  x <- as.numeric(x)
  or <- order(x)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  panel.polygon(c(x[or], x[rev(or)]),
                c(ly[or], uy[rev(or)]),
                col=1, alpha=0.05, border=NA)
  panel.lines(x[or], ly[or], lty=3, lwd=0.5, col=1)
  panel.lines(x[or], uy[or], lty=3, lwd=0.5, col=1)
  panel.xyplot(x, y, subscripts=subscripts, ...)
}

# exemplo de uso
#require(lattice); require(latticeExtra)
#da <- expand.grid(A=gl(2,4), B=0:10)
#set.seed(123); da$y <- with(da, rnorm(A, mean=as.numeric(A)*B, sd=2))
#m0 <- lm(y~A*B, da)
#da <- cbind(da, predict(m0, interval="confidence"))
#str(da)

## sem grupos
#xyplot(y~B|A, da)+
#  as.layer(xyplot(fit~B|A, da, type="l", ly=da$lwr, uy=da$upr,
#                  prepanel=prepanel.ciH, panel=panel.ciH))

## com grupos
#xyplot(y~B, groups=A, da)+
#  as.layer(xyplot(fit~B, groups=A, da, type="l",
#                  ly=da$lwr, uy=da$upr, prepanel=prepanel.ciH,
#                  panel=panel.superpose, panel.groups=panel.ciH))

#------------------------------------------------------------------------------------------
