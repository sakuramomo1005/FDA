
# Please run the main code first
######################################################
# Draw Figure 5 

setwd(store_doc)

png('figure5.png')
plot(x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,10), ylim = c(-0,10),asp = 1,
     main = "Fluoxetine-Treated Subjects - Clustering")

# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}

points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')

points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')

text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
dev.off()

######################################################
# Draw figure 6
png('figure6_control.png')
plot(asp = 1, x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,5), ylim = c(0,10),
     main = "Placebo Subjects - Clustering")
# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}
# draw the cluster center
points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')
points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')
text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
head(data_trans)
cluster4 = data_trans[data_trans$cluster == 1 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 1 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'red')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'red')
legend('topright', legend = c('Responsed','Non-responsed'), pch = c(20, 1), col = c('red','red'))
dev.off()

png('figure6_treat.png')
plot(asp = 1, x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,5), ylim = c(0,10),
     main = "Fluoxetine-Treated Subjects - Clustering")
# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}
# draw the cluster center
points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')
points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')
text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
cluster4 = data_trans[data_trans$cluster == 1 & data_trans$group == 2, ]
cluster4_r = data_trans[data_trans$cluster == 1 & data_trans$group == 2 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'purple')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'black')
legend('topright', legend = c('Responsed','Non-responsed'), pch = c(20, 1), col = c('purple','black'))
dev.off()

setwd(load_doc)