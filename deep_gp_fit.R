library(deepgp)

#lims data
data<-read.csv("22444089_half.csv")
x<-data$norm_time
y<-data$norm_vaf
data$signal<-y
x_new <- seq(0, max(data$months_elapsed), length = 100)
fit<-fit_one_layer(x,y,cov="exp2",verb=FALSE,nmcmc = 10000,g_0 = 0.1,theta_0 = 0.5,
                   settings = list(l = 1, u = 10, alpha = list(g = 1.1, theta = 1.),
                                   beta = list(g = 3.9, theta = 0.1)))
fit<-fit_two_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 20000,g_0 = 0.01,theta_w_0 = 0.1,theta_y_0 = 0.5,true_g = 0.06,
                   settings = list(l = 1, u = 5, alpha = list(g = 1.5, theta_w = 1.2, theta_y = 1.1),
                                   beta = list(g = 3.9, theta_w = 2.95 , theta_y = 0.175 )))
fit<-fit_three_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 10000,
                     g_0 = 0.01,theta_w_0 = 0.5,theta_y_0 = 0.1,
                     settings = list(l = 1, u = 2, alpha = list(g = 1.5, theta_w = 10.5, theta_z = 5.5,theta_y = 1.5),
                                     beta = list(g = 3.9, theta_w = 20., theta_z = 10.,theta_y = 3.9/6)))
fit<-trim(fit,burn=1000)
fit <- trim(fit, 2000)
plot(fit)
plot(fit,hidden = TRUE)
fit <- predict(fit, x_new, store_latent = TRUE)
plot_predict(fit,data)
plot_abs_w(fit,data)
hist(fit$theta)
hist(fit$tau2)
hist(fit$g)
hist(fit$theta_w)
hist(fit$theta_w, breaks=c(0.,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6,max(fit$theta_w)),xlim=c(0,0.6))
noise<-fit$tau2*fit$g
hist(noise)
mean(noise)
c(mean(fit$tau2),mean(fit$g),mean(fit$theta_w),mean(fit$theta_y))

#simulations
setwd("~/gp_process_model/simulations")
sim_data<-read.csv("piecewise_liner_w_noise_biorhythm.csv")
sim_data<-read.csv("liner_w_noise_biorhythm.csv")
sim_data<-read.csv("liner_w_noise_biorhythm_50percent.csv")
sim_data<-read.csv("liner_w_noise_biorhythm_30percent.csv")
x<-sim_data$norm_time
y<-sim_data$norm_response
fit<-fit_one_layer(x,y,cov="exp2",verb=FALSE,nmcmc = 10000,g_0 = 0.1,theta_0 = 0.5,
                   settings = list(l = 1, u = 10, alpha = list(g = 1.1, theta = 1.1),
                                   beta = list(g = 1.9, theta = 3.9)))
fit<-fit_two_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 20000,g_0 = 0.1,theta_w_0 = 0.1,theta_y_0 = 3.0,true_g = NULL,
                   settings = list(l = 1, u = 5, alpha = list(g = 1.1, theta_w = 1.1, theta_y = 6.5),
                                   beta = list(g = 3.9, theta_w = 6.95 , theta_y = 0.6 )))
fit<-fit_three_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 10000,
                     g_0 = 0.01,theta_w_0 = 0.5,theta_y_0 = 0.1,
                     settings = list(l = 1, u = 2, alpha = list(g = 1.5, theta_w = 10.5, theta_z = 5.5,theta_y = 1.5),
                                     beta = list(g = 3.9, theta_w = 20., theta_z = 10.,theta_y = 3.9/6)))
plot(fit,hidden = TRUE)
plot(fit)
fit <- trim(fit, 2000, 2)
fit <- trim(fit, 2000)
x_new <- seq(0, 1.0, length = 100)
fit <- predict(fit, x_new, store_latent = TRUE)
hist(fit$theta_w)
hist(fit$theta_w, breaks=c(0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.,max(fit$theta_w)),xlim=c(0,3.))
hist(fit$theta_y)
hist(fit$theta_y,breaks=c(0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.,max(fit$theta_y)),xlim=c(0,3.))
hist(fit$theta)
hist(fit$tau2)
hist(fit$tau2, breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,max(fit$tau2)),xlim=c(0,13))
sig_variation <- sqrt(1./fit$tau2)
hist(sig_variation)
mean(sig_variation)
hist(fit$g)
hist(fit$g,breaks=c(0.,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.3,max(fit$g)),xlim=c(0,0.3))
noise<-fit$tau2*fit$g
noise<-sqrt(noise)
hist(noise,breaks=c(0.,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,max(noise)),xlim=c(0.1,0.5))
mean(noise)
c(mean(fit$tau2),mean(fit$g),mean(fit$theta_w),mean(fit$theta_y))

plot_predict(fit,sim_data)
plot_abs_w(fit,sim_data)
plot_abs_z(fit,sim_data)
mean(fit$tau2)
mean(fit$g)

#signatera noise
sig_data<-read.csv("595939.csv") #noise
sig_data<-read.csv("645945.csv") #noise
sig_data<-read.csv("619087.csv") #signal
sig_data<-read.csv("508582.csv") #signal
x_new <- seq(0, 1, length = 100)
x<-sig_data$norm_time
y<-sig_data$log_mtm
y<-sig_data$log_mean_VAF
sig_data$signal<-y
fit<-fit_two_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 50000,settings = list(l = 1, u = 3, alpha = list(g = 1.1, theta_w = 1.5, theta_y = 1.5),beta = list(g = 6.9, theta_w = 12.9, theta_y = 8.9)))
fit<-fit_two_layer(x,y,cov="exp2",verb=FALSE, nmcmc = 50000,settings = list(l = 1, u = 3, true_g = 0.003,alpha = list(g = 1.1, theta_w = 1.5, theta_y = 1.5),beta = list(g = 3.9, theta_w = 3.9, theta_y = 3.9)))
fit <- trim(fit, 4000)
plot(fit)
plot(fit,hidden = TRUE)
plot_abs_w(fit,sig_data)
fit <- predict(fit, x_new, store_latent = TRUE)
plot_predict(fit,sig_data)

plot_abs_w <- function(x,data){
  num_plots = 1000
  indx <- floor(seq(from = 1, to = x$nmcmc, length = num_plots))
  if (indx[1] == 0) 
    indx[1] = 1
  col <- heat.colors(num_plots + 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
  o <- order(x$x)
  scale_max = log(max(abs(unlist(x$w[indx]))))
  scale_min = (min(abs(unlist(x$w[indx]))))
  if(scale_min > 0){
    scale_min = log(scale_min)
  }else{
    scale_min  = log(0.00001)
  }
  #plot(x$x[o], log(abs(x$w[[indx[1]]][o] - mean(unlist(x$w[[indx[1]]])))), 
  #     type = "l", col = col[1],  ylim = c(scale_min, scale_max), xlab = "X", ylab = "W")
  #mw <- matrix(log(abs(unlist(x$w[[indx[1]]][o]) - mean(unlist(x$w[[indx[1]]])))),nrow=1,ncol=nrow(x$x))
  mw <- matrix((abs(unlist(x$w[[indx[1]]][o]) - mean(unlist(x$w[[indx[1]]])))),nrow=1,ncol=nrow(x$x))
  for (j in 2:length(indx)) {
    #lines(x$x[o], log(abs(x$w[[indx[j]]][o] - mean(unlist(x$w[[indx[j]]])))), col = col[j])
    #mw<-rbind(mw,c(log(abs(unlist(x$w[[indx[j]]][o]) - mean(unlist(x$w[[indx[j]]]))))))
    mw<-rbind(mw,c((abs(unlist(x$w[[indx[j]]][o]) - mean(unlist(x$w[[indx[j]]]))))))
  }
  mw<-colMeans(mw)
  scale_max = max(mw)
  scale_min = min(mw)
  plot(x$x[o],mw,type = "l",  ylim = c(scale_min, scale_max), xlab = "X", ylab = "W",col='red')
  #lines(x$x[o],(data$signal[o] - min(data$signal))/(max(data$signal) - min(data$signal))*(scale_max - scale_min) + scale_min,col='black')
}

plot_abs_z <- function(x,data){
  num_plots = 1000
  indx <- floor(seq(from = 1, to = x$nmcmc, length = num_plots))
  if (indx[1] == 0) 
    indx[1] = 1
  col <- heat.colors(num_plots + 100)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
  o <- order(x$x)
  scale_max = log(max(abs(unlist(x$z[indx]))))
  scale_min = (min(abs(unlist(x$z[indx]))))
  if(scale_min > 0){
    scale_min = log(scale_min)
  }else{
    scale_min  = log(0.00001)
  }
  plot(x$x[o], log(abs(x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]))), 
       type = "l", col = col[1],  ylim = c(scale_min, scale_max))
  for (j in 2:length(indx)) {
    lines(x$x[o], log(abs(x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]))), 
          col = col[j])
  }
  lines(x$x[o],(data$signal[o] - min(data$signal))/(max(data$signal) - min(data$signal))*(scale_max - scale_min) + scale_min,col='black')
}

plot_predict <- function(x,data){
  par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
  q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
  q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
  o <- order(x$x_new)
  plot(x$x_new[o], x$mean[o], type = "l", xlab = "X", 
       ylab = "Y", ylim = c(min(q1), max(q3)), col = "blue", 
       main = "Posterior Mean and 95% PI")
  lines(x$x_new[o], q1[o], col = "blue", lty = 2)
  lines(x$x_new[o], q3[o], col = "blue", lty = 2)
  points(x$x, x$y, pch = 20)
  lines(x$x,data$norm_signal,col='black')
}

#rank plots
theta_w_111 <-fit$theta_w
df_01<-data.frame(data=theta_w_01)
df_01$label = "theta_w_01"
df<-rbind(df_01,df_02,df_03,df_04,df_05)
sorted_df<-df[order(df$data),]
sorted_df$rank<-1:nrow(sorted_df)
tmp<-sorted_df[which(label=='theta_w_01')]
hist(tmp$rank)

