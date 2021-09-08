library(nlme)
library(mgcv)
#library(rhdf5) # for loading hdf4 
#library(rlist)
#library(rfast)# commented out
#library(plotfunctions)
#library(itsadug)
#library(ggplot2)
#library(writexl)
#library(itsadug)
#library(ramify)
# load utils:
source("utils.R")

##############################

run_analysis<-function(cone_type, description, df0,df1,se=2.58){
  # for comparing two experimental datasets 
  # se =  1.96 results in 95%CI and 2.58 in 99%CI
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  
  # normalize to mean
  df0<-normalize_to_mean(df0)
  df1<-normalize_to_mean(df1)
  
  # restructure 
  df0_re <- restructure_df(df0, df_led,0)
  df1_re <- restructure_df(df1, df_led,1)
  # join data 
  df = rbind(df0_re,df1_re)
  print(cone_type)
  
  # fit model
  r_gam<-gam(value~s(led, by=group, k=13) + group,
             data=df)
  
  # print details to .txt file
  sink("gam_summary.txt", append = TRUE)
  print(cone_type)
  print(description)
  print(summary(r_gam))
  sink() 
  
  # plot predictions 
  if(cone_type=='G'){
    ylims=c(-1.3,1)
  }else{
    ylims=c(-1,1)
  }
  plot_predictions_data(r_gam,cone_type, description, df0,df1,df_led,save_jpg=TRUE, ylims=ylims)
  plot_predictions_data(r_gam,cone_type, description, df0,df1,df_led,save_eps=TRUE, ylims=ylims)
  
  # save plot of differences
  out<-save_diff_plots(r_gam,cone_type,description, 
                       save_jpg = TRUE, 
                       save_eps = TRUE,
                       se=se) #options svg and jpg and eps
  
  # find differences 
  x <- find_difference(out$est, out$CI, f=1, xVals=out$led, as.vector = FALSE)
  
  return(x)
}

##############################################################################
# analyse hc_block vs control condition
###############################################################################
cone_types = c('R','G','B','U')

description = 'control_vs_hc_block' 
for (cone_type in cone_types){
  ## load all data from csv files
  folderpath='../data/all_cone_recordings/csvs/'
  
  # choose data
  if (description=='control_vs_hc_block') {
    # load control
    filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
    df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
    df0<-df0[complete.cases(df0),]
    # load other
    filename_1 = paste(folderpath, cone_type, '-Cone recordings - HCblock.csv', sep='' )
    df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
    df1<-df1[complete.cases(df1),]
    }
  
  x<-run_analysis(cone_type, description, df0,df1)
  
  # store to dataframe
  if(cone_type=='R'){
    df_results<-data.frame( 'cone'=cone_type, 
                         'description'=description, 
                         'sign_difference_start'=x$start,
                         'sign_difference_end'=x$end,
                         stringsAsFactors = FALSE)
  }else{
    for(j in c(1:length(x$start))){
      #print(paste('appendend index.. to results_df',j))
      newrow = c(cone_type, description,x$start[j],x$end[j])
      df_results<-rbind(df_results, newrow)
    }
  }
  print(paste('finished cone type',cone_type))
}

# write to excel file
#write_xlsx(df_results,"gam_results.xlsx")


###########################################################################
# "Other" additional experiments (Fig.2 d-h)
##########################################################################
folderpath='../data/all_cone_recordings/csvs/'

df_results<-data.frame( 'cone'='Dummy', 
                        'description'='Dummy', 
                        'sign_difference_start'=0,
                        'sign_difference_end'=0,
                        stringsAsFactors = FALSE)

### U low intensity ###
cone_type = 'U'

description = 'control_vs_low_intensity'
# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - low intensity.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]


x<- run_analysis(cone_type, description, df0, df1)

# store to dataframe
for(j in c(1:length(x$start))){
    #print(paste('appendend index.. to results_df',j))
    newrow = c(cone_type, description,x$start[j],x$end[j])
    df_results<-rbind(df_results, newrow)
}


### RED UV_ablation ###
cone_type = 'R'

# control
description = 'control_vs_uv_ablation'
# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - UVablation.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]

x<- run_analysis(cone_type, description, df0, df1)

# store to dataframe
for(j in c(1:length(x$start))){
  #print(paste('appendend index.. to results_df',j))
  newrow = c(cone_type, description,x$start[j],x$end[j])
  df_results<-rbind(df_results, newrow)
}


# HC block
description = 'hc_block_vs_uv_ablation'
# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - HCblock.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - UVablation and HCblock.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]

x<- run_analysis(cone_type, description, df0, df1)

# store to dataframe
for(j in c(1:length(x$start))){
  #print(paste('appendend index.. to results_df',j))
  newrow = c(cone_type, description,x$start[j],x$end[j])
  df_results<-rbind(df_results, newrow)
}


### Red low intensity ###
cone_type = 'R'

# control
description = 'control_vs_low_intensity'
# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - low intensity.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]

x<- run_analysis(cone_type, description, df0, df1)

# store to dataframe
for(j in c(1:length(x$start))){
  #print(paste('appendend index.. to results_df',j))
  newrow = c(cone_type, description,x$start[j],x$end[j])
  df_results<-rbind(df_results, newrow)
}


### BLUE UV ablated ###
cone_type = 'B'

description = 'control_vs_uv_ablation'
# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - UVablate.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]

x<- run_analysis(cone_type, description, df0, df1)

# store to dataframe
for(j in c(1:length(x$start))){
  #print(paste('appendend index.. to results_df',j))
  newrow = c(cone_type, description,x$start[j],x$end[j])
  df_results<-rbind(df_results, newrow)
}

# write to excel file
write_xlsx(df_results,"gam_results_suppl_exp_99CI.xlsx")


###################################################################################
# Check "oponency" for UV cones
###################################################################################

df_results_opponency<-data.frame( 'cone'='Dummy', 
                        'region'='Dummy', 
                        'sign_difference_start'=0,
                        'sign_difference_end'=0,
                        stringsAsFactors = FALSE)
cone_type = 'U'

plotting=TRUE
# only one is possible at a time:
save_jpg = FALSE
save_eps = TRUE

# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]

# load led
filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')

# normalize to mean
df0<-normalize_to_mean(df0)

# restructure 
df0_re <- restructure_df(df0, df_led,0, include_region=TRUE)

# fit model
r_gam<-gam(value~s(led,by=region,k=13)+region,
           data=df0_re)

# check the gam
#plot(r_gam)

# predict new data (newdata is dataframe with same predictors as initial df)
res = 0.01
lower = 360.*(1./res)
upper = 655.*(1./res)
x_new<- c(lower :upper )*res
regions = c('Nasal','Dorsal','SZ','Ventral')
for (region in regions){
  df_new = data.frame('led'=x_new, region=rep(region,dim(df_new)[1]))
  y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
  
  # look where difference is significant from 0
  x<-find_difference(y_fit_se$fit,y_fit_se$se.fit,  xVals=x_new, f=1.96)
  
  for(j in c(1:length(x$start))){
    newrow = c(cone_type, region,x$start[j],x$end[j])
    df_results_opponency<-rbind(df_results_opponency, newrow)
  }
  
  if(plotting){
    folderpath_plots = "../R_analysis/plots/"
    if(save_jpg){
      filename = paste(folderpath_plots,'jpgs/', "opponency",cone_type,'_',region,".jpeg", sep="")
      jpeg(file = filename, quality = 100)
    }
    if(save_eps){
      filename = paste(folderpath_plots,'eps/', "opponency",cone_type,'_',region,".eps", sep="")
      cairo_ps(file = filename)
    }
    
    ylims=c(-1.2,1.2)
    plot(x_new, y_fit_se$fit,type='l',ylim =ylims )
    lines(x_new, y_fit_se$fit+1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
    lines(x_new, y_fit_se$fit-1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)# compute zero crossing
    
    # plot raw data
    temp = df0[which(df0$region==region),]
    lines(df_led$Wavelength,colMeans(temp[1:13]), 
          type='p', col='blue',lty = 1,lwd = 2)
    # add lines 
    abline(h=0)
    abline(v=x$start, col='red',lwd=1, lty=2)
    abline(v=x$end, col='red',lwd=1, lty=2)
    segments(x0=x$start, y0=rep(ylims[1],length(x$start)), 
             x1=x$end, y1=rep(ylims[1],length(x$start)),
             col='red')
    # add title
    title(main= paste(region,': From 0 different'))
    dev.off()
    
  }
}

# write to excel file
write_xlsx(df_results,"gam_results_opponency_from0different.xlsx")


###################################################################################
# Opsin vs rest on log transormed opsins
#################################################################################

run_analysis_opsin<-function(cone_type, description, df0,df1){
  # df0: fitted opsin data
  # df1: experimental data
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  # normalize to mean df1
  df1<-normalize_to_mean(df1)
  # restructure df1
  df1_re <- restructure_df(df1, df_led,1)
  df<- df1_re
  # fit model
  r_gam<-gam(value~s(led, k=13) ,
             data=df)
  
  # predict with resolution of 1
  # predict new data (newdata is dataframe with same predictors as initial df)
  x_new<-c(360:655)
  groups <- rep(1,296)
  df_new <- data.frame('led'=x_new, 'group'=groups)
  y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
  
  # substract fitted_opsin from mean
  fitted_opsin <- df0[paste(cone_type,'_fitted_opsin', sep='')][,][df0$X>359 & df0$X<656]
  difference <- fitted_opsin - y_fit_se$fit
  # look where difference is significant from 0
  x<-find_difference(difference,y_fit_se$se.fit,  xVals=x_new,f=1.96)
  
  plot_opsin_difference(difference, x, y_fit_se, x_new,cone_type,description, save_jpg=TRUE, save_eps=FALSE)
  plot_opsin_difference(difference, x, y_fit_se, x_new,cone_type,description, save_jpg=FALSE, save_eps=TRUE)
  
  plot_prediction_opsin(fitted_opsin,y_fit_se,x_new,df_led,df1, save_jpg=TRUE, save_eps=FALSE )
  plot_prediction_opsin(fitted_opsin,y_fit_se,x_new,df_led,df1, save_jpg=FALSE, save_eps=TRUE )
  
  return(x)
}


#folderpath='../data/all_cone_recordings/csvs/'

# opsin vs control
for (cone_type in cone_types) {
  description = 'opsin_vs_control'
  # choose data
  if (description=='opsin_vs_control') {
    # load control
    filename_0 = paste(folderpath, 'opsin_fitted.csv', sep='')
    df0<-read.csv(filename_0, header=TRUE,sep=';',dec='.',stringsAsFactors=FALSE )
    #df0<-df0[360:656,]
    # load other
    filename_1 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
    df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
    df1<-df1[complete.cases(df1),]
  }
  
  # run analysis
  x<- run_analysis_opsin(cone_type, description, df0,df1)
  
  # store to dataframe
  #for(j in c(1:length(x$start))){
  #  #print(paste('appendend index.. to results_df',j))
  #  newrow = c(cone_type, description,x$start[j],x$end[j])
  #  df_results<-rbind(df_results, newrow)
  #}
}

# opsin vs block
for (cone_type in cone_types) {
  description = 'opsin_vs_hc_block'
  # choose data
  if (description=='opsin_vs_hc_block') {
    # load control
    filename_0 = paste(folderpath, 'opsin_fitted.csv', sep='')
    df0<-read.csv(filename_0, header=TRUE,sep=';',dec='.',stringsAsFactors=FALSE )
    #df0<-df0[360:656,]
    # load other
    filename_1 = paste(folderpath, cone_type, '-Cone recordings - HCblock.csv', sep='' )
    df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
    df1<-df1[complete.cases(df1),]
  }
  
  # run analysis
  x<- run_analysis_opsin(cone_type, description, df0,df1)
  
  # store to dataframe
  #for(j in c(1:length(x$start))){
  #  #print(paste('appendend index.. to results_df',j))
  #  newrow = c(cone_type, description,x$start[j],x$end[j])
  #  df_results<-rbind(df_results, newrow)
  #}
}
  
# write to excel file
write_xlsx(df_results,"gam_results.xlsx")


###############################################
# Opsin vs HC block v2
folderpath='../data/all_cone_recordings/csvs/'

run_analysis_opsin_2<-function(cone_type, description, df0,df1){
  # df0:  opsin data
  # df1: fitted experimental data (HC block)
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  # normalize to mean df1
  df1<-normalize_to_mean(df1)
  # restructure df1
  df1_re <- restructure_df(df1, df_led,1)
  df<- df1_re
  # fit model
  r_gam<-gam(value~s(led, k=13) ,
             data=df)
  plot(r_gam)
  
  # predict with resolution of 1
  # predict new data (newdata is dataframe with same predictors as initial df)
  x_new<-c(360:655)
  groups <- rep(1,296)
  df_new <- data.frame('led'=x_new, 'group'=groups)
  y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
  
  # substract raw opsin from mean
  opsin_raw <- df0[cone_type][,][df0$wavelength>359 & df0$wavelength<656]
  opsin<- opsin_raw/max(opsin_raw)
  difference <- opsin - y_fit_se$fit
  # look where difference is significant from 0
  x<-find_difference(difference,y_fit_se$se.fit,  xVals=x_new,f=1.96)
  
  plot_opsin_difference(difference, x, y_fit_se, x_new, cone_type,description, save_jpg=TRUE, save_eps=FALSE)
  plot_opsin_difference(difference, x, y_fit_se, x_new, cone_type,description, save_jpg=FALSE, save_eps=TRUE)
  
  plot_prediction_opsin(opsin, y_fit_se, x_new, save_jpg=TRUE, save_eps=FALSE )
  plot_prediction_opsin(opsin, y_fit_se, x_new, save_jpg=FALSE, save_eps=TRUE )
  
  return(x)
}


# opsin vs block
for (cone_type in cone_types) {
  description = 'opsin_vs_hc_block_v2'
  # choose data
  if (description=='opsin_vs_hc_block_v2') {
    # load opsin
    filename_0 = paste(folderpath, 'opsin wavelength.csv', sep='' )
    df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
    #load fitted opsin
    filename_1 = paste(folderpath, 'opsin_fitted_v2.csv', sep='')
    df1<-read.csv(filename_1, header=TRUE,sep=',',dec='.',stringsAsFactors=FALSE )
    df1 = subset(df1, select = -c(X) )
    df1 = df1[which(df1$cone_type == cone_type), ]
  }
  
  # run analysis
  x<- run_analysis_opsin_2(cone_type, description, df0,df1)
  
  # store to dataframe
  #for(j in c(1:length(x$start))){
  #  #print(paste('appendend index.. to results_df',j))
  #  newrow = c(cone_type, description,x$start[j],x$end[j])
  #  df_results<-rbind(df_results, newrow)
  #}
}


###############################################################
# one specific run

cone_type='R'

description = 'hc_block_vs_CNQX'

# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - HCblock.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]
# load other
filename_1 = paste(folderpath, cone_type, '-Cone recordings - CNQX.csv', sep='' )
df1<-read.csv(filename_1, header=TRUE,sep=';',dec=',')
df1<-df1[complete.cases(df1),]

## load all data from csv files
#folderpath='\\\\cin-storage/cschroeder/Documents/Biophysical Models/Cone_HC_interaction/data/all_cone_recordings/csvs/'
folderpath='../data/all_cone_recordings/csvs/'

# load led
filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')

# normalize to mean
df0<-normalize_to_mean(df0)
df1<-normalize_to_mean(df1)



# restructure 
df0_re <- restructure_df(df0, df_led,0)
df1_re <- restructure_df(df1, df_led,1)
# join data 
df = rbind(df0_re,df1_re)
# fit model
r_gam<-gam(value~s(led, by=group) + group,
           data=df)



plot(r_gam)

plot_predictions_data(r_gam,cone_type,description,df0,df1,df_led, save_jpg = FALSE)



out<-save_diff_plots(r_gam,cone_type,description,save_jpg = TRUE, save_svg = FALSE )


out<-plot_diff(r_gam,
               view='led',
               comp=list(group=c(0,1)))


# predict new data (newdata is dataframe with same predictors as initial df)
x_new<- append(c(360:650),c(360:650))
groups <- append(rep(0,291),rep(1,291))
df_new = data.frame('led'=x_new, 'group'=groups)
y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)

max(y_fit_se$fit[1:291])
max(y_fit_se$fit[292:582])

# plot raw mean predicitons +-se
plot(x_new[292:582], y_fit_se$fit[292:582],type='l',ylim =c(-1,1) )

################################################################################
# compute zero crossing and +-... for control G and B
##############################################################################
plotting=TRUE

# df to store results
df_zero_cross<-data.frame( 'cone'='Dummy', 
                           'zero_cross_0'=0,
                           'zero_cross_1'=0,
                           'zero_cross_2'=0,
                           'conf_interval_0'= 0,
                           'conf_interval_1' = 0,
                           'region'='Dummy',
                           stringsAsFactors = FALSE)

cone_types = c('G','B')
#cone_type = 'B'

for (cone_type in cone_types){
  ## load all data from csv files
  folderpath='../data/all_cone_recordings/csvs/'
  
  # choose data
  # load control
  filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
  df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
  df0<-df0[complete.cases(df0),]
  
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  
  # normalize to mean
  df0<-normalize_to_mean(df0)
  
  # restructure 
  df0_re <- restructure_df(df0, df_led,0)
  
  # fit model
  r_gam<-gam(value~s(led, k=13),
             data=df0_re)
  
  # check the gam
  plot(r_gam)
  
  # make predictions
  
  # predict new data (newdata is dataframe with same predictors as initial df)
  res = 0.01
  lower = 360.*(1./res)
  upper = 655.*(1./res)
  x_new<- c(lower :upper )*res
  df_new = data.frame('led'=x_new)
  y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
  
  # compute 0-crossing
  temp = y_fit_se$fit
  temp = temp[x_new<600]
  zero_cross_1 = x_new[which(temp==max(temp[temp<0]))] 
  
  # compute +- standard error
  temp = y_fit_se$fit-1.96*y_fit_se$se.fit
  temp = temp[x_new<600]
  zero_cross_0 = x_new[which(temp==max(temp[temp<0]))] 
  
  temp = y_fit_se$fit+1.96*y_fit_se$se.fit
  temp = temp[x_new<600]
  zero_cross_2 = x_new[which(temp==max(temp[temp<0]))] 
  
  if(plotting){
    ylims=c(-1.2,1.2)
    plot(x_new, y_fit_se$fit,type='l',ylim =ylims )
    lines(x_new, y_fit_se$fit+1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
    lines(x_new, y_fit_se$fit-1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)# compute zero crossing
    
    # plot raw data
    lines(df_led$Wavelength,colMeans(df0[1:13]), 
          type='p', col='blue',lty = 1,lwd = 2)
    abline(v=zero_cross_1,col='red')
    abline(h=0)
  }
  
  # add to df
  newrow = c(cone_type,
             zero_cross_0,
             zero_cross_1,
             zero_cross_2,
             zero_cross_1-zero_cross_0,
             zero_cross_2-zero_cross_1,
             'all')
  df_zero_cross<-rbind(df_zero_cross, newrow)
}

#### zone specific

for (cone_type in cone_types){
    
  ## load all data from csv files
  folderpath='../data/all_cone_recordings/csvs/'
  
  # choose data
  # load control
  filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
  df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
  df0<-df0[complete.cases(df0),]
  
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  
  # normalize to mean
  df0<-normalize_to_mean(df0)
  
  # restructure 
  df0_re <- restructure_df(df0, df_led,0, include_region=TRUE)
  
  # fit model
  r_gam<-gam(value~s(led,by=region,k=13)+region,
             data=df0_re)
  
  # check the gam
  #plot(r_gam)
  
  # make predictions
  
  # predict new data (newdata is dataframe with same predictors as initial df)
  res = 0.01
  lower = 360.*(1./res)
  upper = 655.*(1./res)
  x_new<- c(lower :upper )*res
  regions = c('Nasal','Dorsal','SZ','Ventral')
  for (region in regions){
    df_new = data.frame('led'=x_new, region=rep(region,dim(df_new)[1]))
    y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
    
    # compute 0-crossing
    temp = y_fit_se$fit
    temp = temp[x_new<550 & x_new>450]
    x_new_temp = x_new[x_new<550 & x_new>450]
    zero_cross_1 = x_new_temp[which(temp==max(temp[temp<0]))] 
    
    # compute +- standard error
    temp = y_fit_se$fit-1.96*y_fit_se$se.fit
    temp = temp[x_new<550 & x_new>450]
    x_new_temp = x_new[x_new<550 & x_new>450]
    zero_cross_0 = x_new_temp[which(temp==max(temp[temp<0]))] 
    
    
    temp = y_fit_se$fit+1.96*y_fit_se$se.fit
    temp = temp[x_new<550 & x_new>450]
    x_new_temp = x_new[x_new<550 & x_new>450]
    zero_cross_2 = x_new_temp[which(temp==max(temp[temp<0]))] 
    
    if(plotting){
      ylims=c(-1.2,1.2)
      plot(x_new, y_fit_se$fit,type='l',ylim =ylims )
      lines(x_new, y_fit_se$fit+1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
      lines(x_new, y_fit_se$fit-1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)# compute zero crossing
      
      # plot raw data
      temp = df0[which(df0$region==region),]
      lines(df_led$Wavelength,colMeans(temp[1:13]), 
            type='p', col='blue',lty = 1,lwd = 2)
      abline(v=zero_cross_0,col='red',lty = 2, lwd = 1)
      abline(v=zero_cross_1,col='red')
      abline(v=zero_cross_2,col='red',lty = 2, lwd = 1)
      
      abline(h=0)
      title(main= paste(cone_type, region,': zero-crossing'))
      
    }
    # add to df
    newrow = c(cone_type,
               zero_cross_0,
               zero_cross_1,
               zero_cross_2,
               zero_cross_1-zero_cross_0,
               zero_cross_2-zero_cross_1,
               region)
    df_zero_cross<-rbind(df_zero_cross, newrow)
  }
}


# write to excel file
write_xlsx(df_zero_cross,"zero_crossings_control_byzone.xlsx")

##############################################################################
# identify if zone specific differences in control
##############################################################################

#Functions
plot_predictions_zone_differences <- function(r_gam,cone_type, description,df0,df_led, save_jpg=FALSE,save_eps=FALSE, ylims=c(-1,1)){
  # only saves one format at a time!
  folderpath_plots = "../R_analysis/plots/regions/"
  ylims=c(-1,1)
  # predict new data (newdata is dataframe with same predictors as initial df)
  x_new<- append(append(append(c(360:655),c(360:655)),c(360:655)),c(360:655))
  regions <- append(append(append(rep(1,296),rep(2,296)),rep(3,296)),rep(4,296))
  retions <- factor(regions)
  df_new = data.frame('led'=x_new, 'region'=regions)
  y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)
  
  if(save_jpg){
    filename = paste(folderpath_plots,'jpgs/', "predictions_",cone_type,'_',description,".jpeg", sep="")
    jpeg(file = filename, quality = 100)
  }
  if(save_eps){
    filename = paste(folderpath_plots,'eps/', "predictions_",cone_type,'_',description,".eps", sep="")
    #setEPS()
    #postscript(file = filename)
    cairo_ps(file = filename)
  }
  
  # plot raw mean predicitons and 95% conf interval (assuming normal distr., calculated as 1.96*SE)
  colors = c('blue','red','green','violet') #colors are NOT automatically updated in the title!
  plot(x_new[1:296], y_fit_se$fit[1:296],type='l',ylim =c(-1,1),col=colors[1])
  lines(x_new[1:296], y_fit_se$fit[1:296]+1.96*y_fit_se$se.fit[1:296],type='l',lty = 2, lwd = 1,col=colors[1])
  lines(x_new[1:296], y_fit_se$fit[1:296]-1.96*y_fit_se$se.fit[1:296],type='l',lty = 2, lwd = 1,col=colors[1])
  
  lines(x_new[1:296], y_fit_se$fit[297:592],type='l',ylim =ylims ,col=colors[2])
  lines(x_new[1:296], y_fit_se$fit[297:592]+1.96*y_fit_se$se.fit[297:592],type='l',lty = 2, lwd = 1,col=colors[2])
  lines(x_new[1:296], y_fit_se$fit[297:592]-1.96*y_fit_se$se.fit[297:592],type='l',lty = 2, lwd = 1,col=colors[2])
  
  lines(x_new[1:296], y_fit_se$fit[593:888],type='l',ylim =ylims ,col=colors[3])
  lines(x_new[1:296], y_fit_se$fit[593:888]+1.96*y_fit_se$se.fit[593:888],type='l',lty = 2, lwd = 1,col=colors[3])
  lines(x_new[1:296], y_fit_se$fit[593:888]-1.96*y_fit_se$se.fit[593:888],type='l',lty = 2, lwd = 1,col=colors[3])
  
  lines(x_new[1:296], y_fit_se$fit[889:1184],type='l',ylim =ylims ,col=colors[4])
  lines(x_new[1:296], y_fit_se$fit[889:1184]+1.96*y_fit_se$se.fit[889:1184],type='l',lty = 2, lwd = 1,col=colors[4])
  lines(x_new[1:296], y_fit_se$fit[889:1184]-1.96*y_fit_se$se.fit[889:1184],type='l',lty = 2, lwd = 1,col=colors[4])
  
  # plot raw data
  lines(df_led$Wavelength,colMeans(df0[df0$region=='SZ',][1:13]),
        type='p', col=colors[3],lty = 1,lwd = 2)
  lines(df_led$Wavelength,colMeans(df0[df0$region=='Nasal',][1:13]),
        type='p', col=colors[2], lty = 1,lwd = 2)
  lines(df_led$Wavelength,colMeans(df0[df0$region=='Dorsal',][1:13]),
        type='p', col=colors[1], lty = 1,lwd = 2)
  lines(df_led$Wavelength,colMeans(df0[df0$region=='Ventral',][1:13]),
        type='p', col=colors[4], lty = 1,lwd = 2)
  
  # add title
  title(main= '1:blue:Dorsal , 2:red:Nasal, 3:green:AZ, 4:violet:Ventral')
  
  if(save_jpg | save_eps){
    dev.off()
  }
}

save_diff_plots_regions <- function(r_gam,cone_type,description,ylims=c(-0.4,0.8),groups,save_jpg=TRUE, save_svg=FALSE, save_eps=FALSE,se=1.96){
  # plot difference
  # se =  1.96 results in 95%CI and 2.58 in 99%CI
  # groups: groups to compare. for ex c(1,2)
  folderpath_plots = "../R_analysis/plots/regions/"
  if(save_svg){
    filename = paste(folderpath_plots, "differences_",cone_type,'_',description,'_',toString(groups),".svg", sep="")
    svg(file = filename)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(region=groups),
                   se=se,
                   ylim=ylims
    )
    dev.off()
  }
  if(save_jpg){
    filename = paste(folderpath_plots,'jpgs/', "differences_",cone_type,'_',description,'_',toString(groups),".jpeg", sep="")
    jpeg(file = filename, quality = 100)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(region=groups),
                   se=se,
                   ylim=ylims
    )
    dev.off()
  }
  if(save_eps){
    filename = paste(folderpath_plots,'eps/',"differences_",cone_type,'_',description,'_',toString(groups),".eps", sep="")
    #setEPS()
    #postscript(file = filename)
    cairo_ps(file = filename)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(region=groups),
                   se=se,
                   ylim=ylims
    )
    dev.off()
  }
  return(out)
}

normalize_to_mean_regionwise <- function(df){
  # normalize to have mean of 1 per region
  for (item in unique(df$region)){
    df[df$region==item, 1:13] <- df[df$region==item,1:13] / max(abs(colMeans(df[df$region==item,1:13])))
  }
  return(df)
}

## pipeline
cone_types = c('R','G','B','U')
description = 'zone_differences'
count=1
for (cone_type in cone_types){
  #cone_type='R'  
  
  ## load all data from csv files
  folderpath='../data/all_cone_recordings/csvs/'
  
  # choose data
  # load control
  filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
  df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
  df0<-df0[complete.cases(df0),]
  
  # load led
  filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
  df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')
  
  # normalize to mean
  # no normalization at the moment!!!
  #df0<-normalize_to_mean_regionwise(df0)
  
  # restructure 
  df0_re <- restructure_df(df0, df_led,0, include_region=TRUE)
  
  # fit model ignoring region
  #r_gam1<-gam(value~s(led,k=13),
  #            data=df0_re)
  
  # fit model by region
  r_gam<-gam(value~s(led,by=region,k=13)+region,
             data=df0_re)
  
  # statistical test if there are differences at all / which model is the better one
  #anova(r_gam1,r_gam,test="Chisq")
  #BIC(r_gam1)
  #BIC(r_gam)
  
  
  # print details to .txt file
  sink("gam_summary_regions.txt", append = TRUE)
  print(cone_type)
  print(description)
  print(summary(r_gam))
  print('ANOVA on GAM: anova.gam results:')
  print(anova.gam(r_gam))
  print('------')
  sink() 
  print(anova.gam(r_gam))
  
  ## plot predictions 
  
  # set ylims
  if(cone_type=='G'){
    ylims=c(-1.3,1)
  }else{
    ylims=c(-1,1)
  }
  
  # plot predictions and save them 
  plot_predictions_zone_differences(r_gam,cone_type, description, df0,df_led,save_jpg=TRUE, ylims=ylims)
  
  
  # iterate over all possible comparisons
  possible_combinations <- list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  for (groups in possible_combinations  ){
    out<-save_diff_plots_regions(r_gam,cone_type,description,groups=groups, 
                                 save_jpg = TRUE, 
                                 save_eps = FALSE,
                                 se=2.58,
                                 ylims=c(-0.6,0.8)) #options svg and jpg and eps
    
    # find differences 
    x <- find_difference(out$est, out$CI, f=1, xVals=out$led, as.vector = FALSE)
    
    # store to dataframe
    print(groups)
    if(count==1){
      df_results<-data.frame( 'cone'=cone_type,
                              'groups'= toString(groups),
                              'description'=description, 
                              'sign_difference_start'=x$start,
                              'sign_difference_end'=x$end,
                              stringsAsFactors = FALSE)
    }else{
      if (is.null(x)){
        newrow = c(cone_type, toString(groups), description, 'not sign.', 'not sign.')
        df_results<-rbind(df_results, newrow)
      }else{
        for(j in c(1:length(x$start))){
          #print(paste('appendend index.. to results_df',j))
          newrow = c(cone_type, toString(groups), description, x$start[j], x$end[j])
          df_results<-rbind(df_results, newrow)
        }
      }
    }
    count=count+1
    }
  print(paste('finished cone type',cone_type))
  
}


# write to excel file
write_xlsx(df_results,"gam_results_regions.xlsx")

############################################################
# some tests
cone_type='R'  

## load all data from csv files
folderpath='../data/all_cone_recordings/csvs/'

# choose data
# load control
filename_0 = paste(folderpath, cone_type, '-Cone recordings - control_merged.csv', sep='' )
df0<-read.csv(filename_0, header=TRUE,sep=';',dec=',')
df0<-df0[complete.cases(df0),]

# load led
filename_led = paste(folderpath, 'LED wavelength.csv', sep='' )
df_led <-read.csv(filename_led, header=TRUE,sep=';',dec=',')

# normalize to mean
# no normalization at the moment!!!
#df0<-normalize_to_mean_regionwise(df0)

# restructure 
df0_re <- restructure_df(df0, df_led,0, include_region=TRUE)

# fit model ignoring region
#r_gam1<-gam(value~s(led,k=13),
#            data=df0_re)

# fit model by region
r_gam<-gam(value~s(led,by=region,k=13)+region,
           data=df0_re)

# statistical test if there are differences at all / which model is the better one
#anova(r_gam1,r_gam,test="Chisq")
#BIC(r_gam1)
#BIC(r_gam)

# predict new data (newdata is dataframe with same predictors as initial df)
x_new<- append(append(append(c(360:655),c(360:655)),c(360:655)),c(360:655))
regions <- append(append(append(rep('SZ',296),rep('Nasal',296)),rep('Dorsal',296)),rep('Ventral',296))
regions <- factor(regions)
df_new = data.frame('led'=x_new, 'region'=regions)
y_fit_se <- predict(r_gam, type="response", newdata=df_new, se.fit = TRUE)




#######################################################
# random stuff

cairo_ps(file = "test.eps", onefile = FALSE, fallback_resolution = 600)
plot(rnorm(100), main="Hey Some Data")
dev.off()
# print to file
sink("test.txt", append = TRUE)
print(temp)
sink() 

# clear workspace: rm(list = ls()) 