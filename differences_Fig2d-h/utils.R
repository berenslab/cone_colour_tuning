library(nlme)
library(mgcv)
#library(rhdf5) # for loading hdf4 
#library(rlist)
#library(rfast) # commented out
#library(plotfunctions)
library(ggplot2)
library(writexl)
library(itsadug)

# restructure dfs function
restructure_df <- function(df_data, df_led, group_id, include_region=FALSE) {
  # prepares the df for the GAM function 
  # df_data: dataframe of the cone data
  # df_led: dataframe containig wavelenghths for LEDs
  # group_id: specify group id (integer) for later comparison
  
  # flatten df
  y0 <- vector()
  x0 <- vector()
  ids0 <- vector()
  regions <- vector()
  for (name in colnames(df_data)){
    if (grepl('LED',name)){
      y_temp<-df_data[,name]
      y0<-append(y0,y_temp)
      x_temp<- df_led[which(df_led$Stimulus == name), ]$Wavelength
      x0<-append(x0, array(data=x_temp,dim=length(y_temp)))
      # ids (cone ids)
      ids0<-append(ids0, c(1:length(y_temp))  )       
    }
    if(include_region){
      if (grepl('region',name)){
        regions_temp<-df_data[,name]
        regions<-append(regions,regions_temp)
      }
    }
    }
  
  # group (control vs test)
  group0 = factor(rep(group_id, each=length(y0) ))
  # construct new df
  df_new = data.frame('led'=x0 ,'value'=as.vector(y0), 'group'=group0, 'id'=ids0 )
  
  if(include_region){
    regions<-factor(regions)
    df_new = data.frame('led'=x0 ,'value'=as.vector(y0), 'group'=group0, 'id'=ids0,'region'=regions )
  }
  
  return(df_new)
}

normalize_to_mean <- function(df){
  # normalize to have mean of 1
  df[1:13]<-df[1:13] / max(colMeans(df[1:13]))
  return(df)
}

save_diff_plots <- function(r_gam,cone_type,description,save_jpg=TRUE, save_svg=FALSE, save_eps=FALSE,se=1.96){
  # plot difference
  # se =  1.96 results in 95%CI and 2.58 in 99%CI
  folderpath_plots = "../R_analysis/plots/"
  if(save_svg){
    filename = paste(folderpath_plots, "differences_",cone_type,'_',description,".svg", sep="")
    svg(file = filename)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(group=c(0,1)),
                   se=se
    )
    dev.off()
  }
  if(save_jpg){
    filename = paste(folderpath_plots,'jpgs/', "differences_",cone_type,'_',description,".jpeg", sep="")
    jpeg(file = filename, quality = 100)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(group=c(0,1)),
                   se=se
    )
    dev.off()
  }
  if(save_eps){
    filename = paste(folderpath_plots,'eps/',"differences_",cone_type,'_',description,".eps", sep="")
    #setEPS()
    #postscript(file = filename)
    cairo_ps(file = filename)
    out<-plot_diff(r_gam,
                   view='led',
                   comp=list(group=c(0,1)),
                   se=se
    )
    dev.off()
  }
  return(out)
}

plot_predictions_data <- function(r_gam,cone_type, description,df0,df1,df_led, save_jpg=FALSE,save_eps=FALSE, ylims=c(-1,1)){
  # only saves one format at a time!
  folderpath_plots = "../R_analysis/plots/"
  
  # predict new data (newdata is dataframe with same predictors as initial df)
  x_new<- append(c(360:655),c(360:655))
  groups <- append(rep(0,296),rep(1,296))
  df_new = data.frame('led'=x_new, 'group'=groups)
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
  plot(x_new[1:296], y_fit_se$fit[297:592],type='l',ylim =ylims )
  lines(x_new[1:296], y_fit_se$fit[297:592]+1.96*y_fit_se$se.fit[297:592],type='l',lty = 2, lwd = 1)
  lines(x_new[1:296], y_fit_se$fit[297:592]-1.96*y_fit_se$se.fit[297:592],type='l',lty = 2, lwd = 1)
  
  lines(x_new[1:296], y_fit_se$fit[1:296],type='l',ylim =c(-1,1) )
  lines(x_new[1:296], y_fit_se$fit[1:296]+1.96*y_fit_se$se.fit[1:296],type='l',lty = 2, lwd = 1)
  lines(x_new[1:296], y_fit_se$fit[1:296]-1.96*y_fit_se$se.fit[1:296],type='l',lty = 2, lwd = 1)
  
  # plot raw data
  lines(df_led$Wavelength,colMeans(df0[1:13]), 
        type='p', col='blue',lty = 1,lwd = 2)
  lines(df_led$Wavelength,colMeans(df1[1:13]), 
        type='p', col='red', lty = 1,lwd = 2)
  
  title(main= 'blue:0 , red:1')
  
  dev.off()
}


## for opsin analysis:

plot_opsin_difference<-function(difference, x, y_fit_se, x_new,cone_type,description, save_jpg=TRUE, save_eps=FALSE){
  # plot prediction
  folderpath_plots = "../R_analysis/plots/"
  if(save_jpg){
    filename = paste(folderpath_plots,'jpgs/', "difference_",cone_type,'_',description,".jpeg", sep="")
    jpeg(file = filename, quality = 100)
  }
  if(save_eps){
    filename = paste(folderpath_plots,'eps/', "difference_",cone_type,'_',description,".eps", sep="")
    #setEPS()
    #postscript(file = filename)
    cairo_ps(file = filename)
  }
  
  # plot predicited difference and 95% conf interval (assuming normal distr., calculated as 1.96*SE)
  ylims=c(min(difference)-0.1, max(difference)+0.1)
  plot(x_new, difference,type='l',ylim =ylims )
  lines(x_new, difference+1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
  lines(x_new, difference-1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
  lines(x_new, rep(0, length(x_new)) )
  # add lines 
  abline(v=x$start, col='red',lwd=1, lty=2)
  abline(v=x$end, col='red',lwd=1, lty=2)
  segments(x0=x$start, y0=rep(ylims[1],length(x$start)), 
           x1=x$end, y1=rep(ylims[1],length(x$start)),
           col='red')
  
  title(main= 'Difference: opsin - other')
  
  dev.off()
}


plot_prediction_opsin<-function(fitted_opsin,y_fit_se,x_new,df_led,df0, save_jpg=TRUE, save_eps=FALSE ){
  # plot prediction
  folderpath_plots = "../R_analysis/plots/"
  if(save_jpg){
    filename = paste(folderpath_plots,'jpgs/', "predictions_",cone_type,'_',description,".jpeg", sep="")
    jpeg(file = filename, quality = 100)
  }
  if(save_eps){
    filename = paste(folderpath_plots,'eps/', "difference_",cone_type,'_',description,".eps", sep="")
    #setEPS()
    #postscript(file = filename)
    cairo_ps(file = filename)
  }
  # plot prediction and opsin and 95% conf interval (assuming normal distr., calculated as 1.96*SE)
  ylims=c(-1,1)
  plot(x_new, y_fit_se$fit,type='l',ylim =ylims )
  lines(x_new, y_fit_se$fit+1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
  lines(x_new, y_fit_se$fit-1.96*y_fit_se$se.fit,type='l',lty = 2, lwd = 1)
  lines(x_new,fitted_opsin,col='red' )
  title(main= 'red: opsin - black: other')
  # plot raw data
  lines(df_led$Wavelength,colMeans(df0[1:13]), 
        type='p', col='blue',lty = 1,lwd = 2)
  
  dev.off()
  
}
