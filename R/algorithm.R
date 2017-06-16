#' @title newfunction
#' @description newfunction
#' This function allows you to express your love of cats.
#' @param input json format.
#' @export
#' @import jsonlite
#' cat_function()


algorithm <- function(input){
  options(warn=-1)
  # dat<-fromJSON(input)
  dat<-input
  dat[,1]<-as.Date(as.character(dat[,1]), format='%Y-%m-%d')
  dat<-dat[order(dat[,1]),]

  t_c<- timeseries_conversion(dat[,2], order.by=dat[,1])
  tsm<-timeseries(t_c)

  indx <- 1:length(dat[,1])
  var.use <- dat[,2]
  lo <- loess(var.use~indx)
  predicted<-predict(lo)


  temp<-var.use

  analysis<-outlier_detection(tsm,maxit.iloop = 1)
  if(nrow(analysis)==0){
    diff<-abs(predicted-var.use)
    per_diff<-diff/predicted *100
    per_limt<-min(per_diff[analysis$ind])
    limit<-max(per_diff[- analysis$ind])
    limit<-mean(limit,per_limt)
    limit_s<-summary(per_diff)

  } else {
    diff<-abs(predicted-var.use)
    per_diff<-diff/predicted *100
    limit<-90
    limit_s<-summary(per_diff)

  }

  n<-length(temp)
  for(i in analysis$ind){
    if(n<=6){
      if(i<4){
        rng<-(i+1):(n)
      } else rng<-c((i-2):(i-1),(n-1):(n))

    }else{rng<-c((i-3):(i-1),(i+1):(i+3))}



    temp[i]<-mean(temp[rng])
  }



  temp<-predict(loess(temp~indx))


  outlier<-list()
  outlier[1:length(dat[,1])]<-"0"
  outlier[analysis$ind]<-"outlier"
  tmp<-data.frame(dat[,1],predicted,var.use,temp +(temp*limit)/100 ,temp -(temp*limit)/100,temp +(temp*limit_s[5]*2)/100 ,temp -(temp*limit_s[5]*2)/100,unlist(outlier))
  colnames(tmp)<-c("Date","Predicted","Actual","upperOffOutlierLimit","lowerOffOutlierLimit","upperSafe3rdQuantLimit","lowerSafe3rdQuantLimit","Sensitivity")

  chk<-tmp[tmp$Sensitivity=="outlier",]
  if(nrow(chk)!=0){
    for(i in nrow(chk)){
      if(abs(chk$Actual[i]-chk$upperOffOutlierLimit[i]) <abs(chk$Actual[i]-chk$lowerOffOutlierLimit[i])){

        if(!(chk$Actual[i]>chk$upperOffOutlierLimit[i]))
        {
          chk$upperOffOutlierLimit[i]<- chk$Actual[i]*0.30
        }
        if(!(chk$Actual[i]>chk$upperSafe3rdQuantLimit[i]))
        {
          chk$upperSafe3rdQuantLimit[i]<- chk$Actual[i]*0.30
        }
      }else{
        if((chk$Actual[i]>chk$lowerOffOutlierLimit[i]))
        {
          chk$lowerOffOutlierLimit[i]<- chk$Actual[i]*1.2
        }
        if((chk$Actual[i]>chk$upperSafe3rdQuantLimit[i]))
        {
          chk$lowerSafe3rdQuantLimit[i]<- chk$Actual[i]*1.2
        }
      }

    }
    chk->tmp[tmp$Sensitivity=="outlier",]}
  ###############################################
  options(warn=0)
  # if(plot_value==T){
  #   plot_dia<-ggplot(tmp,aes(Date))+
  #     geom_line(data = tmp,aes(Date,Predicted),size=2,col="steelblue")+
  #     geom_point(data = tmp,aes(Date,Actual,col=Sensitivity))+
  #     geom_ribbon(aes(ymin=lowerOffOutlierLimit,ymax=upperOffOutlierLimit),alpha=.25,fill="steelblue")+
  #     labs(title='Outlier Detection')+
  #     theme(axis.text= element_text(face = "bold",size = 7.5)     ,
  #           axis.title.y =element_blank() ,
  #           axis.title.x =element_blank() ,
  #           plot.title = element_text(size =11, face = "bold") ,
  #           plot.subtitle = element_text(size = 12),legend.position = "None"
  #     )+
  #     scale_color_manual(values=c("0" = "#333BFF", "outlier" = "#CC6600"))
  #
  #   return(plot_dia)
  # } else return(tmp)

  json<-toJSON(tmp)
  json
}


