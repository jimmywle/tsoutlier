#' @title Algorithm
#' @description Function deals with outlier detection functionality for timeseries data
#'  USE Some Prebild function from packages
#' @param input json format.
#' @export
#' @import jsonlite
#' algorithm()


algorithm <- function(input){
  options(warn=-1)
  # dat<-fromJSON(input)
  input[,2]<-as.numeric(input[,2])
  dat<-input

  if(nchar(dat[,1][1])>10){
    seq_len<-nrow(dat)
    input[,1]<-seq(as.Date("2000/1/1"), by = "day", length.out = seq_len)
    input[,1]<-as.Date(as.character(input[,1]), format='%Y-%m-%d')
    dat[,1]<-as.POSIXct(dat[,1])
  }else {
    input[,1]<-as.Date(as.character(input[,1]), format='%Y-%m-%d')
    input<-input[order(input[,1]),]

    dat[,1]<-as.Date(as.character(dat[,1]), format='%Y-%m-%d')
    dat<-dat[order(dat[,1]),]
  }
  t_c<- timeseries_conversion(input[,2], order.by=input[,1])
  tsm<-timeseries(t_c)

  indx <- 1:length(dat[,1])
  var.use <- dat[,2]




  temp<-var.use

  #analysis<-outlier_detection(tsm,maxit.iloop = 1)
  result<-tso(tsm)
  analysis<-result$outliers
  lo <- loess(result$yadj~indx)
  predicted<-predict(lo)
####
  if(nrow(analysis)!=0){
  dat[,2][result$times]<-result$yadj[result$times]
  outlier<-list()
  outlier[1:length(dat[,1])]<-0
  outlier[analysis$ind]<-1

  raw_limit<-sum(abs(var.use-dat[,2])/var.use)/length(result$times)
  limit<-raw_limit*0.5
  chk<-unique((var.use-dat[,2])/var.use)[-1]

  if(sum(chk<limit)==0){
    limit_s<-limit*0.5
  } else{
    limit<-limit*0.75
    limit_s<-limit*0.5
  }

  if(sum(chk<limit)!=0){
    limit<-limit*0.75
    limit_s<-limit*0.75
  }

  temp<-dat[,2]
  tmp<-data.frame(dat[,1],predicted,var.use,temp +(temp*limit) ,temp -(temp*limit),temp +(temp*limit_s),temp -(temp*limit_s),unlist(outlier))

  } else{

    temp<-dat[,2]
    outlier<-list()
    outlier[1:length(dat[,1])]<-0

    limit<-0.5
    limit_s<-limit*0.6

    tmp<-data.frame(dat[,1],predicted,var.use,temp +(temp*limit) ,temp -(temp*limit),temp +(temp*limit_s),temp -(temp*limit_s),unlist(outlier))

  }

    colnames(tmp)<-c("Date","Predicted","Actual","upperOffOutlierLimit","lowerOffOutlierLimit","upperSafe3rdQuantLimit","lowerSafe3rdQuantLimit","Sensitivity")

  # chk<-tmp[tmp$Sensitivity==1,]
  # if(nrow(chk)!=0){
  #   for(i in nrow(chk)){
  #     if(abs(chk$Actual[i]-chk$upperOffOutlierLimit[i]) <abs(chk$Actual[i]-chk$lowerOffOutlierLimit[i])){
  #
  #       if(!(chk$Actual[i]>chk$upperOffOutlierLimit[i]))
  #       {
  #         chk$upperOffOutlierLimit[i]<- chk$Actual[i]*0.30
  #       }
  #       if(!(chk$Actual[i]>chk$upperSafe3rdQuantLimit[i]))
  #       {
  #         chk$upperSafe3rdQuantLimit[i]<- chk$Actual[i]*0.30
  #       }
  #     }else{
  #       if((chk$Actual[i]>chk$lowerOffOutlierLimit[i]))
  #       {
  #         chk$lowerOffOutlierLimit[i]<- chk$Actual[i]*1.2
  #       }
  #       if((chk$Actual[i]>chk$upperSafe3rdQuantLimit[i]))
  #       {
  #         chk$lowerSafe3rdQuantLimit[i]<- chk$Actual[i]*1.2
  #       }
  #     }
  #
  #   }
  #   chk->tmp[tmp$Sensitivity==1,]}
  # ###############################################
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

  # json<-toJSON(tmp)
  # json
  tmp
}


