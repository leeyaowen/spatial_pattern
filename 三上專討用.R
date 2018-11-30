library(DBI)
library(RPostgreSQL)
library(ecespa)

library(sqldf)
library(dplyr) #過濾資料用
library(spatstat)
library(magrittr)
library(ggplot2)
library(ggthemes)
drv<-dbDriver("PostgreSQL") #宣告資料庫名稱
con<-dbConnect(drv,user="postgres",password="2717484",dbname="LJC") #建立連線字串
sp15<-dbGetQuery(con,"select sp from ljc2018 group by sp having count(sp)>=15") #產生株數大於15的名錄
sp15$sp %<>% iconv(.,"UTF-8","CP950")
splist<-as.matrix(sp15) #物種清單

db<-read.csv("./ljc2018.csv")
sp15<-sqldf("select sp from db group by sp having count(sp)>=15")
splist<-as.matrix(sp15)

#ggplot
dt<-sqldf("select * from db")
dttemp<-filter(dt,sp=="南仁灰木")
mypattern<-ppp(dttemp$x4,dttemp$y4,c(0,300),c(0,100))
L<-envelope(mypattern,Lest,nsim = 99,correction="best")
Ltemp<-mutate(L,obs=obs-r,theo=theo-r,lo=lo-r,hi=hi-r)
Ltemp<-as.data.frame(Ltemp)
#Lplot<-ggplot(Ltemp,aes(r,obs))+geom_line(colour=c("#000000"))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.5, colour=c("#e0e0e0"))+xlab("Distance r (m)") +ylab("L(r)")+geom_hline(yintercept=0, linetype = "dashed", colour=c("#ff0000"))+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.title.y = element_text(face = "italic"))
Lplot<-ggplot(Ltemp,aes(r,obs))+geom_line(colour=c("#000000"))+
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.5, colour=c("#e0e0e0"))+
  xlab("Distance r (m)") +ylab("L(r)")+
  geom_hline(yintercept=0, linetype = "dashed", colour=c("#ff0000"))+theme_bw()+
  theme(axis.title.y = element_text(face = "italic"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(color = "black"))+
  geom_rug(data=Ltemp[Ltemp$obs > Ltemp$hi,], sides="b", colour=c("#d73027"))+
  geom_rug(data=Ltemp[Ltemp$obs < Ltemp$lo,], sides="b", colour=c("#ffffbf"))+
  geom_rug(data=Ltemp[Ltemp$obs >= Ltemp$lo & Ltemp$obs <= Ltemp$hi,], sides="b", color=c("#91bfdb"))
Lplot



#單變量單一物種
unionesp<-function(spname,datatype,year,dbhlow,dbhhigh){
  if(datatype=="all"){
    onedt<-dbGetQuery(con,paste("select * from ljc2018 where ba",year,">0",sep = ""))
    onedt$sp<-iconv(onedt$sp,"UTF-8","CP950")
    dbhlow=""
    dbhhigh=""
  }else if(datatype=="dbh"){
    onedt<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=",dbhlow,"and dbh",year,"<",dbhhigh,sep = ""))
    onedt$sp<-iconv(onedt$sp,"UTF-8","CP950")
  }else{
    stop("datatype wrong")
  }
  spdt<-filter(onedt,sp==spname)
  if(nrow(spdt)>=15){
    mypattern<-ppp(spdt$x4,spdt$y4,c(0,300),c(0,100))
    L<-envelope(mypattern, Lest, nsim = 99,correction="best")
    png(filename = paste(datatype,year,"_",spname,"_",dbhlow,"-",dbhhigh,".png",sep = ""),width = 800, height = 645, units = "px")
    par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
    plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = spname,legend=FALSE)
    dev.off()
  }else{
    print(paste(spname,"有",nrow(spdt),"株","can't run"))}
}



#死亡植株
#year1,上次調查 
#year2,當次調查
deadtree<-function(year1,year2,dbh){
  deaddt<-sqldf(paste("select * from db where dbh",year1,">=",dbh," and dbh",year2,"=0",sep = ""))
  if(nrow(deaddt)>=15){
    mypattern<-ppp(deaddt$x4,deaddt$y4,c(0,300),c(0,100))
    L<-envelope(mypattern, Lest, nsim = 999,correction="best")
    year3<-""
    if(year2=="97"){
      year3<-"1997"
    }else if(year2=="05"){
      year3<-"2005"
    }else if(year2=="13"){
      year3<-"2013"
    }
    png(filename = paste("dead",year3,"_",dbh,".png",sep = ""),width = 800, height = 645, units = "px")
    par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
    plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
    dev.off()
  }else{
    print(paste("只有",nrow(deaddt),"株","can't run"))}
}

#新增植株
newtree<-function(year1,year2,species){
  newdt<-sqldf(paste("select * from db where sp='",species,"' and dbh",year1,"=0 and dbh",year2,">0",sep = ""))
  if(nrow(newdt)>=15){
    mypattern<-ppp(newdt$x4,newdt$y4,c(0,300),c(0,100))
    L<-envelope(mypattern, Lest, nsim = 999,correction="best")
    png(filename = paste(species,"new",year2,".png",sep = ""),width = 800, height = 645, units = "px")
    par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
    year3<-""
    if(year2=="97"){
      year3<-"1997"
    }else if(year2=="05"){
      year3<-"2005"
    }else if(year2=="13"){
      year3<-"2013"
    }
    plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
    mtext(paste("n=",nrow(newdt),sep = ""),line = -2.5,cex = 3,adj = 0.02)
    print(paste("正在產生",species,sep=""))
    dev.off()
  }else{
    print(paste(species,"只有",nrow(newdt),"株","can't run"))}
}

newyear1<-matrix(c("05"))
newyear2<-matrix(c("13"))
for(i in 1:length(newyear1)){
  for(j in 1:length(splist)){
      newtree(newyear1[i],newyear2[i],splist[j])
  }
}

#新增所有
newall<-sqldf("select * from db where dbh05=0 and dbh13>0")
if(nrow(newall)>=15){
  mypattern<-ppp(newall$x4,newall$y4,c(0,300),c(0,100))
  L<-envelope(mypattern, Lest, nsim = 999,correction="best")
  png(filename = "newall2013.png",width = 800, height = 645, units = "px")
  par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
  plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
  mtext("(a)",adj = 0,cex = 2,line = 1)
  dev.off()
  }else{}


#存活所有
newall<-sqldf("select * from db where dbh05>0 and dbh13>8")
if(nrow(newall)>=15){
  mypattern<-ppp(newall$x4,newall$y4,c(0,300),c(0,100))
  L<-envelope(mypattern, Lest, nsim = 999,correction="best")
  png(filename = "survive2013.png",width = 800, height = 645, units = "px")
  par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
  plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
  mtext("(a)",adj = 0,cex = 2,line = 1)
  dev.off()
}else{}


#新增與存活
dt1<-sqldf("select * from db where dbh05=0 and dbh13>0")
dt2<-sqldf("select * from db where dbh05>0 and dbh13>=8")
bidt<-rbind(dt1,dt2)
bidt$status05<-ifelse(bidt$dbh05==0,"new","survive")
m<-as.factor(bidt$status05)
mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
L<-envelope(mypattern,Lcross,nsim = 999,i="new",j="survive")
png(filename = "newsurvive2013.png",width = 800,height = 645,units = "px")
par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = "",legend=FALSE)
dev.off()

#新增與死亡
dt1<-sqldf("select * from db where dbh05=0 and dbh13>0")
dt2<-sqldf("select * from db where dbh05>=8 and dbh13=0")
bidt<-rbind(dt1,dt2)
bidt$status13<-ifelse(bidt$dbh13==0,"new","dead")
m<-as.factor(bidt$status13)
mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
L<-envelope(mypattern,Lcross,nsim = 999,i="new",j="dead")
png(filename = "newdead2013.png",width = 800,height = 645,units = "px")
par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = "",legend=FALSE)
dev.off()


#單變量分析迴圈
#datatype,資料種類 *all,所有存活 *dbh,不同徑級
#year,資料年分 *91,97,05,13
#dbhlow,徑級下界
#dbhhigh,徑級上界
uniloop<-function(datatype,year,dbhlow=0,dbhhigh=100,yearplot=""){
  if(datatype=="all"){
    loopdt<-dbGetQuery(con,paste("select * from ljc2018 where ba",year,">0",sep = ""))
    loopdt$sp %<>% iconv(.,"UTF-8","CP950")
  }else if(datatype=="dbh"){
    loopdt<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=",dbhlow," and dbh",year,"<",dbhhigh,sep = ""))
    loopdt$sp %<>% iconv(.,"UTF-8","CP950")
  }else{
    stop("datatype wrong")
  }
  ksp<-matrix()
  k<-1
  for(i in 1:length(splist)){
    spatialdt<-filter(loopdt,sp==splist[i])
    if(nrow(spatialdt)>=15){
      print(paste(splist[i],"分析中",",分析第",k,"個物種"))
      ksp[k]<-splist[i]  
      k=k+1
      if(datatype=="all"){
        png(filename = paste(datatype,year,"_",splist[i],".png",sep = ""),width = 800, height = 645, units = "px")
      }else if(dbhlow==8){
      png(filename = paste(datatype,year,"_",splist[i],"DBH≧",dbhlow,"cm",".png",sep = ""),width = 800, height = 645, units = "px")
      }else{
      png(filename = paste(datatype,year,"_",splist[i],"DBH ",dbhlow,"-",dbhhigh," cm",".png",sep = ""),width = 800, height = 645, units = "px")
      }
      par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
      mypattern<-ppp(spatialdt$x4,spatialdt$y4,c(0,300),c(0,100))
      L<-envelope(mypattern,Lest, nsim = 99,correction="best")
      if(year=="91"){
        yearplot="1991"
      }else if(year=="97"){
        yearplot="1997"
      }else if(year=="05"){
        yearplot="2005"
      }else if(year=="13"){
        yearplot="2013"
      }
      if(datatype=="all" & dbhlow==0){
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,splist[i],sep = "") ,legend=FALSE)
        dev.off()
      }else if(dbhlow==8){
        #輸出原始數據
        #write.table(L, file = paste(yearplot,splist[i],"DBH≧ ",dbhlow,"cm.csv",sep = "") , sep = ",")
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,splist[i],"DBH≧ ",dbhlow,"cm",sep = "") ,legend=FALSE)
        dev.off()
      }else if(datatype=="dbh" & dbhlow>0 & dbhlow!=8){
        #write.table(L, file = paste(yearplot,splist[i],"DBH ",dbhlow,"-",dbhhigh," cm.csv",sep = ""), sep = ",")
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,splist[i],"DBH ",dbhlow,"-",dbhhigh," cm",sep = "") ,legend=FALSE)
        dev.off()
      }else{
        dev.off()
        stop("datatype wrong")
      }
 
    }else{print(paste(splist[i],"有",nrow(spatialdt),"株","can't run"))}
  }
  writeLines(paste(".\n共分析了",k-1,"個物種,包含:"))
  ksp
}



#單變量dbh迴圈2
unidbh1<-matrix(c("1","2","4","8"))
unidbh2<-matrix(c("2","4","8","100"))
yearlist<-matrix(c("91","97","05","13"))
for(i in 1:length(yearlist)){
  for(j in 1:length(unidbh1)){
    uniloop("dbh",yearlist[i],unidbh1[j],unidbh2[j])
  }
}


#test independence
dt1<-dbGetQuery(con,"select * from ljc2018 where dbh91=0 and dbh97=0 and dbh05=0 and dbh13>0")
dt1$sp %<>% iconv(.,"UTF8","CP950")
spdt<-filter(dt1,sp=="長尾栲")
dt2<-dbGetQuery(con,"select * from ljc2018 where dbh05>=8 and dbh13=0")
dt2$sp %<>% iconv(.,"UTF-8","CP950")
bidt<-rbind(spdt,dt2)
bidt$status13<-ifelse(bidt$status13>0,"new","dead")
mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = as.factor(bidt$status13))
pattern<-as.ppp(mypattern)
X<-rshift(pattern,which = "new")
#plot(X)
L<-envelope(X,Lcross,nsim = 99,i="new",j="dead")
L<-envelope(mypattern,Lcross,nsim = 99,i="new",j="dead",envir.simul = rshift(pattern,which = "new"))
plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",legend=FALSE)



#雙變量(類別)
bivar<-function(datatype1,datatype2,year1,year2,spname,dbhbig){
  if(datatype1=="new"){
    if(datatype2=="dead"){
      year3<-""
      dt1<-sqldf(paste("select * from db where dbh",year1,"=0 and dbh",year2,">0",sep = ""))
      spdt<-filter(dt1,sp==spname)
      dt2<-sqldf(paste("select * from db where dbh",year1,">=",dbhbig," and dbh",year2,"=0",sep = ""))
      if(nrow(spdt)>=15 & nrow(dt2)>=15){
        bidt<-rbind(spdt,dt2)
        if(year2=="05"){
          year3<-"2005"
          bidt$status05<-ifelse(bidt$dbh05>0,"new","dead")
          m<-as.factor(bidt$status05)
        }else if(year2=="13"){
          year3<-"2013"
          bidt$status13<-ifelse(bidt$dbh13>0,"new","dead")
          m<-as.factor(bidt$status13)
        }else if(year2=="97"){
          year3<-"1997"
          bidt$status97<-ifelse(bidt$dbh97>0,"new","dead")
          m<-as.factor(bidt$status97)
        }
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 999,i="new",j="dead")
        png(filename = paste(spname,year3,"newdead_",dbhbig,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
        dev.off()
      }else{print(paste(spname,"或死亡大樹資料不足"))}
    }else if(datatype2=="survive"){
      year3<-""
      dt1<-sqldf(paste("select * from db where dbh",year1,"=0 and dbh",year2,">0",sep = ""))
      spdt1<-filter(dt1,sp==spname)
      dt2<-sqldf(paste("select * from db where dbh",year1,">0 and dbh",year2,">=",dbhbig,sep = ""))
      #spdt2<-filter(dt2,sp==spname)
      if(nrow(spdt1)>=15 & nrow(dt2)>=15){
        bidt<-rbind(spdt1,dt2)
        if(year2=="05"){
          year3<-"2005"
          bidt$status05<-ifelse(bidt$dbh97==0,"new","survive")
          m<-as.factor(bidt$status05)
        }else if(year2=="13"){
          year3<-"2013"
          bidt$status13<-ifelse(bidt$dbh05==0,"new","survive")
          m<-as.factor(bidt$status13)
        }else if(year2=="97"){
          year3<-"1997"
          bidt$status97<-ifelse(bidt$dbh91==0,"new","survive")
          m<-as.factor(bidt$status97)
        }
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 999,i="new",j="survive")
        png(filename = paste(spname,year3,"newsurviveall_",dbhbig,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main ="",legend=FALSE)
        dev.off()
      }else{print(paste(spname,"資料不足"))}
    }
  }else{}
}



#雙變量類別迴圈
#datatype *"dead","survive"
newdeadtempsplist<-matrix(splist)
newdeadyearlist1<-matrix(c("05"))
newdeadyearlist2<-matrix(c("13"))
dbhlist<-matrix(c("8"))
bivarloop<-function(datatype){
  for(i in 1:length(newdeadyearlist2)){
    for(j in 1:length(newdeadtempsplist)){
      for(k in 1:length(dbhlist)){
        bivar("new",datatype,newdeadyearlist1[i],newdeadyearlist2[i],newdeadtempsplist[j],dbhlist[k])
      }
    }
  }
}



#雙變量(物種)
bisp<-function(datatype1,datatype2,year,sp1,sp2,dbhlow,dbhhigh,yearplot=""){
  if(datatype1=="all10"){
    if(datatype2=="all10"){
      bidt<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=10",sep = ""))
      bidt$sp<-iconv(bidt$sp,"UTF-8","CP950")
      sp1dt<-filter(bidt,sp==sp1)
      sp2dt<-filter(bidt,sp==sp2)
      if(nrow(sp1dt)>=15 & nrow(sp2dt)>=15){
        bispdt<-rbind(sp1dt,sp2dt)
        m<-as.factor(bispdt$sp)
        mypattern<-ppp(bispdt$x4,bispdt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i=sp1,j=sp2)
        png(filename = paste(year,"all10_",sp1,"-",sp2,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,sp1,"與",sp2,sep=""),legend=FALSE)
        dev.off()
      }else{print("資料不足")}
    }else{}
  }else if(datatype1=="dbh"){
    if(datatype2=="all10"){
      dt1<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=",dbhlow," and dbh",year,"<",dbhhigh,sep = ""))
      dt1$sp<-iconv(dt1$sp,"UTF-8","CP950")
      spdt1<-filter(dt1,sp==sp1)
      dt2<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=10",sep = ""))
      dt2$sp<-iconv(dt2$sp,"UTF-8","CP950")
      spdt2<-filter(dt2,sp==sp2)
      if(sp1!=sp2 & nrow(spdt1)>=15 & nrow(spdt2)>=15){
        bidt<-rbind(spdt1,spdt2)
        m<-as.factor(bidt$sp)
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i=sp1,j=sp2)
        png(filename = paste(year,sp1,dbhlow,"-",dbhhigh,"_",sp2,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(sp1,"DBH ",dbhlow,"-",dbhhigh," cm與",sp2,sep=""),legend=FALSE)
        dev.off()
      }else{print("資料不足")}
    }else if(datatype2=="alive"){
      dt1<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">=",dbhlow," and dbh",year,"<",dbhhigh,sep = ""))
      dt1$sp<-iconv(dt1$sp,"UTF-8","CP950")
      spdt1<-filter(dt1,sp==sp1)
      dt2<-dbGetQuery(con,paste("select * from ljc2018 where dbh",year,">0 except select * from ljc2018 where dbh",year,">=",dbhlow," and dbh",year,"<",dbhhigh,sep = ""))
      dt2$sp<-iconv(dt2$sp,"UTF-8","CP950")
      spdt2<-filter(dt2,sp==sp2)
      if(sp1==sp2 & nrow(spdt1)>=15 & nrow(spdt2)>=15){
        bidt<-rbind(spdt1,spdt2)
        if(year=="91"){
          yearplot<-"1991"
          bidt$stasus91<-ifelse(bidt$dbh91>=as.numeric(dbhlow) & bidt$dbh91<as.numeric(dbhhigh),"sample","alive")
          m<-as.factor(bidt$stasus91)
        }else if(year=="97"){
          yearplot<-"1997"
          bidt$stasus97<-ifelse(bidt$dbh97>=as.numeric(dbhlow) & bidt$dbh97<as.numeric(dbhhigh),"sample","alive")
          m<-as.factor(bidt$stasus97)
        }else if(year=="05"){
          yearplot  <-"2005"
          bidt$stasus05<-ifelse(bidt$dbh05>=as.numeric(dbhlow) & bidt$dbh05<as.numeric(dbhhigh),"sample","alive")
          m<-as.factor(bidt$stasus05)
        }else if(year=="13"){
          yearplot="2013"
          bidt$stasus13<-ifelse(bidt$dbh13>=as.numeric(dbhlow) & bidt$dbh13<as.numeric(dbhhigh),"sample","alive")
          m<-as.factor(bidt$stasus13)
        }
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i="sample",j="alive")
        png(filename = paste(year,sp1,dbhlow,"-",dbhhigh,"_",sp2,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,6,5,3),cex.lab=2,cex.axis=2,cex.main=2)
        if(dbhlow!="25"){
          plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,sp1," ",dbhlow,"-",dbhhigh,"cm",sep=""),legend=FALSE)
          dev.off()
        }else{
          plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(yearplot,sp1," ≧ ",dbhlow,"cm",sep=""),legend=FALSE)
          dev.off()
        }
      }else{print(paste(year,sp1,"資料不足"))}
    }
  }
}



#雙變量大樹間
yearlist<-matrix(c("91","97","05","13"))
tempsplist<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
yearplotlist<-matrix(c("1991","1997","2005","2013"))
for(i in 1:length(yearlist)){
  for(j in 1:length(tempsplist)){
    for(k in j+1:length(tempsplist)){
      if(k>length(tempsplist)){
        
      }else if(tempsplist[j]!=tempsplist[k]){
        bisp("all10","all10",yearlist[i],tempsplist[j],tempsplist[k],yearplot = yearplotlist[i])
      }else{}
    }
  }
}




#雙變量迴圈(物種)(先載入bisp)
bidbhsploop<-function(year){
  tempsplist1<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
  tempsplist2<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
  dbh1<-matrix(c("1","2","4"))
  dbh2<-matrix(c("2","4","6"))
  for(i in 1:length(tempsplist1)){
    for(j in 1:length(tempsplist2)){
      for(k in 1:length(dbh1)){
        bisp("dbh","all10",year,tempsplist1[i],tempsplist2[j],dbh1[k],dbh2[k])
      }
    }
  }
}



#雙變量同物種特定徑級與其他存活
bidbh1<-matrix(c("1","2","4","6","8","10","25"))
bidbh2<-matrix(c("2","4","6","8","10","25","100"))
yearlist<-matrix(c("91","97","05","13"))
tempsplist<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
for(i in 1:length(yearlist)){
  for(j in 1:length(tempsplist)){
    for(k in 1:length(bidbh1)){
      bisp("dbh","alive",yearlist[i],tempsplist[j],tempsplist[j],bidbh1[k],bidbh2[k])
    }
  }
}


#畫點分布圖
env<-read.csv("./欖仁溪環境資料.csv")
newall<-sqldf("select * from db where dbh05=0 and dbh13>0")
mypattern<-ppp(newall$x4,newall$y4,c(0,300),c(0,100))
x<-mypattern[[3]]
y<-mypattern[[4]]
xy<-cbind(x,y)
xy<-as.data.frame(xy)
envf<-filter(env,Y>=29)
envf %<>% mutate(.,X=(X-26)*10,Y=(Y-29)*10)

  #p<-ggplot(xy,aes(x,y))+geom_point()+coord_fixed()+geom_contour(data = envf,aes(x=X,y=Y,z=envf$Eleva))+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
  #p<-ggplot(xy,aes(x,y))+geom_point()+coord_fixed()+geom_contour(data = envf,aes(x=X,y=Y,z=envf$Eleva))+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+scale_x_continuous(sec.axis = sec_axis(~ . /10+26,breaks = seq(26,56,1)),breaks = seq(0,300,50))+scale_y_continuous(sec.axis = sec_axis(~./10+29,breaks = seq(29,39,1)))

p<-ggplot(xy,aes(x,y))+geom_point()+coord_fixed()+
  geom_contour(data = envf,aes(x=X,y=Y,z=envf$Eleva))+geom_rect(data = xy,aes(xmin=0,xmax=300,ymin=0,ymax=100),alpha=0,color="black",size=0.2)+
  scale_x_continuous(sec.axis = sec_axis(~ . /10+26,breaks = seq(26,56,1),name = "Quadrat"),breaks = seq(0,300,50))+
  scale_y_continuous(sec.axis = sec_axis(~./10+29,breaks = seq(29,39,1),name = "Quadrat"))+
  labs(x="Distance(m)",y="Distance(m)")+
  theme_grey()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
png("new2013plot.png",width = 800,height = 400)
p
dev.off()


#05-13物種名錄>=15
spname<-read.csv("./欖仁溪物種名錄.csv")
spname$sp<-as.character(spname$sp)
db$sp<-as.character(db$sp)
sp0513_15<-db %>% filter(.,dbh05==0,dbh13>0) %>% group_by(sp) %>% summarise(num=n()>=15) %>% filter(.,num==TRUE) %>% select(.,sp) %>% left_join(select(spname,sp,species),by="sp") %>% arrange(.,sp)
