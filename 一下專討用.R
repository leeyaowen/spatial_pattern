library(DBI)
library(RPostgreSQL)
library(dplyr) #過濾資料用
library(spatstat)
drv<-dbDriver("PostgreSQL") #宣告資料庫名稱 
con<-dbConnect(drv,user="postgres",password="2717484",dbname="LJC") #建立連線字串
sp15<-dbGetQuery(con,"select sp from ljc3haok group by sp having count(sp)>=15") #產生株數大於15的名錄
sp15$sp<-iconv(sp15$sp,"UTF-8","CP950")
splist<-as.matrix(sp15) #物種清單


#選存活資料
dtall<-dbGetQuery(con,"select * from ljc3haok") #選取所有資料
dtall$sp<-iconv(dtall$sp,"UTF-8","CP950") #sp中文編碼記得轉換

dtall13<-dbGetQuery(con,"select * from ljc3haok where ba13>0") #選2013年存活資料
dtall13$sp<-iconv(dtall13$sp,"UTF-8","CP950") 

dtall05<-dbGetQuery(con,"select * from ljc3haok where ba05>0") #選2005年存活資料
dtall05$sp<-iconv(dtall05$sp,"UTF-8","CP950") 

dtall97<-dbGetQuery(con,"select * from ljc3haok where ba97>0") #選1997年存活資料
dtall97$sp<-iconv(dtall97$sp,"UTF-8","CP950")


#選不同徑級資料
dbh130102<-dbGetQuery(con,"select * from ljc3haok where dbh13>=1 and dbh13<2") #2013年0102
dbh130102$sp<-iconv(dbh130102$sp,"UTF-8","CP950")

dbh130204<-dbGetQuery(con,"select * from ljc3haok where dbh13>=2 and dbh13<4") #2013年0204
dbh130204$sp<-iconv(dbh130204$sp,"UTF-8","CP950")



#單變量單一物種
unionesp<-function(spname,datatype,year,dbhlow,dbhhigh){
  if(datatype=="all"){
    onedt<-dbGetQuery(con,paste("select * from ljc3haok where ba",year,">0",sep = ""))
    onedt$sp<-iconv(onedt$sp,"UTF-8","CP950")
    dbhlow=""
    dbhhigh=""
  }else if(datatype=="dbh"){
    onedt<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year,">=",dbhlow,"and dbh",year,"<",dbhhigh,sep = ""))
    onedt$sp<-iconv(onedt$sp,"UTF-8","CP950")
  }else{
    stop("datatype wrong")
  }
  spdt<-filter(onedt,sp==spname)
  if(nrow(spdt)>=15){
    mypattern<-ppp(spdt$x4,spdt$y4,c(0,300),c(0,100))
    L<-envelope(mypattern, Lest, nsim = 99,correction="best")
    png(filename = paste(datatype,year,"_",spname,"_",dbhlow,"-",dbhhigh,".png",sep = ""),width = 800, height = 645, units = "px")
    par(mar=c(5,5,5,3))
    plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = spname,legend=FALSE)
    dev.off()
  }else{
    print(paste(spname,"有",nrow(spdt),"株","can't run"))}
}



#死亡植株
#year1,上次調查 *97,05(91只有1株死亡大樹)
#year2,當次調查 *05,13
deadtree<-function(year1,year2){
  deaddt<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year1,">=25 and dbh",year2,"=0",sep = ""))
  deaddt$sp<-iconv(deaddt$sp,"UTF-8","CP950")
  if(nrow(deaddt)>=15){
    mypattern<-ppp(deaddt$x4,deaddt$y4,c(0,300),c(0,100))
    L<-envelope(mypattern, Lest, nsim = 99,correction="best")
    png(filename = paste("dead",year2,".png",sep = ""),width = 800, height = 645, units = "px")
    par(mar=c(5,5,5,3))
    year3<-""
    if(year2=="97"){
      year3<-"1997"
    }else if(year2=="05"){
      year3<-"2005"
    }else if(year2=="13"){
      year3<-"2013"
    }
    plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main =paste(year3,"死亡",sep = "") ,legend=FALSE)
    dev.off()
  }else{
    print(paste("只有",nrow(deaddt),"株","can't run"))}
}



#單變量分析迴圈
#datatype,資料種類 *all,所有存活 *dbh,不同徑級
#year,資料年分 *91,97,05,13
#dbhlow,徑級下界
#dbhhigh,徑級上界
uniloop<-function(datatype,year,dbhlow,dbhhigh){
  if(datatype=="all"){
    loopdt<-dbGetQuery(con,paste("select * from ljc3haok where ba",year,">0",sep = ""))
    loopdt$sp<-iconv(loopdt$sp,"UTF-8","CP950")
    dbhlow=""
    dbhhigh=""
  }else if(datatype=="dbh"){
    loopdt<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year,">=",dbhlow,"and dbh",year,"<",dbhhigh,sep = ""))
    loopdt$sp<-iconv(loopdt$sp,"UTF-8","CP950")
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
      png(filename = paste(datatype,year,"_",splist[i],"_",dbhlow,"-",dbhhigh,".png",sep = ""),width = 800, height = 645, units = "px")
      par(mar=c(5,5,5,3))
      mypattern<-ppp(spatialdt$x4,spatialdt$y4,c(0,300),c(0,100))
      L<-envelope(mypattern,Lest, nsim = 99,correction="best")
      plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = splist[i],legend=FALSE)
      dev.off()
    }else{print(paste(splist[i],"有",nrow(spatialdt),"株","can't run"))}
  }
  writeLines(paste(".\n共分析了",k-1,"個物種,包含:"))
  ksp
}



#雙變量(類別)
bivar<-function(datatype1,datatype2,year1,year2,spname){
  if(datatype1=="new"){
    if(datatype2=="dead"){
      year3<-""
      dt1<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year1,"=0 and dbh",year2,">0",sep = ""))
      dt1$sp<-iconv(dt1$sp,"UTF-8","CP950")
      spdt<-filter(dt1,sp==spname)
      dt2<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year1,">=25 and dbh",year2,"=0",sep = ""))
      dt2$sp<-iconv(dt2$sp,"UTF-8","CP950")
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
        }
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i="new",j="dead")
        png(filename = paste(year3,"newdead_",spname,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,5,5,3))
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(spname,"新增與死亡大樹",sep=""),legend=FALSE)
        dev.off()
      }else{print(paste(spname,"或死亡大樹資料不足"))}
    }else if(datatype2=="survive"){
      year3<-""
      dt1<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year1,"=0 and dbh",year2,">0",sep = ""))
      dt1$sp<-iconv(dt1$sp,"UTF-8","CP950")
      spdt1<-filter(dt1,sp==spname)
      dt2<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year1,">0 and dbh",year2,">10",sep = ""))
      dt2$sp<-iconv(dt2$sp,"UTF-8","CP950")
      spdt2<-filter(dt2,sp==spname)
      if(nrow(spdt1)>=15 & nrow(spdt2)>=15){
        bidt<-rbind(spdt1,spdt2)
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
        L<-envelope(mypattern,Lcross,nsim = 99,i="new",j="survive")
        png(filename = paste(year3,"newsurvive_",spname,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,5,5,3))
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(spname,"新增與存活",sep=""),legend=FALSE)
        dev.off()
      }else{print(paste(spname,"資料不足"))}
    }
  }else{}
}



#雙變量(物種)
bisp<-function(datatype1,datatype2,year,sp1,sp2,dbhlow,dbhhigh){
  if(datatype1=="all10"){
    if(datatype2=="all10"){
      bidt<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year,">=10",sep = ""))
      bidt$sp<-iconv(bidt$sp,"UTF-8","CP950")
      sp1dt<-filter(bidt,sp==sp1)
      sp2dt<-filter(bidt,sp==sp2)
      if(nrow(sp1dt)>=15 & nrow(sp2dt)>=15){
        bispdt<-rbind(sp1dt,sp2dt)
        m<-as.factor(bispdt$sp)
        mypattern<-ppp(bispdt$x4,bispdt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i=sp1,j=sp2)
        png(filename = paste(year,"all10_",sp1,"-",sp2,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,5,5,3))
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(sp1,"與",sp2,sep=""),legend=FALSE)
        dev.off()
      }else{print("資料不足")}
    }else{}
  }else if(datatype1=="dbh"){
    if(datatype2=="all10"){
      dt1<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year,">=",dbhlow," and dbh",year,"<",dbhhigh,sep = ""))
      dt1$sp<-iconv(dt1$sp,"UTF-8","CP950")
      spdt1<-filter(dt1,sp==sp1)
      dt2<-dbGetQuery(con,paste("select * from ljc3haok where dbh",year,">=10",sep = ""))
      dt2$sp<-iconv(dt2$sp,"UTF-8","CP950")
      spdt2<-filter(dt2,sp==sp2)
      if(sp1!=sp2 & nrow(spdt1)>=15 & nrow(spdt2)>=15){
        bidt<-rbind(spdt1,spdt2)
        m<-as.factor(bidt$sp)
        mypattern<-ppp(bidt$x4,bidt$y4,c(0,300),c(0,100),marks = m)
        L<-envelope(mypattern,Lcross,nsim = 99,i=sp1,j=sp2)
        png(filename = paste(year,sp1,dbhlow,"-",dbhhigh,"_",sp2,".png",sep = ""),width = 800, height = 645, units = "px")
        par(mar=c(5,5,5,3))
        plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = paste(sp1,"0",dbhlow,"0",dbhhigh,"與",sp2,sep=""),legend=FALSE)
        dev.off()
      }else{print("資料不足")}
    }else{}
  }
}



#雙變量迴圈(物種)
biloop<-function(year){
  tempsplist1<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
  tempsplist2<-matrix(c("紅花八角","長尾栲","奧氏虎皮楠","南仁灰木"))
  dbh1<-matrix(c("1","2","4"))
  dbh2<-matrix(c("2","4","6"))
  for(i in 1:length(tempsplist1)){
    for(j in 1:length(tempsplist2)){
      for(k in 1:length(dbh1)){
        bisp("dbh","all10",year,tempsplist[i],tempsplist[j],dbh1[k],dbh2[k])
      }
    }
  }
}