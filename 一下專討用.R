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
plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = spname,legend=FALSE)
}else{
  print(paste(spname,"有",nrow(spdt),"株","can't run"))}
}



#分析迴圈
#datatype,資料種類 *all,所有存活 *dbh,不同徑級
#year,資料年分 *91,97,05,13
#dbhlow,徑級下界
#dbhhigh,徑級上界
uniloop<-function(datatype,year,dbhlow,dbhhigh){
  if(datatype=="all"){
    loopdt<-dbGetQuery(con,paste("select * from ljc3haok where ba",year,">0",sep = ""))
    loopdt$sp<-iconv(loopdt$sp,"UTF-8","CP950")
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
  mypattern<-ppp(spatialdt$x4,spatialdt$y4,c(0,300),c(0,100))
  L<-envelope(mypattern,Lest, nsim = 99,correction="best")
  plot(L,.-r~r,ylab = expression(L(r)),xlab = "d(m)",main = splist[i],legend=FALSE)
}else{print(paste(splist[i],"有",nrow(spatialdt),"株","can't run"))}
}
writeLines(paste(".\n共分析了",k,"個物種,包含:"))
ksp
}  