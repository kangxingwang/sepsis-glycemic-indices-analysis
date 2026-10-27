fit_cox<-function(d,rhs,response="Surv(Day28_LOS,Day28_Outcome)"){

  f<-as.formula(paste(response,"~",paste(rhs,collapse="+"))); z<-withCallingHandlers(survival::coxph(f,data=d,ties="efron",cluster=stay_id,robust=TRUE,x=TRUE,y=TRUE,model=TRUE,na.action=na.fail,singular.ok=FALSE,control=survival::coxph.control(iter.max=50)),warning=function(w) stop(conditionMessage(w),call.=FALSE))

  if(any(!is.finite(coef(z)))||any(!is.finite(vcov(z)))) stop("Cox估计无效。"); z

}

pool_mi<-function(fits){

  nm<-names(coef(fits[[1]])); m<-length(fits); k<-length(nm); stopifnot(m>=2,all(vapply(fits,function(f) identical(names(coef(f)),nm),logical(1))))

  Q<-matrix(unlist(lapply(fits,coef)),k,dimnames=list(nm,NULL)); U<-Reduce("+",lapply(fits,vcov))/m; B<-matrix(cov(t(Q)),k,k,dimnames=list(nm,nm)); df<-min(vapply(fits,function(f) f$nevent-length(coef(f)),numeric(1)))

  if(df<=1) stop("模型事件数不足。"); list(q=rowMeans(Q),U=U,B=B,m=m,df=df,names=nm)

}

pool_contrast<-function(P,C=diag(length(P$q))){

  C<-as.matrix(C); b<-drop(C%*%P$q); u<-rowSums((C%*%P$U)*C); B<-pmax(rowSums((C%*%P$B)*C),0); r<-ifelse(u>0,(1+1/P$m)*B/u,0)

  df0<-ifelse(r>0,(P$m-1)*(1+1/r)^2,Inf); df1<-(P$df+1)/(P$df+3)*P$df/(1+r); df<-1/(1/df0+1/df1); se<-sqrt(pmax(u+(1+1/P$m)*B,0))

  data.frame(HR=exp(b),Lower=exp(b-qt(.975,df)*se),Upper=exp(b+qt(.975,df)*se),P_value=ifelse(se>0,2*pt(-abs(b/se),df),NA_real_))

}

joint_p<-function(P,terms){

  j<-match(terms,P$names); k<-length(j); if(!k||anyNA(j)) stop("检验项匹配失败。"); if(k==1) return(pool_contrast(P,diag(length(P$q))[j,,drop=FALSE])$P_value)

  # 多参数总体P：D1合并Wald检验，大样本近似。

  u<-P$U[j,j,drop=FALSE]; b<-P$B[j,j,drop=FALSE]; q<-P$q[j]; r<-max(0,(1+1/P$m)*sum(diag(solve(u,b)))/k); F<-as.numeric(crossprod(q,solve(u,q)))/(k*(1+r)); t<-k*(P$m-1)

  df2<-if(r<1e-12) Inf else if(t>4) 4+(t-4)*(1+(1-2/t)/r)^2 else t*(1+1/k)*(1+1/r)^2/2; pf(F,k,df2,lower.tail=FALSE)

}

cox_table<-function(fits,e,form,model,interval){

  P<-pool_mi(fits); term<-if(form=="continuous") "Exposure" else "Exposure_Q"; ix<-grep(paste0("^",term),P$names); z<-pool_contrast(P)[ix,,drop=FALSE]

  z$Comparison<-if(form=="continuous") unit_labels[e] else paste(sub("Exposure_Q","",P$names[ix],fixed=TRUE),"vs Q1"); if(form=="quartile") z<-rbind(data.frame(HR=1,Lower=1,Upper=1,P_value=NA_real_,Comparison="Q1 (reference)"),z)

  f<-fits[[1]]; z$N<-f$n; z$Events<-f$nevent; if(form=="quartile"){g<-f$model$Exposure_Q; qs<-substr(z$Comparison,1,2); z$N<-as.integer(table(g)[qs]); z$Events<-as.integer(tapply(f$y[,ncol(f$y)],g,sum)[qs])}

  z$P_overall<-joint_p(P,P$names[ix]); cbind(Exposure=e,Form=form,Model=model,Interval=interval,N_at_risk=f$n,Events_interval=f$nevent,z,row.names=NULL)

}

add_holm<-function(z,by){

  keys<-c("Exposure","Form","Model","Interval"); a<-unique(z[c(keys,"P_overall")]); a$P_overall_Holm<-ave(a$P_overall,do.call(interaction,a[by]),FUN=function(p) p.adjust(p,"holm"))

  key<-function(d) do.call(paste,c(d[keys],sep="__")); z$P_overall_Holm<-a$P_overall_Holm[match(key(z),key(a))]; z

}



