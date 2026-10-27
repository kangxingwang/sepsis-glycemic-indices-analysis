if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent=TRUE)
Sys.setenv(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
# 接续01特征工程：MIMIC-IV训练 → 内部测试 / NWICU外部验证 → XGBoost SHAP。
# 从第1节开始，每次选中完整一节运行。正式运行较慢；参数、路径和图片尺寸统一在第1节修改。
# 导出CSV和RDS；图形由后续R-only步骤统一生成。RDS及逐折文件提供正式断点续跑。

# ==================== 1. 导入包、路径和参数 ====================
pkgs<-c("data.table","mlr3","mlr3learners","glmnet","Boruta","ranger","xgboost","e1071","kknn","pROC","ggplot2","ggbeeswarm","viridisLite","ragg","scales","openxlsx")
miss<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]; if(length(miss)) stop("Install the required R packages before running: ",paste(miss,collapse=", "))
suppressPackageStartupMessages({library(data.table);library(mlr3);library(mlr3learners);library(ggplot2)})
fe_file<-file.path(Sys.getenv("SEPSIS_ML_OUTPUT"), "private_checkpoints/特征工程结果.rds"); out_dir<-Sys.getenv("SEPSIS_ML_OUTPUT"); png_dir<-file.path(out_dir,"figures")
dir.create(png_dir,recursive=TRUE,showWarnings=FALSE); model_file<-file.path(out_dir,"private_checkpoints/机器学习模型.rds"); excel_file<-file.path(out_dir,"机器学习结果.xlsx")
quick_test<-FALSE; resume<-TRUE; seed<-1117L; threads<-4L; set.seed(seed)   # quick_test=TRUE仅用于检查代码，不能用于论文结果
parallel_tuning<-TRUE; max_workers<-3L                     # 并行执行不同CV折；树模型每个进程仍使用threads线程
cv_folds<-if(quick_test) 3L else 10L; tune_n<-if(quick_test) c(xgboost=2,ranger=2,svm=2,kknn=2) else c(xgboost=30,ranger=20,svm=20,kknn=20)
rf_trees<-if(quick_test) 80L else 1500L; bootstrap_n<-if(quick_test) 30L else 1000L; use_platt<-TRUE
fig_width<-12; fig_height<-10; fig_dpi<-600; font_family<-"Arial"; shap_top_n<-30L # 全部单张PNG尺寸（英寸）
roc_width<-fig_width; roc_height<-fig_height; shap_bar_width<-fig_width; shap_width<-fig_width; shap_height<-fig_height
plot_font<-font_family; roc_font_size<-14; shap_font_size<-roc_font_size; shap_point_size<-1.1 # 统一白底、浅灰网格、Arial
labels<-c(log_reg="Logistic regression",xgboost="XGBoost",ranger="Random forest",svm="SVM-RBF",kknn="k-NN")
cols<-setNames(viridisLite::plasma(length(labels)),unname(labels)); roc_cols<-cols; shap_bar_col<-unname(cols["XGBoost"])
logmsg<-function(...){cat(format(Sys.time(),"%H:%M:%S"),...,"\n");flush.console()}; clip<-function(p) pmin(pmax(as.numeric(p),1e-6),1-1e-6)
save_png<-function(p,name,w=fig_width,h=fig_height){ggsave(file.path(png_dir,name),p,device=ragg::agg_png,width=w,height=h,units="in",dpi=fig_dpi,bg="white");if(interactive()) print(p)}

# ==================== 2. 读取上一阶段数据及全训练集最终交集变量 ====================
if(!file.exists(fe_file)) stop("请先完成01特征工程：",fe_file)
z<-readRDS(fe_file); selected<-z$selected; if(!length(selected)) stop("LASSO与Boruta交集为空，不能建模")
stopifnot(!any(z$raw$train$stay_id%in%z$raw$internal$stay_id),all(c("train","internal","nwicu")%in%names(z$data)))
ids<-lapply(z$raw,`[[`,"stay_id"); truth<-lapply(z$raw,function(d) as.integer(as.character(d$Outcome)))
stopifnot(all(vapply(truth,function(y) !anyNA(y)&&all(y%in%0:1)&&length(unique(y))==2,logical(1))))
candidates<-z$candidates; factor_vars<-names(z$prep$levels); fs<-z$settings
if(quick_test){fs$boruta_runs<-25L;fs$boruta_trees<-60L;fs$nfolds<-3L}
cfg<-list(version="ML_SHAP_first_patient_v2",quick_test=quick_test,seed=seed,cv_folds=cv_folds,tune_n=tune_n,rf_trees=rf_trees,threads=threads,bootstrap_n=bootstrap_n,use_platt=use_platt,feature_settings=fs)
cat("最终建模变量（",length(selected),"项）：",paste(selected,collapse=", "),"\n");print(z$cohort_table,row.names=FALSE)
if(!"HGI"%in%selected) cat("HGI未通过上一阶段交集筛选，本版模型不包含HGI。\n")

# ==================== 3. 每个训练折独立预处理和特征筛选 ====================
# 保留上一阶段规则。此处每折从43个候选中重新筛选；全训练集交集变量仅用于最终模型。
# 分折后重新估计HGI公式、插补、缩尾和编码；内部测试/NWICU不参与这些估计。
make_folds<-function(y,k,s){if(min(table(y))<k) stop("较少结局类别不足以分折");set.seed(s);f<-integer(length(y));for(g in unique(y)){ii<-sample(which(y==g));f[ii]<-rep(seq_len(k),length.out=length(ii))};f}
mode_value<-function(x){tt<-table(x);if(!length(tt)) stop("训练列全部缺失");names(tt)[which.max(tt)]}
fit_prep<-function(d){
  h<-coef(lm(HbA1c~Glucose,data=d));if(any(!is.finite(h))) stop("HGI公式无法估计");d$HGI<-d$HbA1c-(h[1]+h[2]*d$Glucose)
  nums<-setdiff(candidates,factor_vars);med<-vapply(d[nums],median,numeric(1),na.rm=TRUE);if(any(!is.finite(med))) stop("某训练折存在全缺失列")
  modes<-vapply(d[factor_vars],mode_value,character(1));lev<-lapply(d[factor_vars],function(x) sort(unique(na.omit(as.character(x)))))
  for(v in names(lev)){ref<-switch(v,Gender="F",Race="White",ICU_Type="Medical ICU","No");if(ref%in%lev[[v]]) lev[[v]]<-c(ref,setdiff(lev[[v]],ref))}
  list(hgi=h,median=med,mode=modes,levels=lev,bounds=lapply(d[fs$winsor_vars],quantile,probs=fs$winsor_probs,na.rm=TRUE,names=FALSE),vars=candidates)
}
apply_prep<-function(d,p){
  d$HGI<-d$HbA1c-(p$hgi[1]+p$hgi[2]*d$Glucose)
  for(v in names(p$median)){d[[v]][is.na(d[[v]])]<-p$median[[v]];if(v%in%names(p$bounds)) d[[v]]<-pmin(pmax(d[[v]],p$bounds[[v]][1]),p$bounds[[v]][2])}
  for(v in names(p$levels)){a<-as.character(d[[v]]);a[is.na(a)|!a%in%p$levels[[v]]]<-p$mode[[v]];d[[v]]<-factor(a,levels=p$levels[[v]])}
  d<-d[p$vars];if(anyNA(d)) stop("预处理后仍有缺失");d
}
encode<-function(d){tt<-terms(~.,data=d);mm<-model.matrix(tt,d);a<-attr(mm,"assign");keep<-a!=0;parent<-attr(tt,"term.labels")[a[keep]];mm<-mm[,keep,drop=FALSE];storage.mode(mm)<-"double";attr(mm,"parent")<-setNames(parent,colnames(mm));mm}
fit_scale<-function(x){ss<-apply(x,2,sd);keep<-names(ss)[is.finite(ss)&ss>0];if(!length(keep)) stop("没有可用设计列");list(columns=keep,center=colMeans(x[,keep,drop=FALSE]),scale=ss[keep])}
apply_scale<-function(x,p){stopifnot(all(p$columns%in%colnames(x)));a<-sweep(sweep(x[,p$columns,drop=FALSE],2,p$center,"-"),2,p$scale,"/");stopifnot(all(is.finite(a)));a}
select_fold<-function(it,iv,k){
  d<-z$raw$train[it,,drop=FALSE];p<-fit_prep(d);a<-apply_prep(d,p);p$vars<-names(a)[vapply(a,function(v) length(unique(v))>1,logical(1))]
  a<-apply_prep(d,p);b<-apply_prep(z$raw$train[iv,,drop=FALSE],p);if(ncol(a)<2) stop("第",k,"折可用变量不足")
  x<-encode(a);pm<-attr(x,"parent");f<-make_folds(d$Outcome,fs$nfolds,seed+100*k);set.seed(seed+100*k)
  cv<-glmnet::cv.glmnet(x,d$Outcome,family="binomial",alpha=1,standardize=TRUE,foldid=f,type.measure="deviance")
  co<-as.matrix(coef(cv,s=fs$lasso_rule))[,1];ls<-unique(unname(pm[setdiff(names(co)[co!=0],"(Intercept)")]))
  set.seed(seed+100*k);bo<-Boruta::Boruta(x=a,y=factor(d$Outcome,levels=c(0,1)),pValue=fs$boruta_p,mcAdj=fs$boruta_mcAdj,maxRuns=fs$boruta_runs,
    doTrace=0,holdHistory=TRUE,getImp=Boruta::getImpRfZ,num.trees=fs$boruta_trees,num.threads=fs$boruta_threads,respect.unordered.factors=fs$boruta_factors)
  yes<-names(bo$finalDecision)[bo$finalDecision=="Confirmed"];sel<-p$vars[p$vars%in%intersect(ls,yes)]
  if(!length(sel)) stop("第",k,"折交集为空；停止本次CV，请审核，不会自动放宽筛选规则")
  xt<-encode(a[sel]);xv<-encode(b[sel]);sc<-fit_scale(xt);stopifnot(identical(colnames(xt),colnames(xv)))
  list(it=it,iv=iv,x_train=apply_scale(xt,sc),x_valid=apply_scale(xv,sc),selected=sel,prep=p,scale=sc,
    selection=data.frame(Fold=k,Variable=p$vars,LASSO=p$vars%in%ls,Boruta=as.character(bo$finalDecision[p$vars]),Intersection=p$vars%in%sel))
}

# ==================== 4. 进度保存和训练集分折 ====================
# 只有此RDS用于续跑；修改流程或想完全重跑时，把第1节resume改为FALSE。
key<-list(input_md5=unname(tools::md5sum(fe_file)),cfg=cfg,prep=deparse(fit_prep),apply=deparse(apply_prep),selection=deparse(select_fold),encode=deparse(encode))
state<-list(key=key,folds=list(),tuning=list(),models=list(),calibrators=list(),thresholds=list(),oof=list(),predictions=list(),best_params=list())
if(resume&&file.exists(model_file)){old<-readRDS(model_file);if(identical(old$key,key)){state<-old;logmsg("读取已完成进度")} else logmsg("输入/参数已变化，将重新计算")}
save_state<-function(){tmp<-paste0(model_file,".tmp");saveRDS(state,tmp);if(file.exists(model_file))unlink(model_file);stopifnot(file.rename(tmp,model_file))}
foldid<-make_folds(truth$train,cv_folds,seed+1L)
fold_cache_dir<-file.path(out_dir,"private_checkpoints/折内筛选断点");dir.create(fold_cache_dir,showWarnings=FALSE)
for(k in seq_len(cv_folds)){
  cache<-file.path(fold_cache_dir,paste0("fold",k,".rds"))
  if(is.null(state$folds[[paste0("fold",k)]])&&file.exists(cache)){
    cached<-readRDS(cache);if(identical(cached$key,key))state$folds[[paste0("fold",k)]]<-cached$value
  }
}
pending_folds<-which(vapply(seq_len(cv_folds),function(k)is.null(state$folds[[paste0("fold",k)]]),logical(1)))
if(length(pending_folds)){
  logmsg("正式折内筛选：",length(pending_folds),"折；最多3进程，每进程4线程；逐折缓存")
  fold_cl<-parallel::makePSOCKcluster(min(3L,length(pending_folds)),outfile="")
  tryCatch({
    parallel::clusterCall(fold_cl,function(lib){if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent=TRUE);.libPaths(lib);Sys.setenv(OMP_NUM_THREADS="4",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1");suppressPackageStartupMessages(library(glmnet));suppressPackageStartupMessages(library(Boruta));NULL},.libPaths())
    parallel::clusterExport(fold_cl,c("z","candidates","factor_vars","fs","seed","fit_prep","apply_prep","mode_value","encode","fit_scale","apply_scale","make_folds","select_fold","foldid","fold_cache_dir","key"),envir=environment())
    rr<-parallel::parLapplyLB(fold_cl,pending_folds,function(k){
      cat(format(Sys.time(),"%H:%M:%S"),"start fold",k,"\n");flush.console()
      value<-select_fold(which(foldid!=k),which(foldid==k),k)
      saveRDS(list(key=key,value=value),file.path(fold_cache_dir,paste0("fold",k,".rds")))
      cat(format(Sys.time(),"%H:%M:%S"),"finished fold",k,"selected",length(value$selected),"\n");flush.console()
      list(k=k,value=value)
    })
    for(a in rr){state$folds[[paste0("fold",a$k)]]<-a$value;save_state()}
  },finally={parallel::stopCluster(fold_cl)})
}
state$folds<-state$folds[paste0("fold",seq_len(cv_folds))]
stopifnot(length(state$folds)==cv_folds,all(vapply(state$folds,function(f)!is.null(f),logical(1))))
fold_selection<-do.call(rbind,lapply(state$folds,`[[`,"selection"));stability<-data.frame(Variable=candidates,Selected_folds=vapply(candidates,function(v) sum(fold_selection$Variable==v&fold_selection$Intersection),integer(1)),Total_folds=cv_folds)
print(stability,row.names=FALSE)

# ==================== 5. mlr3学习器与调参范围 ====================
# 显式循环CV，让每一折使用自己的预处理和变量清单；学习器均使用mlr3。
# 延续模板的PR-AUC调参方向，具体用Average Precision（AP）；CV值仅用于训练内选择。
ap_score<-function(y,p){o<-order(p,decreasing=TRUE);yy<-y[o];pp<-p[o];end<-which(!duplicated(pp,fromLast=TRUE));tp<-cumsum(yy)[end];sum(diff(c(0,tp/sum(y)))*tp/end)}
auc_score<-function(y,p) as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))
make_task<-function(x,y,id){d<-as.data.frame(x,check.names=FALSE);d$Outcome<-factor(y,levels=c(0,1));mlr3::TaskClassif$new(id,d,target="Outcome",positive="1")}
new_learner<-function(id,p){base<-switch(id,
  log_reg=list(maxit=200),xgboost=list(nthread=threads,objective="binary:logistic",eval_metric="aucpr",verbose=0),
  ranger=list(num.trees=rf_trees,num.threads=threads,importance="permutation"),
  svm=list(kernel="radial",type="C-classification",scale=FALSE,cachesize=1000),kknn=list(scale=FALSE,store_model=TRUE))
  do.call(mlr3::lrn,c(list(paste0("classif.",id),predict_type="prob"),modifyList(base,p)))
}
fit_model<-function(id,p,x,y,s){set.seed(s);m<-new_learner(id,p);m$train(make_task(x,y,paste0("fit_",id)));m}
predict_prob<-function(m,x){p<-as.numeric(m$predict_newdata(as.data.frame(x,check.names=FALSE))$prob[,"1"]);if(any(!is.finite(p))||any(p<0|p>1)) stop("预测概率无效");p}
draw_params<-function(id){
  if(id=="log_reg") return(list(list()));set.seed(seed+match(id,names(labels)));n<-as.integer(tune_n[id]);rlog<-function(a,b) exp(runif(1,log(a),log(b)))
  lapply(seq_len(n),function(i) switch(id,
    xgboost=list(nrounds=if(quick_test) 50L else sample(100:1200,1),eta=rlog(.005,.30),max_depth=sample(2:8,1),min_child_weight=rlog(1,20),subsample=runif(1,.5,1),colsample_bytree=runif(1,.5,1),gamma=rlog(1e-5,5),lambda=rlog(1e-4,10),alpha=rlog(1e-4,10)),
    ranger=list(mtry.ratio=runif(1,.1,1),min.node.size=sample(1:40,1),sample.fraction=runif(1,.5,1),splitrule=sample(c("gini","extratrees"),1)),
    svm=list(cost=rlog(.001,1000),gamma=rlog(1e-5,10)),kknn=list(k=sample(seq(3,101,2),1),distance=runif(1,1,2),kernel=sample(c("optimal","triangular","gaussian"),1))))
}
eval_cv_job<-function(job){f<-job$fold;p<-job$params;if(job$id=="ranger") p$importance<-"none" # CV不计算未使用的排列重要性；最终模型仍保留
  m<-fit_model(job$id,p,f$x_train,job$y_train,job$seed);pp<-predict_prob(m,f$x_valid)
  list(iv=f$iv,p=pp,ap=ap_score(job$y_valid,pp),auc=auc_score(job$y_valid,pp))}

# ==================== 6. 五个模型在训练集内调参 ====================
# 不按内部测试集或NWICU AUC选择参数/算法；这里使用各折AP的平均值。
for(id in names(labels)){
  pars<-draw_params(id);if(is.null(state$tuning[[id]])) state$tuning[[id]]<-vector("list",length(pars))
  cores<-parallel::detectCores();if(is.na(cores)) cores<-1L
  workers<-if(!parallel_tuning||id=="log_reg") 1L else min(max_workers,cv_folds,max(1L,floor(cores/if(id%in%c("xgboost","ranger")) threads else 1L)))
  cl<-NULL;tryCatch({
    if(workers>1L){cl<-parallel::makePSOCKcluster(workers);parallel::clusterCall(cl,function(lib){.libPaths(lib);suppressPackageStartupMessages(library(mlr3learners));NULL},.libPaths())
      parallel::clusterExport(cl,c("threads","rf_trees","new_learner","make_task","fit_model","predict_prob","ap_score","auc_score","eval_cv_job"),envir=environment())}
    for(j in seq_along(pars)){
      if(!is.null(state$tuning[[id]][[j]])) next
      jobs<-lapply(seq_len(cv_folds),function(k){f<-state$folds[[k]];list(id=id,params=pars[[j]],fold=f,y_train=truth$train[f$it],y_valid=truth$train[f$iv],seed=seed+10000*match(id,names(labels))+100*j+k)})
      rr<-if(is.null(cl)) lapply(jobs,eval_cv_job) else parallel::parLapply(cl,jobs,eval_cv_job)
      pred<-rep(NA_real_,nrow(z$raw$train));for(a in rr) pred[a$iv]<-a$p;ap<-vapply(rr,`[[`,numeric(1),"ap");au<-vapply(rr,`[[`,numeric(1),"auc")
      state$tuning[[id]][[j]]<-list(params=pars[[j]],oof=pred,fold_ap=ap,fold_auc=au,mean_ap=mean(ap),mean_auc=mean(au));save_state()
      logmsg(labels[id],": 参数",j,"/",length(pars),"；CV AP=",round(mean(ap),4),"；CV AUC=",round(mean(au),4))
    }
  },finally={if(!is.null(cl)) parallel::stopCluster(cl)})
  best<-which.max(vapply(state$tuning[[id]],`[[`,numeric(1),"mean_ap"));state$best_params[[id]]<-state$tuning[[id]][[best]]$params
  state$oof[[id]]<-state$tuning[[id]][[best]]$oof;save_state()
}
cv_summary<-do.call(rbind,lapply(names(labels),function(id){a<-state$tuning[[id]];b<-which.max(vapply(a,`[[`,numeric(1),"mean_ap"));data.frame(Model=unname(labels[id]),Best_trial=b,CV_AP=a[[b]]$mean_ap,CV_AUC=a[[b]]$mean_auc)}))
best_model<-names(labels)[which.max(cv_summary$CV_AP)];cat("按训练内CV AP选择的模型：",labels[best_model],"\n");print(cv_summary,row.names=FALSE)

# ==================== 7. 最终训练、OOF校准、内外部预测 ====================
# 最终模型仅拟合MIMIC-IV 70%训练集；同一个固定模型应用到30%内部测试和全部NWICU。
# 沿用模板的训练OOF Platt校准；阈值也从训练OOF确定；测试数据只用于评价。
scale_final<-fit_scale(z$x_train);xx<-list(train=apply_scale(z$x_train,scale_final),internal=apply_scale(z$x_internal,scale_final),nwicu=apply_scale(z$x_nwicu,scale_final))
fit_cal<-function(y,p){if(!use_platt) return(c(intercept=0,slope=1));lp<-qlogis(clip(p));if(sd(lp)<1e-12){warning("OOF预测为常量，校准结果为训练集事件率");return(c(intercept=qlogis(mean(y)),slope=0))};a<-coef(glm(y~lp,family=binomial()));if(any(!is.finite(a))) stop("训练OOF校准失败");setNames(unname(a),c("intercept","slope"))}
cal_prob<-function(a,p) if(use_platt) plogis(a[1]+a[2]*qlogis(clip(p))) else as.numeric(p)
best_cut<-function(y,p){r<-pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE);v<-pROC::coords(r,"best",best.method="youden",ret="threshold",transpose=FALSE)$threshold;v<-v[is.finite(v)];if(!length(v)) stop("训练OOF阈值不可用");sort(v)[ceiling(length(v)/2)]}
for(id in names(labels)){
  if(!is.null(state$models[[id]])&&!is.null(state$predictions[[id]])) next
  logmsg("最终训练：",labels[id]);m<-fit_model(id,state$best_params[[id]],xx$train,truth$train,seed+30000+match(id,names(labels)))
  ca<-fit_cal(truth$train,state$oof[[id]]);cut<-best_cut(truth$train,cal_prob(ca,state$oof[[id]]))
  state$models[[id]]<-m;state$calibrators[[id]]<-ca;state$thresholds[[id]]<-cut
  state$predictions[[id]]<-do.call(rbind,lapply(c("internal","nwicu"),function(co){pr<-predict_prob(m,xx[[co]])
    data.frame(stay_id=ids[[co]],Cohort=co,Model=id,Outcome=truth[[co]],Raw_probability=pr,Probability=cal_prob(ca,pr),Threshold=cut)}));save_state()
}
predictions<-do.call(rbind,state$predictions);rownames(predictions)<-NULL
state$feature_data<-z;state$scale_final<-scale_final;state$x<-xx;state$best_model<-best_model;state$cv_summary<-cv_summary;state$fold_selection<-fold_selection;save_state()

# ==================== 8. 内外部性能、95%CI及校准统计 ====================
# AUC采用DeLong 95%CI；AP和Brier采用结局分层Bootstrap。仅重采样预测，不重训练模型。
safe_div<-function(a,b) if(b>0) a/b else NA_real_
metric_one<-function(d,s){
  y<-d$Outcome;p<-d$Probability;cl<-p>=d$Threshold[1];tp<-sum(cl&y==1);tn<-sum(!cl&y==0);fp<-sum(cl&y==0);fn<-sum(!cl&y==1)
  r<-pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE);ac<-as.numeric(pROC::ci.auc(r,method="delong"));lp<-qlogis(clip(p))
  ci_fit<-glm(y~offset(lp),family=binomial());cs_fit<-glm(y~lp,family=binomial())
  set.seed(s);pos<-which(y==1);neg<-which(y==0);bt<-replicate(bootstrap_n,{ii<-c(sample(pos,length(pos),TRUE),sample(neg,length(neg),TRUE));c(AP=ap_score(y[ii],p[ii]),Brier=mean((p[ii]-y[ii])^2))})
  qc<-apply(bt,1,quantile,probs=c(.025,.975),na.rm=TRUE)
  data.frame(Cohort=d$Cohort[1],Model=d$Model[1],N=length(y),Deaths=sum(y),AUC=as.numeric(pROC::auc(r)),AUC_lower=ac[1],AUC_upper=ac[3],
    AP=ap_score(y,p),AP_lower=qc[1,"AP"],AP_upper=qc[2,"AP"],Brier=mean((p-y)^2),Brier_lower=qc[1,"Brier"],Brier_upper=qc[2,"Brier"],
    Calibration_intercept=unname(coef(ci_fit)[1]),Calibration_slope=unname(coef(cs_fit)[2]),Threshold=d$Threshold[1],
    Sensitivity=safe_div(tp,tp+fn),Specificity=safe_div(tn,tn+fp),PPV=safe_div(tp,tp+fp),NPV=safe_div(tn,tn+fn),F1=safe_div(2*tp,2*tp+fp+fn),Accuracy=mean(cl==y))
}
groups<-split(predictions,interaction(predictions$Cohort,predictions$Model,drop=TRUE));performance<-do.call(rbind,lapply(seq_along(groups),function(i) metric_one(groups[[i]],seed+40000+i)))
rownames(performance)<-NULL;print(performance[,c("Cohort","Model","AUC","AUC_lower","AUC_upper","AP","Brier")],row.names=FALSE)

cohort_labels<-c(internal="MIMIC-IV internal validation",nwicu="NWICU external validation")
plot_sources<-list()
# ==================== 10. XGBoost精确SHAP：按原始临床变量汇总 ====================
# 按模板预先指定XGBoost解释；最优算法另由训练CV决定，不根据外部AUC改选解释模型。
# 多水平分类变量先逐患者加总虚拟列SHAP，再计算mean(abs(SHAP))，避免按虚拟列重复排名。
# 如启用Platt：校准log-odds=a+b×原始margin，故SHAP乘b、基线变为a+b×BIAS。
extract_booster<-function(m){if(inherits(m,"xgb.Booster")) return(m);if(inherits(m,"Learner")) return(extract_booster(m$model));if(is.list(m)) for(a in m){ans<-try(extract_booster(a),silent=TRUE);if(!inherits(ans,"try-error")) return(ans)};stop("未找到XGBoost模型")}
booster<-extract_booster(state$models$xgboost);shap_results<-list();shap_checks<-list()
for(co in names(cohort_labels)){
  x<-xx[[co]];if(!setequal(booster$feature_names,colnames(x))) stop("XGBoost和输入设计列不一致");x<-x[,booster$feature_names,drop=FALSE];dm<-xgboost::xgb.DMatrix(x)
  sh0<-as.matrix(predict(booster,dm,predcontrib=TRUE,approxcontrib=FALSE));bias<-which(colnames(sh0)=="BIAS");if(length(bias)!=1) stop("SHAP基线列异常")
  raw_margin<-as.numeric(predict(booster,dm,outputmargin=TRUE));raw_error<-max(abs(rowSums(sh0)-raw_margin));if(raw_error>1e-4) stop("SHAP原始加和校验未通过")
  sh<-sh0[,setdiff(seq_len(ncol(sh0)),bias),drop=FALSE];sh<-sh[,colnames(x),drop=FALSE];ca<-state$calibrators$xgboost
  # 极端概率在校准函数中会裁剪；仅在线性区间解释Platt变换，否则解释原始模型并明确标记。
  raw_prob<-predict_prob(state$models$xgboost,x);linear_cal<-use_platt&&all(raw_prob>1e-6&raw_prob<1-1e-6)
  if(linear_cal){base<-ca[1]+ca[2]*sh0[,bias];sh<-sh*ca[2];explained_margin<-ca[1]+ca[2]*raw_margin;scale_label<-"calibrated model log-odds"
  } else {base<-sh0[,bias];explained_margin<-raw_margin;scale_label<-"raw model log-odds"}
  pm<-z$parent_map[colnames(x)];if(anyNA(pm)) stop("SHAP设计列无法映射到原始变量")
  group_shap<-sapply(selected,function(v) rowSums(sh[,pm==v,drop=FALSE]));group_shap<-as.matrix(group_shap)
  grouped_error<-max(abs(rowSums(group_shap)+base-explained_margin));if(grouped_error>1e-4) stop("SHAP原始变量汇总加和校验未通过")
  imp<-data.frame(Variable=selected,MeanAbsSHAP=colMeans(abs(group_shap)));imp<-imp[order(-imp$MeanAbsSHAP),];imp$Rank<-seq_len(nrow(imp));rownames(imp)<-NULL
  shap_results[[co]]<-list(importance=imp,values=group_shap,base=base,stay_id=ids[[co]],data=z$data[[co]],scale_label=scale_label)
  shap_checks[[co]]<-data.frame(Cohort=co,Raw_additivity_error=raw_error,Grouped_additivity_error=grouped_error,Scale=scale_label)
}

shap_table<-do.call(rbind,lapply(names(shap_results),function(co) data.frame(Cohort=co,shap_results[[co]]$importance)))
marker_table<-do.call(rbind,lapply(names(shap_results),function(co) do.call(rbind,lapply(c("GV","SHR","HGI"),function(v){a<-shap_results[[co]]$importance;i<-match(v,a$Variable)
  data.frame(Cohort=co,Variable=v,Included=v%in%selected,Rank=if(is.na(i)) NA_integer_ else a$Rank[i],MeanAbsSHAP=if(is.na(i)) NA_real_ else a$MeanAbsSHAP[i])}))))
print(marker_table,row.names=FALSE)

# ==================== 12. 汇总保存：一个Excel、一个模型RDS ====================
tuning_table<-do.call(rbind,lapply(names(labels),function(id) do.call(rbind,lapply(seq_along(state$tuning[[id]]),function(i){a<-state$tuning[[id]][[i]]
  data.frame(Model=id,Trial=i,CV_AP=a$mean_ap,CV_AUC=a$mean_auc,Parameters=paste(names(a$params),unlist(a$params),sep="=",collapse="; "))}))))
parameter_table<-do.call(rbind,lapply(names(labels),function(id) data.frame(Model=id,Parameters=paste(names(state$best_params[[id]]),unlist(state$best_params[[id]]),sep="=",collapse="; "),Platt_intercept=state$calibrators[[id]][1],Platt_slope=state$calibrators[[id]][2],Threshold=state$thresholds[[id]])))
sheets<-list("性能与置信区间"=performance,"训练CV选择"=cv_summary,"最终变量"=data.frame(Variable=selected),"最优参数"=parameter_table,"调参记录"=tuning_table,
  "分折筛选频次"=stability,"SHAP排名"=shap_table,"血糖指标SHAP"=marker_table,"SHAP加和校验"=do.call(rbind,shap_checks))
for(n in names(sheets))write.csv(sheets[[n]],file.path(out_dir,"results",paste0("机器学习_",n,".csv")),row.names=FALSE,fileEncoding="UTF-8")
write.csv(predictions,file.path(out_dir,"private_model_inputs/固定模型预测.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(fold_selection,file.path(out_dir,"results/每折完整筛选记录.csv"),row.names=FALSE,fileEncoding="UTF-8")
capture.output(sessionInfo(),file=file.path(out_dir,"audit/sessionInfo.txt"))
state$performance<-performance;state$predictions_table<-predictions;state$plot_sources<-plot_sources;state$shap<-shap_results;state$shap_checks<-shap_checks;state$marker_table<-marker_table;state$session<-sessionInfo();state$completed<-TRUE;save_state()
cat("\n完成：",if(quick_test) "代码快速测试（不是正式分析）" else "正式机器学习与SHAP分析","\n结果目录：",out_dir,"\n")
cat("按训练CV AP选择：",labels[best_model],"；SHAP预先指定解释XGBoost。\n")
