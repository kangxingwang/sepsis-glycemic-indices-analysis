if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent=TRUE)
options(stringsAsFactors=FALSE,warn=1,width=180)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggbeeswarm);library(mlr3learners)})
out_dir<-Sys.getenv("SEPSIS_ML_OUTPUT")
state<-readRDS(file.path(out_dir,'private_checkpoints/机器学习模型.rds'))
z<-readRDS(file.path(out_dir,'private_checkpoints/特征工程结果.rds'))
fe<-readRDS(file.path(out_dir,'private_checkpoints/特征工程图表对象.rds'))
stopifnot(isTRUE(state$completed),!state$key$cfg$quick_test)
labels<-c(log_reg='Logistic regression',xgboost='XGBoost',ranger='Random forest',svm='SVM-RBF',kknn='k-NN')
cols<-c('Logistic regression'='#356B9A','XGBoost'='#318C83','Random forest'='#D08A36','SVM-RBF'='#8D69A3','k-NN'='#59636F')
cohort_labels<-c(internal='MIMIC-IV internal validation',nwicu='NWICU external validation')
pred<-state$predictions_table
clip<-function(p)pmin(pmax(as.numeric(p),1e-6),1-1e-6)
# Finalize the prespecified Youden rule from locked training OOF scores only.
# An actual returned maximizer is selected; averaging two separated maxima is prohibited.
threshold_audit<-list()
safe_div<-function(a,b)if(b>0)a/b else NA_real_
for(id in names(labels)){
  y<-as.integer(z$raw$train$Outcome)
  ca<-state$calibrators[[id]]
  op<-plogis(ca[1]+ca[2]*qlogis(clip(state$oof[[id]])))
  rr<-pROC::roc(y,op,levels=c(0,1),direction='<',quiet=TRUE)
  vv<-pROC::coords(rr,'best',best.method='youden',ret='threshold',transpose=FALSE)$threshold
  vv<-sort(vv[is.finite(vv)]);stopifnot(length(vv)>0)
  cut<-vv[ceiling(length(vv)/2)]
  all_coords<-pROC::coords(rr,'all',ret=c('sensitivity','specificity'),transpose=FALSE)
  max_j<-max(all_coords$sensitivity+all_coords$specificity-1)
  old_cut<-state$thresholds[[id]]
  actual_j<-mean(op[y==1]>=cut)+mean(op[y==0]<cut)-1
  old_j<-mean(op[y==1]>=old_cut)+mean(op[y==0]<old_cut)-1
  stopifnot(abs(actual_j-max_j)<1e-12,cut%in%vv)
  threshold_audit[[id]]<-data.frame(Model=id,Training_N=length(y),Training_events=sum(y),Number_best_thresholds=length(vv),Prior_threshold=old_cut,Final_threshold=cut,Prior_Youden_J=old_j,Final_Youden_J=actual_j,Maximum_Youden_J=max_j,Actual_maximizer=TRUE,Threshold_changed=old_cut!=cut,Rule='sort(best_finite_thresholds)[ceiling(length/2)]')
  state$thresholds[[id]]<-cut
  state$predictions[[id]]$Threshold<-cut
  pred$Threshold[pred$Model==id]<-cut
  for(co in names(cohort_labels)){
    d<-pred[pred$Model==id&pred$Cohort==co,];cl<-d$Probability>=cut;yy<-d$Outcome
    tp<-sum(cl&yy==1);tn<-sum(!cl&yy==0);fp<-sum(cl&yy==0);fn<-sum(!cl&yy==1)
    index<-which(state$performance$Model==id&state$performance$Cohort==co)
    state$performance[index,c('Threshold','Sensitivity','Specificity','PPV','NPV','F1','Accuracy')]<-list(cut,safe_div(tp,tp+fn),safe_div(tn,tn+fp),safe_div(tp,tp+fp),safe_div(tn,tn+fn),safe_div(2*tp,2*tp+fp+fn),mean(cl==yy))
  }
}
state$predictions_table<-pred
state$threshold_rule<-'Actual central maximizer of training-OOF Youden J; for an even number, take the lower central returned threshold.'
temp_state<-file.path(out_dir,'private_checkpoints/机器学习模型.rds.tmp')
saveRDS(state,temp_state)
unlink(file.path(out_dir,'private_checkpoints/机器学习模型.rds'))
stopifnot(file.rename(temp_state,file.path(out_dir,'private_checkpoints/机器学习模型.rds')))
write.csv(state$performance,file.path(out_dir,'results/机器学习_性能与置信区间.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(pred,file.path(out_dir,'private_model_inputs/固定模型预测.csv'),row.names=FALSE,fileEncoding='UTF-8')
parameter_file<-file.path(out_dir,'results/机器学习_最优参数.csv')
parameter_table<-read.csv(parameter_file,fileEncoding='UTF-8-BOM')
parameter_table$Threshold<-vapply(parameter_table$Model,function(id)state$thresholds[[id]],numeric(1))
write.csv(parameter_table,parameter_file,row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,threshold_audit),file.path(out_dir,'audit/OOF_Youden实际最优阈值核验.csv'),row.names=FALSE,fileEncoding='UTF-8')
first_threshold_record<-file.path(out_dir,'audit/OOF_Youden实际最优阈值核验_首次.csv')
if(!file.exists(first_threshold_record))write.csv(do.call(rbind,threshold_audit),first_threshold_record,row.names=FALSE,fileEncoding='UTF-8')
oof_fold<-integer(nrow(z$raw$train));for(k in seq_along(state$folds))oof_fold[state$folds[[k]]$iv]<-k
oof_table<-do.call(rbind,lapply(names(labels),function(id){
  ca<-state$calibrators[[id]];raw<-as.numeric(state$oof[[id]])
  data.frame(stay_id=z$raw$train$stay_id,Training_row=seq_len(nrow(z$raw$train)),Fold=oof_fold,Model=id,Outcome=z$raw$train$Outcome,Raw_probability=raw,Probability=plogis(ca[1]+ca[2]*qlogis(clip(raw))))
}))
write.csv(oof_table,file.path(out_dir,'private_model_inputs/训练OOF预测_锁定.csv'),row.names=FALSE,fileEncoding='UTF-8')
fold_scores<-do.call(rbind,lapply(names(labels),function(id)do.call(rbind,lapply(seq_along(state$tuning[[id]]),function(j){
  a<-state$tuning[[id]][[j]]
  data.frame(Model=id,Trial=j,Fold=seq_len(10L),AP=a$fold_ap,AUC=a$fold_auc)
}))))
write.csv(fold_scores,file.path(out_dir,'audit/全部调参逐折AP与AUC.csv'),row.names=FALSE,fileEncoding='UTF-8')
theme_fig<-function(){theme_classic(base_family='Arial',base_size=7)+theme(axis.line=element_line(linewidth=.35,colour='black'),axis.ticks=element_line(linewidth=.35),axis.ticks.length=grid::unit(1.4,'mm'),axis.text=element_text(size=6.5,colour='black'),axis.title=element_text(size=7),legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),plot.title=element_text(size=8,face='bold',hjust=0),strip.text=element_text(size=7,face='bold'),plot.tag=element_text(size=9,face='bold'),legend.position='bottom',plot.margin=grid::unit(c(3,3,3,3),'mm'))}
theme_set(theme_fig())
figure_manifest<-list()
save_fig<-function(p,name,width=183,height=125,source=''){
  file<-file.path(out_dir,'figures',paste0(name,'.png'))
  ragg::agg_png(file,width=width,height=height,units='mm',res=300,background='white')
  print(p);dev.off()
  figure_manifest[[length(figure_manifest)+1]]<<-data.frame(Figure=name,Width_mm=width,Height_mm=height,PNG_DPI=300,Backend='R',Source=source)
}

write_result<-function(d,name)write.csv(d,file.path(out_dir,'results',paste0(name,'.csv')),row.names=FALSE,fileEncoding='UTF-8',na='')
write_source<-function(d,name)write.csv(d,file.path(out_dir,'private_plot_data',paste0(name,'.csv')),row.names=FALSE,fileEncoding='UTF-8',na='')

# Independent checks rebuild preprocessing parameters using only the designated training rows.
checks<-list()
check<-function(name,pass,detail=''){checks[[length(checks)+1]]<<-data.frame(Check=name,Pass=isTRUE(pass),Detail=detail);if(!isTRUE(pass))stop('Audit failed: ',name)}
check('Original source files unchanged',identical(unname(z$source_md5),unname(tools::md5sum(z$paths))))
check('Training checkpoint uses this feature-engineering object',identical(state$key$input_md5,unname(tools::md5sum(file.path(out_dir,'private_checkpoints/特征工程结果.rds')))))
source_data<-lapply(z$paths,function(p)read.csv(p,fileEncoding='UTF-8-BOM',check.names=FALSE,colClasses=c(stay_id='character',subject_id='character',hadm_id='character')))
for(co in c('train','internal','nwicu')){
  d<-z$raw[[co]];source<-source_data[[if(co=='nwicu')'NWICU'else'MIMICIV']]
  index<-match(d$stay_id,source$stay_id)
  check(paste(co,'original IDs and 28-day labels'),!anyNA(index)&&all(d$subject_id==source$subject_id[index])&&all(d$hadm_id==source$hadm_id[index])&&all(d$Outcome==source$day28_outcome[index]))
  write.csv(z$data[[co]],file.path(out_dir,'private_model_inputs',paste0(co,'_最终预处理输入.csv')),row.names=FALSE,fileEncoding='UTF-8')
}
check('MIMIC unique patient split',!any(z$raw$train$subject_id%in%z$raw$internal$subject_id))
check('Training N/events',nrow(z$raw$train)==3635&&sum(z$raw$train$Outcome)==525)
check('Internal N/events',nrow(z$raw$internal)==1560&&sum(z$raw$internal$Outcome)==226)
check('NWICU N/events',nrow(z$raw$nwicu)==1140&&sum(z$raw$nwicu$Outcome)==260)
check('Formal trials and resampling',identical(as.integer(state$key$cfg$tune_n),c(30L,20L,20L,20L))&&state$key$cfg$cv_folds==10&&state$key$cfg$rf_trees==1500&&state$key$cfg$bootstrap_n==1000)
check('Boruta formal settings',z$settings$boruta_runs==200&&z$settings$boruta_trees==500&&z$settings$boruta_p==.01&&isTRUE(z$settings$boruta_mcAdj))
check('Training-only final HGI',max(abs(coef(lm(HbA1c~Glucose,data=z$raw$train))-z$prep$hgi))<1e-10)
check('No outcome/time/ID as predictors',!any(grepl('Outcome|outcome|death|day28|LOS|los|stay_id|subject_id|hadm_id',z$candidates)))
check('Vital exceptions preserved',all(c('MAP','Heart_Rate','Resp_Rate','Temperature')%in%z$candidates))
check('All other predictors meet missingness rule',all(subset(z$missing_table,Keep&!Vital_exception)$MIMIC_missing_pct<=20)&all(subset(z$missing_table,Keep&!Vital_exception)$NWICU_missing_pct<=20))
check('Only confirmed Boruta intersection',setequal(z$selected,subset(z$selection_table,LASSO&Boruta=='Confirmed')$Variable))
fold_parameters<-list()
validation_coverage<-integer(nrow(z$raw$train))
mode_value<-function(x){a<-table(x);names(a)[which.max(a)]}
audit_prep<-function(d,p,name){
  h<-coef(lm(HbA1c~Glucose,data=d));check(paste(name,'HGI parameters'),max(abs(h-p$hgi))<1e-10)
  d$HGI<-d$HbA1c-(h[1]+h[2]*d$Glucose)
  for(v in names(p$median))check(paste(name,'median',v),abs(median(d[[v]],na.rm=TRUE)-p$median[[v]])<1e-10)
  for(v in names(p$mode))check(paste(name,'mode',v),identical(mode_value(d[[v]]),p$mode[[v]]))
  for(v in names(p$levels))check(paste(name,'training-only categorical levels',v),setequal(unique(na.omit(as.character(d[[v]]))),p$levels[[v]]))
  for(v in names(p$bounds))check(paste(name,'winsor',v),max(abs(quantile(d[[v]],probs=c(.01,.99),na.rm=TRUE,names=FALSE)-p$bounds[[v]]))<1e-10)
}
prepared_matrix<-function(d,p,variables){
  d$HGI<-d$HbA1c-(p$hgi[1]+p$hgi[2]*d$Glucose)
  for(v in names(p$median)){
    d[[v]][is.na(d[[v]])]<-p$median[[v]]
    if(v%in%names(p$bounds))d[[v]]<-pmin(pmax(d[[v]],p$bounds[[v]][1]),p$bounds[[v]][2])
  }
  for(v in names(p$levels)){
    x<-as.character(d[[v]]);x[is.na(x)|!x%in%p$levels[[v]]]<-p$mode[[v]]
    d[[v]]<-factor(x,levels=p$levels[[v]])
  }
  x<-model.matrix(~.,d[,variables,drop=FALSE]);x<-x[,colnames(x)!='(Intercept)',drop=FALSE];storage.mode(x)<-'double';x
}
audit_scale<-function(d,p,selected,scale,stored,name,training=TRUE){
  x<-prepared_matrix(d,p,selected);x<-x[,scale$columns,drop=FALSE]
  if(training){
    check(paste(name,'training-only centering'),max(abs(colMeans(x)-scale$center))<1e-10)
    check(paste(name,'training-only scaling'),max(abs(apply(x,2,sd)-scale$scale))<1e-10)
  }
  x<-sweep(sweep(x,2,scale$center,'-'),2,scale$scale,'/')
  check(paste(name,'reconstructed model matrix'),identical(dim(x),dim(stored))&&identical(colnames(x),colnames(stored))&&max(abs(x-stored))<1e-10)
}
audit_prep(z$raw$train,z$prep,'Final')
audit_scale(z$raw$train,z$prep,z$selected,state$scale_final,state$x$train,'Final train')
audit_scale(z$raw$internal,z$prep,z$selected,state$scale_final,state$x$internal,'Final internal',FALSE)
audit_scale(z$raw$nwicu,z$prep,z$selected,state$scale_final,state$x$nwicu,'Final NWICU',FALSE)
for(k in seq_along(state$folds)){
  f<-state$folds[[k]]
  check(paste('Fold',k,'disjoint'),!length(intersect(f$it,f$iv))&&setequal(c(f$it,f$iv),seq_len(nrow(z$raw$train))))
  validation_coverage[f$iv]<-validation_coverage[f$iv]+1L
  audit_prep(z$raw$train[f$it,],f$prep,paste('Fold',k))
  audit_scale(z$raw$train[f$it,],f$prep,f$selected,f$scale,f$x_train,paste('Fold',k,'train'))
  audit_scale(z$raw$train[f$iv,],f$prep,f$selected,f$scale,f$x_valid,paste('Fold',k,'heldout'),FALSE)
  check(paste('Fold',k,'intersection'),all(f$selection$Intersection==(f$selection$LASSO&f$selection$Boruta=='Confirmed')))
  fold_parameters[[k]]<-data.frame(Fold=k,Training_N=length(f$it),Validation_N=length(f$iv),Training_events=sum(z$raw$train$Outcome[f$it]),Validation_events=sum(z$raw$train$Outcome[f$iv]),Selected_variables=length(f$selected),HGI_intercept=unname(f$prep$hgi[1]),HGI_slope=unname(f$prep$hgi[2]))
}
check('Each training patient held out exactly once',all(validation_coverage==1L))
for(id in names(labels)){
  ts<-state$tuning[[id]]
  expected<-if(id=='log_reg')1L else state$key$cfg$tune_n[[id]]
  check(paste(id,'all formal trials complete'),length(ts)==expected&&all(vapply(ts,function(x)!is.null(x)&&length(x$fold_ap)==10&&!anyNA(x$oof),logical(1))))
  best<-which.max(vapply(ts,`[[`,numeric(1),'mean_ap'))
  check(paste(id,'AP-only hyperparameter selection'),identical(state$best_params[[id]],ts[[best]]$params)&&identical(state$oof[[id]],ts[[best]]$oof))
  oof<-clip(state$oof[[id]]);fit<-glm(z$raw$train$Outcome~qlogis(oof),family=binomial())
  check(paste(id,'OOF-only calibrator'),max(abs(unname(coef(fit))-unname(state$calibrators[[id]])))<1e-8)
  cal<-plogis(state$calibrators[[id]][1]+state$calibrators[[id]][2]*qlogis(oof))
  rr<-pROC::roc(z$raw$train$Outcome,cal,levels=c(0,1),direction='<',quiet=TRUE)
  cut<-pROC::coords(rr,'best',best.method='youden',ret='threshold',transpose=FALSE)$threshold;cut<-sort(cut[is.finite(cut)]);cut<-cut[ceiling(length(cut)/2)]
  check(paste(id,'OOF-only fixed threshold'),abs(cut-state$thresholds[[id]])<1e-10)
}
write_result(do.call(rbind,fold_parameters),'折内样本及HGI方程')

# Independent metric recalculation and paired DeLong comparisons.
ap<-function(y,p){levels<-sort(unique(p),decreasing=TRUE);tp<-cumsum(vapply(levels,function(t)sum(y[p==t]),numeric(1)));n<-cumsum(vapply(levels,function(t)sum(p==t),integer(1)));sum(diff(c(0,tp/sum(y)))*tp/n)}
metric_check<-list();delong<-list()
for(co in names(cohort_labels)){
  dd<-pred[pred$Cohort==co,]
  for(id in names(labels)){
    d<-dd[dd$Model==id,];m<-subset(state$performance,Cohort==co&Model==id)
    check(paste(co,id,'labels and row IDs'),identical(as.character(d$stay_id),as.character(z$raw[[co]]$stay_id))&&identical(as.integer(d$Outcome),as.integer(z$raw[[co]]$Outcome)))
    raw_reprediction<-as.numeric(state$models[[id]]$predict_newdata(as.data.frame(state$x[[co]],check.names=FALSE))$prob[,'1'])
    check(paste(co,id,'locked RDS model reproduces probabilities'),max(abs(raw_reprediction-d$Raw_probability))<1e-8)
    cc<-state$calibrators[[id]];recalibrated<-plogis(cc[1]+cc[2]*qlogis(clip(raw_reprediction)))
    check(paste(co,id,'same training calibrator applied'),max(abs(recalibrated-d$Probability))<1e-8)
    auc<-as.numeric(pROC::auc(pROC::roc(d$Outcome,d$Probability,levels=c(0,1),direction='<',quiet=TRUE)))
    ap_value<-ap(d$Outcome,d$Probability);brier<-mean((d$Probability-d$Outcome)^2)
    metric_check[[length(metric_check)+1]]<-data.frame(Cohort=co,Model=id,AUC_abs_diff=abs(auc-m$AUC),AP_abs_diff=abs(ap_value-m$AP),Brier_abs_diff=abs(brier-m$Brier))
  }
  for(pair in combn(names(labels),2,simplify=FALSE)){
    a<-dd[dd$Model==pair[1],];b<-dd[dd$Model==pair[2],]
    stopifnot(identical(a$stay_id,b$stay_id),identical(a$Outcome,b$Outcome))
    ra<-pROC::roc(a$Outcome,a$Probability,levels=c(0,1),direction='<',quiet=TRUE)
    rb<-pROC::roc(b$Outcome,b$Probability,levels=c(0,1),direction='<',quiet=TRUE)
    dt<-pROC::roc.test(ra,rb,method='delong',paired=TRUE)
    delong[[length(delong)+1]]<-data.frame(Cohort=co,Model_A=pair[1],Model_B=pair[2],AUC_A=as.numeric(pROC::auc(ra)),AUC_B=as.numeric(pROC::auc(rb)),Delta_AUC=as.numeric(pROC::auc(ra)-pROC::auc(rb)),Z=unname(dt$statistic),P=dt$p.value)
  }
}
metric_check<-do.call(rbind,metric_check);check('Independent performance calculation',max(as.matrix(metric_check[,3:5]))<1e-10)
write_result(metric_check,'性能独立重算差值')
delong<-do.call(rbind,delong);delong$Holm_P<-ave(delong$P,delong$Cohort,FUN=function(p)p.adjust(p,'holm'));write_result(delong,'两两AUC_DeLong及Holm')

# Original-NA training/internal/external baseline; missingness is explicit, not hidden by imputation.
numeric_summary<-function(x){x<-x[!is.na(x)];sprintf('%.2f [%.2f, %.2f]',median(x),quantile(x,.25,names=FALSE),quantile(x,.75,names=FALSE))}
baseline<-list()
for(v in c('Outcome','HbA1c','Glucose',z$candidates)){
  a<-z$raw$train[[v]];b<-z$raw$internal[[v]];c<-z$raw$nwicu[[v]]
  if(v=='HGI'){
    # The baseline uses the training-estimated equation uniformly, like the models.
    a<-z$raw$train$HbA1c-(z$prep$hgi[1]+z$prep$hgi[2]*z$raw$train$Glucose)
    b<-z$raw$internal$HbA1c-(z$prep$hgi[1]+z$prep$hgi[2]*z$raw$internal$Glucose)
    c<-z$raw$nwicu$HbA1c-(z$prep$hgi[1]+z$prep$hgi[2]*z$raw$nwicu$Glucose)
  }
  if(is.numeric(a)&&v!='Outcome'){
    denom<-sqrt((var(a,na.rm=TRUE)+var(b,na.rm=TRUE))/2);smd<-if(denom>0)abs(mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/denom else 0
    test<-wilcox.test(a,b,exact=FALSE)$p.value
    baseline[[length(baseline)+1]]<-data.frame(Variable=v,Level='Median [IQR]',Training=numeric_summary(a),Internal=numeric_summary(b),NWICU=numeric_summary(c),Training_missing=sum(is.na(a)),Internal_missing=sum(is.na(b)),NWICU_missing=sum(is.na(c)),Abs_SMD=smd,P=test,Final_selected=v%in%z$selected)
  }else{
    levels<-sort(unique(na.omit(c(as.character(a),as.character(b),as.character(c)))))
    contingency<-table(Cohort=c(rep('Train',length(a)),rep('Internal',length(b))),Value=c(a,b))
    test<-suppressWarnings(chisq.test(contingency,correct=FALSE)$p.value)
    for(lev in levels){
      pa<-mean(as.character(a)==lev,na.rm=TRUE);pb<-mean(as.character(b)==lev,na.rm=TRUE);denom<-sqrt((pa*(1-pa)+pb*(1-pb))/2);smd<-if(denom>0)abs(pa-pb)/denom else 0
      val<-function(x)sprintf('%d (%.1f%%)',sum(as.character(x)==lev,na.rm=TRUE),100*mean(as.character(x)==lev,na.rm=TRUE))
      baseline[[length(baseline)+1]]<-data.frame(Variable=v,Level=lev,Training=val(a),Internal=val(b),NWICU=val(c),Training_missing=sum(is.na(a)),Internal_missing=sum(is.na(b)),NWICU_missing=sum(is.na(c)),Abs_SMD=smd,P=test,Final_selected=v%in%z$selected)
    }
  }
}
write_result(do.call(rbind,baseline),'训练内部及NWICU基线_原NA')

# SHAP probability reconstruction is independently checked from grouped values.
shap_check<-list()
for(co in names(cohort_labels)){
  sh<-state$shap[[co]];margin<-rowSums(sh$values)+sh$base
  dd<-subset(pred,Cohort==co&Model=='xgboost')
  if(grepl('calibrated',sh$scale_label)){
    p<-plogis(margin)
  }else{
    a<-state$calibrators$xgboost;p<-plogis(a[1]+a[2]*qlogis(clip(plogis(margin))))
  }
  err<-max(abs(p-dd$Probability));check(paste(co,'grouped SHAP predicts calibrated probability'),err<1e-4)
  shap_check[[co]]<-data.frame(Cohort=co,Probability_reconstruction_max_abs_error=err,Scale=sh$scale_label)
}
write_result(do.call(rbind,shap_check),'SHAP预测概率重建核验')

# Training selection figures.
cv<-z$lasso
cv_df<-data.frame(Log_lambda=log(cv$lambda),Deviance=cv$cvm,SE=cv$cvsd,N_nonzero=cv$nzero)
write_source(cv_df,'LASSO_CV')
top_index<-unique(round(seq(1,nrow(cv_df),length.out=6)))
cvplot<-ggplot(cv_df,aes(Log_lambda,Deviance))+geom_errorbar(aes(ymin=Deviance-SE,ymax=Deviance+SE),width=.02,linewidth=.25,colour='#A7AFB8')+geom_point(size=.8,colour='#356B9A')+geom_vline(xintercept=log(cv$lambda.min),linetype='dotted',colour='#59636F',linewidth=.35)+geom_vline(xintercept=log(cv$lambda.1se),linetype='dashed',colour='#D08A36',linewidth=.4)+scale_x_continuous(sec.axis=dup_axis(breaks=cv_df$Log_lambda[top_index],labels=cv_df$N_nonzero[top_index],name='Nonzero encoded coefficients'))+labs(x=expression(log(lambda)),y='10-fold CV binomial deviance',title='LASSO selection in the training cohort',subtitle='Dashed: lambda.1se; dotted: lambda.min')+theme(plot.subtitle=element_text(size=6.5))
save_fig(cvplot,'Selection_LASSO_CV',height=105,source='LASSO_CV.csv')
bm<-as.matrix(cv$glmnet.fit$beta)
path<-data.frame(Feature=rep(rownames(bm),ncol(bm)),Lambda=rep(cv$glmnet.fit$lambda,each=nrow(bm)),Coefficient=as.vector(bm))
path$Variable<-unname(z$parent_map[path$Feature]);path$Glycemic<-ifelse(path$Variable%in%c('GV','SHR','HGI'),path$Variable,'Other')
write_source(path,'LASSO_path')
pathplot<-ggplot(path,aes(log(Lambda),Coefficient,group=Feature,colour=Glycemic))+geom_line(linewidth=.3,alpha=.85)+geom_vline(xintercept=log(cv$lambda.min),linetype='dotted',colour='#59636F',linewidth=.35)+geom_vline(xintercept=log(cv$lambda.1se),linetype='dashed',linewidth=.4)+scale_x_continuous(sec.axis=dup_axis(breaks=cv_df$Log_lambda[top_index],labels=cv_df$N_nonzero[top_index],name='Nonzero encoded coefficients'))+scale_colour_manual(values=c(Other='#BBC1C7',GV='#356B9A',SHR='#D08A36',HGI='#318C83'))+labs(x=expression(log(lambda)),y='Coefficient',colour=NULL,title='LASSO coefficient paths')
save_fig(pathplot+theme(plot.margin=grid::unit(c(7,3,3,3),'mm')),'Selection_LASSO_path',height=105,source='LASSO_path.csv')
bo<-fe$boruta_table;bo$Variable<-factor(bo$Variable,levels=bo$Variable[order(bo$medianImp)])
hist<-z$boruta$ImpHistory
bo$Q25<-vapply(as.character(bo$Variable),function(v)quantile(hist[,v],.25,na.rm=TRUE,names=FALSE),numeric(1))
bo$Q75<-vapply(as.character(bo$Variable),function(v)quantile(hist[,v],.75,na.rm=TRUE,names=FALSE),numeric(1))
write_source(bo,'Boruta_importance')
boxdf<-data.frame(Variable=rep(colnames(hist),each=nrow(hist)),Importance=as.vector(hist))
boxdf<-boxdf[is.finite(boxdf$Importance),]
decision_map<-c(setNames(as.character(z$boruta$finalDecision),names(z$boruta$finalDecision)),shadowMax='Shadow',shadowMean='Shadow',shadowMin='Shadow')
boxdf$Decision<-unname(decision_map[boxdf$Variable]);stopifnot(!anyNA(boxdf$Decision))
orders<-aggregate(Importance~Variable,boxdf,median);orders<-orders$Variable[order(orders$Importance)]
boxdf$Variable<-factor(boxdf$Variable,levels=orders)
write_source(boxdf,'Boruta_importance_history')
borutaplot<-ggplot(boxdf,aes(Importance,Variable,fill=Decision))+geom_boxplot(width=.6,linewidth=.25,outlier.size=.3,outlier.alpha=.35)+scale_fill_manual(values=c(Confirmed='#318C83',Tentative='#D08A36',Rejected='#B8BEC4',Shadow='#356B9A'))+labs(x='Boruta permutation importance Z score',y=NULL,fill='Decision',title='Boruta selection in the training cohort')+theme(axis.text.y=element_text(size=6.5))
save_fig(borutaplot,'Selection_Boruta',height=210,source='Boruta_importance_history.csv;Boruta_importance.csv')
selection_counts<-data.frame(Set=c('LASSO only','Intersection','Boruta only'),N=c(sum(z$selection_table$LASSO&z$selection_table$Boruta!='Confirmed'),length(z$selected),sum(!z$selection_table$LASSO&z$selection_table$Boruta=='Confirmed')))
write_source(selection_counts,'Selection_intersection')
intplot<-ggplot(selection_counts,aes(Set,N,fill=Set))+geom_col(width=.6)+geom_text(aes(label=N),vjust=-.35,size=2.6)+scale_fill_manual(values=c('LASSO only'='#356B9A',Intersection='#318C83','Boruta only'='#D08A36'))+scale_y_continuous(expand=expansion(mult=c(0,.12)))+labs(x=NULL,y='Number of clinical variables',title='LASSO lambda.1se ∩ Boruta Confirmed')+theme(legend.position='none')
save_fig(intplot,'Selection_intersection',height=90,source='Selection_intersection.csv')
stability<-data.frame(Variable=z$candidates,Selected_folds=vapply(z$candidates,function(v)sum(state$fold_selection$Variable==v&state$fold_selection$Intersection),integer(1)),Total_folds=10L)
write_source(stability,'Fold_selection_stability')
stability$Variable<-factor(stability$Variable,levels=stability$Variable[order(stability$Selected_folds,stability$Variable)])
stabplot<-ggplot(stability,aes(Selected_folds,Variable))+geom_col(fill='#356B9A',width=.65)+scale_x_continuous(breaks=0:10,limits=c(0,10),expand=expansion(mult=c(0,.02)))+labs(x='Selected in training folds (out of 10)',y=NULL,title='Selection stability across training folds')+theme(axis.text.y=element_text(size=6.5))
save_fig(stabplot,'Selection_fold_stability',height=200,source='Fold_selection_stability.csv')
fig4<-(cvplot|pathplot)/(borutaplot|intplot)+plot_layout(heights=c(1,2))+plot_annotation(tag_levels='a')
save_fig(fig4,'Figure4_LASSO_Boruta',height=230,source='LASSO_CV.csv;LASSO_path.csv;Boruta_importance_history.csv;Selection_intersection.csv')

# Fixed model ROC, calibration, and DCA. All curves retain the same color mapping.
plot_list<-list();roc_all<-cal_all<-dca_all<-list()
for(co in names(cohort_labels)){
  dd<-pred[pred$Cohort==co,]
  co_n<-nrow(z$raw[[co]]);co_y<-sum(z$raw[[co]]$Outcome)
  title<-sprintf('%s (n = %s; deaths = %s)',cohort_labels[co],format(co_n,big.mark=','),co_y)
  rocdf<-do.call(rbind,lapply(names(labels),function(id){d<-dd[dd$Model==id,];r<-pROC::roc(d$Outcome,d$Probability,levels=c(0,1),direction='<',quiet=TRUE);data.frame(Cohort=co,Model=unname(labels[id]),FPR=1-r$specificities,Sensitivity=r$sensitivities)}))
  mt<-state$performance[state$performance$Cohort==co,];mt<-mt[match(names(labels),mt$Model),]
  legend_labels<-setNames(sprintf('%s (AUC %.3f)',unname(labels[mt$Model]),mt$AUC),unname(labels[mt$Model]))
  rocp<-ggplot(rocdf,aes(FPR,Sensitivity,colour=Model))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='#A7AFB8',linewidth=.3)+geom_path(linewidth=.5)+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+scale_x_continuous(breaks=seq(0,1,.2))+scale_y_continuous(breaks=seq(0,1,.2))+scale_colour_manual(values=cols,breaks=unname(labels),labels=legend_labels)+labs(x='False-positive rate',y='Sensitivity',colour=NULL,title=title)+guides(colour=guide_legend(ncol=2))
  save_fig(rocp,paste0(co,'_ROC'),height=130,source=paste0(co,'_ROC.csv'))
  caldf<-do.call(rbind,lapply(names(labels),function(id){d<-dd[dd$Model==id,];br<-unique(quantile(d$Probability,seq(0,1,.1),names=FALSE));bin<-if(length(br)<2)rep(1,nrow(d))else as.integer(cut(d$Probability,br,include.lowest=TRUE));do.call(rbind,lapply(split(seq_len(nrow(d)),bin),function(ii)data.frame(Cohort=co,Model=unname(labels[id]),N=length(ii),Predicted=mean(d$Probability[ii]),Observed=mean(d$Outcome[ii]))))}))
  calp<-ggplot(caldf,aes(Predicted,Observed,colour=Model))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='#A7AFB8',linewidth=.3)+geom_line(linewidth=.45)+geom_point(size=1)+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+scale_x_continuous(breaks=seq(0,1,.2))+scale_y_continuous(breaks=seq(0,1,.2))+scale_colour_manual(values=cols,breaks=unname(labels))+labs(x='Mean predicted probability',y='Observed mortality proportion',colour=NULL,title=title,subtitle='Deciles of predicted risk; calibration fixed from training OOF')+theme(plot.subtitle=element_text(size=6.2))+guides(colour=guide_legend(ncol=3))
  save_fig(calp,paste0(co,'_Calibration'),height=130,source=paste0(co,'_Calibration.csv'))
  thresholds<-seq(.01,.50,.01)
  dcadf<-do.call(rbind,lapply(names(labels),function(id){d<-dd[dd$Model==id,];data.frame(Cohort=co,Strategy=unname(labels[id]),Threshold=thresholds,NetBenefit=vapply(thresholds,function(t){positive<-d$Probability>=t;mean(positive&d$Outcome==1)-mean(positive&d$Outcome==0)*t/(1-t)},numeric(1)))}))
  prevalence<-mean(z$raw[[co]]$Outcome)
  dcadf<-rbind(dcadf,data.frame(Cohort=co,Strategy='Treat all',Threshold=thresholds,NetBenefit=prevalence-(1-prevalence)*thresholds/(1-thresholds)),data.frame(Cohort=co,Strategy='Treat none',Threshold=thresholds,NetBenefit=0))
  dccols<-c(cols,'Treat all'='#6A6A6A','Treat none'='#AAAAAA');lty<-setNames(c(rep('solid',5),'dashed','dotted'),names(dccols))
  dcap<-ggplot(dcadf,aes(Threshold,NetBenefit,colour=Strategy,linetype=Strategy))+geom_line(linewidth=.5)+scale_colour_manual(values=dccols,breaks=names(dccols))+scale_linetype_manual(values=lty,breaks=names(dccols))+scale_x_continuous(labels=scales::percent_format(accuracy=1))+coord_cartesian(ylim=c(-.02,max(dcadf$NetBenefit)+.01))+labs(x='Threshold probability',y='Net benefit',colour=NULL,linetype=NULL,title=title)+guides(colour=guide_legend(ncol=3),linetype=guide_legend(ncol=3))
  save_fig(dcap,paste0(co,'_DCA'),height=125,source=paste0(co,'_DCA.csv'))
  write_source(rocdf,paste0(co,'_ROC'));write_source(caldf,paste0(co,'_Calibration'));write_source(dcadf,paste0(co,'_DCA'))
  plot_list[[co]]<-list(ROC=rocp,Calibration=calp,DCA=dcap)
}

# Exact XGBoost SHAP, grouped by original variable before ranking.
all_categories<-sort(c(if('Race'%in%z$selected)paste('Race',z$prep$levels$Race,sep=': '),if('ICU_Type'%in%z$selected)paste('ICU_Type',z$prep$levels$ICU_Type,sep=': ')))
all_catcols<-setNames(rep(c('#356B9A','#D08A36','#8D69A3','#59636F','#318C83','#AE7B62','#7F8F55','#9C6B87','#626FA8'),length.out=length(all_categories)),all_categories)
for(co in names(cohort_labels)){
  a<-state$shap[[co]];top<-head(a$importance$Variable,30)
  imp<-a$importance;imp$Variable<-factor(imp$Variable,levels=rev(top))
  shapbar<-ggplot(imp,aes(MeanAbsSHAP,Variable))+geom_col(fill=cols[['XGBoost']],width=.65)+labs(x='Mean |SHAP| (log-odds)',y=NULL,title=paste('XGBoost:',cohort_labels[co]))
  write_source(imp,paste0(co,'_SHAP_importance'))
  save_fig(shapbar,paste0(co,'_SHAP_importance'),height=max(130,6*length(top)+35),source=paste0(co,'_SHAP_importance.csv'))
  percent<-function(x){if(is.numeric(x))return((rank(x,ties.method='average')-.5)/length(x));lev<-if(is.factor(x))levels(x)else unique(as.character(x));if(setequal(lev,c('No','Yes')))return(as.numeric(as.character(x)=='Yes'));if(setequal(lev,c('F','M')))return(as.numeric(as.character(x)=='M'));rep(NA_real_,length(x))}
  long<-do.call(rbind,lapply(top,function(v){x<-a$data[[v]];data.frame(Variable=v,SHAP=a$values[,v],Feature_percentile=percent(x),Category=if(v%in%c('Race','ICU_Type'))paste(v,as.character(x),sep=': ')else NA_character_)}))
  long$Variable<-factor(long$Variable,levels=rev(top))
  feature_col_labels<-if('Gender'%in%top)c('Low / No / Female','High / Yes / Male')else c('Low / No','High / Yes')
  categories<-all_categories;catcols<-all_catcols
  catcounts<-table(factor(long$Category,levels=categories));catlabels<-setNames(paste0(sub('ICU_Type: ','',categories),' (n=',as.integer(catcounts),')'),categories)
  shapp<-ggplot(long,aes(SHAP,Variable,group=Variable))+geom_vline(xintercept=0,colour='#B7B7B7',linewidth=.3)+ggbeeswarm::geom_quasirandom(data=long[is.na(long$Category),],aes(colour=Feature_percentile),orientation='y',width=.35,bandwidth=.5,size=.45,alpha=.7)+scale_colour_gradient2(low='#356B9A',mid='#C0C7CC',high='#D08A36',midpoint=.5,limits=c(0,1),breaks=c(0,1),labels=feature_col_labels,name='Feature value')+labs(x=paste0('SHAP value (',a$scale_label,')'),y=NULL,title=paste('XGBoost:',cohort_labels[co]))+theme(axis.text.y=element_text(size=6.5),legend.position='right',legend.box='vertical',legend.title=element_text(size=6.5),legend.text=element_text(size=6.5))
  if(length(categories))shapp<-shapp+ggbeeswarm::geom_quasirandom(data=long[!is.na(long$Category),],aes(fill=Category),orientation='y',width=.35,bandwidth=.5,shape=21,colour='white',stroke=.08,size=.65,alpha=.8,show.legend=TRUE)+scale_fill_manual(values=catcols,breaks=categories,limits=categories,labels=catlabels,drop=FALSE,name='Clinical category')
  shapp<-shapp+scale_y_discrete(limits=rev(top),drop=FALSE)+guides(colour=guide_colourbar(order=1,barheight=grid::unit(20,'mm'),barwidth=grid::unit(3,'mm')),fill=guide_legend(order=2,override.aes=list(size=1.5,alpha=1)))
  write_source(long,paste0(co,'_SHAP_summary'))
  save_fig(shapp,paste0(co,'_SHAP_summary'),width=183,height=max(150,6*length(top)+35),source=paste0(co,'_SHAP_summary.csv'))
  plot_list[[co]]$SHAP<-shapp
  plot_list[[co]]$SHAPbar<-shapbar
  plot_list[[co]]$SHAPcompact<-shapp+labs(fill='Category (n)')+theme(legend.position='bottom',legend.box='vertical',legend.box.just='center',legend.text=element_text(size=6.5),legend.title=element_text(size=6.5,hjust=.5),legend.spacing.y=grid::unit(.5,'mm'),legend.key.height=grid::unit(2.5,'mm'),legend.key.width=grid::unit(2.5,'mm'),axis.text.y=element_text(size=6.5),plot.margin=grid::unit(c(2,2,2,2),'mm'))+guides(colour=guide_colourbar(order=1,title.position='top',direction='horizontal',barwidth=grid::unit(18,'mm'),barheight=grid::unit(2,'mm')),fill=guide_legend(order=2,title.position='top',ncol=2,byrow=TRUE,override.aes=list(size=1,alpha=1)))
}
# Compact manuscript composites retain all five curves and all selected SHAP variables.
compact<-function(p,co)p+labs(title=cohort_labels[co],subtitle=sprintf('n = %s; deaths = %s',format(nrow(z$raw[[co]]),big.mark=','),sum(z$raw[[co]]$Outcome)))+theme(plot.title=element_text(size=7),plot.subtitle=element_text(size=6.5))
short_labels<-c('Logistic regression'='LR','XGBoost'='XGB','Random forest'='RF','SVM-RBF'='SVM-RBF','k-NN'='kNN')
compact_roc<-function(co){
  m<-state$performance[state$performance$Cohort==co,];m<-m[match(names(labels),m$Model),]
  ll<-setNames(sprintf('%s %.3f',short_labels[unname(labels[m$Model])],m$AUC),unname(labels[m$Model]))
  compact(plot_list[[co]]$ROC,co)+scale_colour_manual(values=cols,breaks=unname(labels),labels=ll)+labs(colour='Model (AUC)')+guides(colour=guide_legend(ncol=3,byrow=TRUE,title.position='top'))
}
compact_cal<-function(co)compact(plot_list[[co]]$Calibration,co)+scale_colour_manual(values=cols,breaks=unname(labels),labels=short_labels)+labs(colour='Model')+guides(colour=guide_legend(ncol=3,byrow=TRUE,title.position='top'))
compact_dca<-function(co){
  vv<-c(cols,'Treat all'='#6A6A6A','Treat none'='#AAAAAA');lt<-setNames(c(rep('solid',5),'dashed','dotted'),names(vv));ll<-c(short_labels,'Treat all'='Treat all','Treat none'='Treat none')
  compact(plot_list[[co]]$DCA,co)+scale_colour_manual(values=vv,breaks=names(vv),labels=ll)+scale_linetype_manual(values=lt,breaks=names(vv),labels=ll)+labs(colour='Strategy',linetype='Strategy')+guides(colour=guide_legend(ncol=3,byrow=TRUE,title.position='top'),linetype=guide_legend(ncol=3,byrow=TRUE,title.position='top'))
}
compact_shap<-function(co)compact(plot_list[[co]]$SHAPcompact,co)+labs(title=paste('XGBoost:',cohort_labels[co]))
fig5<-(compact_roc('internal')|compact_roc('nwicu'))/(compact_shap('internal')|compact_shap('nwicu'))+plot_layout(heights=c(1,1.6))+plot_annotation(tag_levels='a')
save_fig(fig5,'Figure5_ROC_SHAP',width=183,height=247,source='Individual ROC and SHAP source CSVs')
figs5<-(compact_cal('internal')|compact_cal('nwicu'))/(compact_dca('internal')|compact_dca('nwicu'))+plot_annotation(tag_levels='a')
save_fig(figs5,'FigureS5_Calibration_DCA',width=183,height=190,source='Individual calibration and DCA source CSVs')
write_result(do.call(rbind,figure_manifest),'图形导出清单')
csv_files<-list.files(file.path(out_dir,'results'),pattern='[.]csv$',full.names=TRUE)
csv_files<-csv_files[basename(csv_files)!='结果CSV索引.csv']
csv_index<-do.call(rbind,lapply(csv_files,function(p){d<-read.csv(p,fileEncoding='UTF-8-BOM',check.names=FALSE);data.frame(File=basename(p),Rows=nrow(d),Columns=ncol(d),MD5=unname(tools::md5sum(p)))}))
write_result(csv_index,'结果CSV索引')
legends<-c(
  'Figure 4. Training-set feature selection. All selection uses the MIMIC-IV training cohort (n = 3,635; deaths = 525). LASSO uses 10-fold binomial deviance and lambda.1se; the CV panel shows mean deviance ±1 standard error. The dotted and dashed lines identify lambda.min and lambda.1se, respectively; upper-axis numbers count nonzero encoded coefficients. Boruta uses permutation importance (500 trees per iteration; maxRuns = 200; pValue = .01 with multiplicity adjustment); boxes show median and IQR of the recorded importance history, whiskers extend to 1.5 IQR, and individual outliers are displayed. Shadow controls are retained. Only Confirmed variables enter the intersection with LASSO. These boxes describe algorithm iterations and are not clinical effect confidence intervals.',
  'Figure 5. Fixed-model discrimination and XGBoost interpretation. Models are trained using the MIMIC-IV training cohort, with preprocessing and feature selection repeated within each of 10 training CV folds. Hyperparameters are selected by mean fold AP. Training OOF Platt calibration and the OOF Youden threshold are fixed before application to the MIMIC-IV internal cohort (n = 1,560; deaths = 226) and NWICU (n = 1,140; deaths = 260). ROC panels show all five algorithms. SHAP summary panels show exact patient-level TreeSHAP contributions after summing dummy-variable contributions to their original variable; rows are ranked by mean absolute SHAP. Continuous/binary feature colors run from low/no to high/yes; categorical colors retain race and ICU-type counts. Standalone mean-absolute importance bars are also exported. SHAP describes model predictions, not causal effects. AUC uncertainty is reported in the performance table using DeLong intervals.',
  paste('SHAP scales:',paste(names(state$shap),vapply(state$shap,`[[`,character(1),'scale_label'),sep=' = ',collapse='; ')),
  'Figure S5. Fixed-model calibration and decision curves. Calibration points show mean predicted probability and observed mortality within prediction deciles, with no uncertainty intervals plotted. Calibration intercepts/slopes are diagnostic estimates in the evaluation cohorts and are not used to refit their prediction models. DCA shows net benefit for fixed predictions and treat-all/treat-none strategies over thresholds .01–.50; the displayed vertical range focuses on net benefit >= -.02. Full source-data values are preserved. These retrospective full-ICU analyses do not establish prospective early prediction utility.',
  'Standalone selection-stability panel: count of 10 training folds in which each original variable satisfied both LASSO lambda.1se and Boruta Confirmed. No selection is performed in the internal or NWICU evaluation data.',
  'Performance uncertainty: AUC uses DeLong 95% intervals. AP and Brier use 1,000 outcome-stratified bootstrap resamples of the fixed evaluation predictions; no model retraining is included. Paired model AUC comparisons use DeLong tests with Holm adjustment of 10 pairwise tests separately within each evaluation cohort.'
  ,'Threshold rule: among all finite thresholds achieving maximum training-OOF Youden J, sort the returned thresholds and select the actual central value (the lower central value when their number is even). The selected threshold is independently verified to attain the global ROC maximum J and is applied unchanged to both evaluation cohorts.'
  ,'Baseline table: original-NA values are summarized before imputation/winsorization; HGI is recomputed using the fixed MIMIC-IV training equation. Continuous values are median [IQR]. Categorical percentages use observed nonmissing denominators and missing counts are separately stated. Training/internal comparisons use Wilcoxon or Pearson chi-square tests descriptively, with per-variable/per-level absolute standardized differences; these P values are not used for selection.'
  ,'SHAP color detail: continuous features are colored by within-evaluation-cohort percentile ranks (average ranks for ties); No/Yes factors map to 0/1. Race and ICU type use fixed common colors and cohort-specific counts. These colors aid within-cohort reading and do not set common absolute clinical cutoffs across databases.'
)
writeLines(legends,file.path(out_dir,'analysis_notes/图注与统计说明.txt'),useBytes=TRUE)
write.csv(do.call(rbind,checks),file.path(out_dir,'audit/训练隔离及数值核查.csv'),row.names=FALSE,fileEncoding='UTF-8')
digest<-data.frame(File=c(file.path(out_dir,'private_checkpoints/特征工程结果.rds'),file.path(out_dir,'private_checkpoints/机器学习模型.rds')),MD5=unname(tools::md5sum(c(file.path(out_dir,'private_checkpoints/特征工程结果.rds'),file.path(out_dir,'private_checkpoints/机器学习模型.rds')))))
write.csv(digest,file.path(out_dir,'audit/最终模型文件指纹.csv'),row.names=FALSE,fileEncoding='UTF-8')
txt<-c('原五算法机器学习更新完成。','',capture.output(print(z$cohort_table,row.names=FALSE)),paste0('候选变量：',length(z$candidates),'；最终LASSO lambda.1se与Boruta Confirmed交集：',length(z$selected),'项。'),paste('最终变量：',paste(z$selected,collapse=', ')),sprintf('训练HGI = HbA1c - (%.9f + %.9f × Glucose)',z$prep$hgi[1],z$prep$hgi[2]),'',paste('按训练10折平均AP选择算法：',unname(labels[state$best_model])),capture.output(print(state$cv_summary,row.names=FALSE)),'',capture.output(print(state$performance[,c('Cohort','Model','N','Deaths','AUC','AUC_lower','AUC_upper','AP','AP_lower','AP_upper','Brier','Calibration_intercept','Calibration_slope')],row.names=FALSE)),'','正式配置：10折；XGB30、RF20、SVM20、kNN20组；RF1500树；Boruta maxRuns200×500树；1000次分层bootstrap。','预处理和特征筛选均由各训练部分估计，固定模型评估内部和NW；4生命体征按预先例外保留，缺失由训练中位数填补。','训练OOF校准及阈值固定；不能根据NW表现重新选择超参、阈值或算法。','AUC区间为DeLong；AP/Brier区间基于固定预测的分层bootstrap，未覆盖完整训练不确定性。','全ICU GV包含病程信息，结果为回顾性区分，不能宣称入ICU时即可获得的早期预测效能。','SHAP解释XGB预测，不代表因果或独立生物学效应。','PNG figures exported; visual review is recorded separately.','完整CSV在results，正式模型及逐折断点在private_checkpoints。')
writeLines(txt,file.path(out_dir,'先看机器学习结果.txt'),useBytes=TRUE)
saveRDS(list(figures=plot_list,figure_manifest=figure_manifest,checks=checks,delong=delong,baseline=baseline),file.path(out_dir,'private_checkpoints/图表与核查对象.rds'))
ready<-list(best_model=unname(state$best_model),seed=state$key$cfg$seed,bootstrap_n=state$key$cfg$bootstrap_n,quick_test=state$key$cfg$quick_test,completed=isTRUE(state$completed),postprocessing_completed=TRUE,Models=5L,Training_N=nrow(z$raw$train),Training_events=sum(z$raw$train$Outcome),CV_folds=state$key$cfg$cv_folds,Tuning_trials=sum(vapply(state$tuning,length,integer(1))),OOF_rows=nrow(oof_table),All_audit_checks_passed=all(vapply(checks,function(x)x$Pass,logical(1))),Time=as.character(Sys.time()),Ready_for_independent_performance_audit=TRUE,Final_visual_review_record='Final visual inspection is recorded separately')
jsonlite::write_json(ready,file.path(out_dir,'audit/正式后处理完成.json'),pretty=TRUE,auto_unbox=TRUE)
cat('FORMAL_ML_FIGURES_AND_AUDIT_COMPLETED\n')
