if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent=TRUE)
Sys.setenv(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
# MIMIC-IV建模 / NWICU外部验证：特征工程（从第1节开始，每次选中完整一节运行）
# 两库共有且缺失率均≤20% → MIMIC-IV分层70:30 → 训练内插补 → LASSO ∩ Boruta。
# 本文件只做特征工程；只有stay_id作为ID。改变前面数据/参数后，后续小节顺序重跑。
# 输出：CSV与RDS；图形由后续R-only图表步骤统一生成；筛选本身不等同于预测效能验证。

# ==================== 1. 导入包、路径和参数 ====================
pkgs<-c("data.table","glmnet","Boruta","ranger","ggplot2","ggridges","viridisLite","ragg","openxlsx")
miss<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]; if(length(miss)) stop("Install the required R packages before running: ",paste(miss,collapse=", "))
suppressPackageStartupMessages({library(data.table);library(glmnet);library(Boruta);library(ggplot2);library(ggridges)})
options(width=160); seed<-1117L; train_ratio<-.70; missing_limit<-.20; set.seed(seed)
paths<-c(MIMICIV=Sys.getenv("SEPSIS_ML_MIMIC_CSV"),NWICU=Sys.getenv("SEPSIS_ML_NW_CSV"))
out_dir<-Sys.getenv("SEPSIS_ML_OUTPUT"); png_dir<-file.path(out_dir,"figures"); dir.create(png_dir,recursive=TRUE,showWarnings=FALSE)
lasso_rule<-"lambda.1se"; nfolds<-10L                      # 严格规则；可手动改成lambda.min，代码不会自动放宽
boruta_p<-.01; boruta_runs<-200L; threads<-4L             # 沿用模板：mcAdj=TRUE、getImpRfZ、holdHistory=TRUE
winsor_probs<-c(.01,.99)                                  # 仅生命体征/检验缩尾；改成c(0,1)关闭缩尾
fig_width<-12; fig_height<-9; fig_dpi<-600; font_family<-"Arial"  # 原模板TNR实际注册为Arial
stopifnot(lasso_rule%in%c("lambda.1se","lambda.min"),train_ratio>0,train_ratio<1,missing_limit>=0,missing_limit<1)

# ==================== 2. 导入数据、重命名、形式转换 ====================
# 明确候选列，防止把ID、住院时长、其他结局或随访字段自动放进预测模型。
rename_map<-c(age="Age",gender="Gender",race="Race",sofa="SOFA",oasis="OASIS",charlson="CCI",
  icd2_hypertension="Hypertension",icd2_diabetes="Diabetes",icd2_copd="COPD",icd2_stroke="Stroke",icd2_hf="HF",
  icd_liver="Liver_disease",icd2_ckd="CKD",icd_malignancy="Malignancy",icu_type="ICU_Type",
  lung_infection="Lung_infection",gi_infection="GI_infection",genito_urinary_infection="GU_infection",other_infection="Other_infection",
  vs_24h_map_min="MAP",vs_24h_heart_rate_max="Heart_Rate",vs_24h_resp_rate_max="Resp_Rate",vs_24h_temp_max="Temperature",
  lab_fst_24h_hemoglobin_first="Hb",lab_fst_24h_wbc_first="WBC",lab_fst_24h_platelet_first="Plt",
  lab_fst_24h_alt_first="ALT",lab_fst_24h_ast_first="AST",lab_fst_24h_tbil_first="Bilirubin",lab_fst_24h_alb_first="Albumin",
  lab_fst_24h_creatinine_first="Creatinine",lab_fst_24h_bun_first="BUN",lab_fst_24h_aptt_first="APTT",lab_fst_24h_lactate_first="Lactate",
  lab_fst_24h_na_first="Sodium",lab_fst_24h_k_first="Potassium",lab_fst_24h_ca_first="Calcium",lab_fst_24h_cl_first="Chloride",
  lab_fst_24h_rbc_first="RBC",lab_fst_24h_hct_first="Hct",lab_fst_24h_mcv_first="MCV",lab_fst_24h_rdw_first="RDW",
  itvtn_24h_vent_tag="MV",itvtn_24h_rrt_tag="RRT",drug_24h_vaso_tag="Vaso",drug_24h_insulin_tag="Insulin",
  gv_cv_icu_mixed="GV",shr_first24h_chem="SHR",hgi_first24h_chem="HGI")
binary_vars<-c("Hypertension","Diabetes","COPD","Stroke","HF","Liver_disease","CKD","Malignancy",
  "Lung_infection","GI_infection","GU_infection","Other_infection","MV","RRT","Vaso","Insulin")
factor_vars<-c("Gender","Race","ICU_Type",binary_vars); glycemic<-c("GV","SHR","HGI")
read_cohort<-function(file){
  if(!file.exists(file)) stop("文件不存在：",file)
  r<-as.data.frame(fread(file,encoding="UTF-8",na.strings=c("","NA","NULL","NaN"))); need<-c("stay_id","subject_id","hadm_id","day28_outcome","hba1c","glucose_first24h_chem")
  if(length(setdiff(need,names(r)))) stop("缺少必要列：",paste(setdiff(need,names(r)),collapse=", "))
  used<-intersect(names(rename_map),names(r)); d<-r[c(need,used)]; names(d)<-c("stay_id","subject_id","hadm_id","Outcome","HbA1c","Glucose",unname(rename_map[used]))
  for(id in c("stay_id","subject_id","hadm_id")) d[[id]]<-as.character(d[[id]])
  for(v in setdiff(names(d),c("stay_id","subject_id","hadm_id"))){
    z<-trimws(as.character(d[[v]])); z[toupper(z)%in%c("","NA","NULL","NAN")]<-NA_character_
    if(v%in%binary_vars){bad<-!is.na(z)&!toupper(z)%in%c("0","1","NO","YES","FALSE","TRUE","N","Y"); if(any(bad)) stop(v,"出现非0/1分类值")
      d[[v]]<-ifelse(is.na(z),NA_character_,ifelse(toupper(z)%in%c("1","YES","TRUE","Y"),"Yes","No"))
    } else if(v%in%factor_vars){d[[v]]<-z
    } else {a<-suppressWarnings(as.numeric(z)); if(any(!is.na(z)&is.na(a))) stop(v,"出现不能转为数值的内容"); a[!is.finite(a)]<-NA_real_; d[[v]]<-a}
  }
  if(any(vapply(d[c("stay_id","subject_id","hadm_id")],function(x)anyNA(x)||any(!nzchar(x))||anyDuplicated(x)>0,logical(1)))) stop("患者/住院/ICU ID缺失或重复")
  if(anyNA(d$Outcome)||!all(d$Outcome%in%0:1)||length(unique(d$Outcome))!=2) stop("Outcome需要为无缺失的0/1，且两类均存在")
  if(anyNA(d[c("HbA1c","Glucose")])) stop("当前队列的HbA1c或初次生化血糖存在缺失，请核对输入")
  d
}
cohort<-lapply(paths,read_cohort); mimic<-cohort$MIMICIV; nwicu<-cohort$NWICU
cat("MIMIC-IV：",nrow(mimic),"人；NWICU：",nrow(nwicu),"人。\n")

# ==================== 3. 两库共同变量：缺失率均≤20% ====================
all_vars<-unname(rename_map); common_vars<-Reduce(intersect,list(all_vars,names(mimic),names(nwicu)))
missing_pct<-function(d,v) if(v%in%names(d)) 100*mean(is.na(d[[v]])) else NA_real_
missing_table<-data.frame(Variable=all_vars,Source_column=names(rename_map),Common=all_vars%in%common_vars,
  MIMIC_missing_pct=vapply(all_vars,function(v) missing_pct(mimic,v),numeric(1)),NWICU_missing_pct=vapply(all_vars,function(v) missing_pct(nwicu,v),numeric(1)))
missing_table$Keep<-with(missing_table,Common & !is.na(MIMIC_missing_pct) & !is.na(NWICU_missing_pct) & MIMIC_missing_pct<=100*missing_limit & NWICU_missing_pct<=100*missing_limit)
vital_exceptions<-c("MAP","Heart_Rate","Resp_Rate","Temperature")
missing_table$Vital_exception<-missing_table$Variable%in%vital_exceptions & missing_table$Common
missing_table$Keep<-missing_table$Keep | missing_table$Vital_exception
missing_table$Reason<-ifelse(missing_table$Keep,ifelse(missing_table$Vital_exception & (missing_table$MIMIC_missing_pct>20|missing_table$NWICU_missing_pct>20),"Prespecified vital-sign exception; training-only imputation","Common; both missingness <=20%"),"Not common or missingness >20%")
candidates<-missing_table$Variable[missing_table$Keep]; factor_vars<-intersect(factor_vars,candidates)
winsor_vars<-intersect(unname(rename_map[grepl("^vs_|^lab_",names(rename_map))]),candidates)
if(length(candidates)<2) stop("符合共同缺失率条件的候选变量不足2项")
print(missing_table,row.names=FALSE); cat("保留",length(candidates),"项；剔除：",paste(setdiff(all_vars,candidates),collapse=", "),"\n")

# ==================== 4. MIMIC-IV分层拆分70:30 ====================
# NWICU全体保留为外部验证；其结局不参与插补、HGI公式估计或特征筛选。
set.seed(seed); train_index<-sort(unlist(lapply(split(seq_len(nrow(mimic)),mimic$Outcome),function(i) sample(i,floor(length(i)*train_ratio)))))
train<-mimic[train_index,,drop=FALSE]; internal<-mimic[-train_index,,drop=FALSE]; external<-nwicu
cohort_table<-data.frame(Cohort=c("MIMIC-IV training","MIMIC-IV internal","NWICU external"),N=c(nrow(train),nrow(internal),nrow(external)),Deaths=c(sum(train$Outcome),sum(internal$Outcome),sum(external$Outcome)))
print(cohort_table,row.names=FALSE)

# ==================== 5. 训练集插补、缩尾、HGI和编码 ====================
# 沿用模板的预测预处理：连续变量中位数、分类变量众数；参数只从训练集估计。
# HGI=HbA1c−预测HbA1c；HbA1c~初次生化血糖仅在训练集拟合，再固定应用到另两组。
# 输入CSV原HGI用于可用性筛查，本节重算建模用HGI；HbA1c/血糖本身不进入候选模型。
mode_value<-function(x){tt<-table(x); if(!length(tt)) stop("训练列全部缺失"); names(tt)[which.max(tt)]}
fit_prep<-function(d){
  h<-coef(lm(HbA1c~Glucose,data=d)); if(length(h)!=2||any(!is.finite(h))) stop("HGI公式拟合失败")
  d$HGI<-d$HbA1c-(h[1]+h[2]*d$Glucose); nums<-setdiff(candidates,factor_vars)
  med<-vapply(d[nums],median,numeric(1),na.rm=TRUE); if(any(!is.finite(med))) stop("训练集存在全缺失数值列")
  lim<-lapply(d[winsor_vars],quantile,probs=winsor_probs,na.rm=TRUE,names=FALSE)
  modes<-vapply(d[factor_vars],mode_value,character(1)); lev<-lapply(d[factor_vars],function(x) sort(unique(na.omit(x))))
  for(v in names(lev)){ref<-switch(v,Gender="F",Race="White",ICU_Type="Medical ICU","No"); if(ref%in%lev[[v]]) lev[[v]]<-c(ref,setdiff(lev[[v]],ref))}
  list(hgi=h,median=med,bounds=lim,mode=modes,levels=lev,vars=candidates)
}
apply_prep<-function(d,p){
  d$HGI<-d$HbA1c-(p$hgi[1]+p$hgi[2]*d$Glucose)
  for(v in names(p$median)){d[[v]][is.na(d[[v]])]<-p$median[[v]]; if(v%in%names(p$bounds)) d[[v]]<-pmin(pmax(d[[v]],p$bounds[[v]][1]),p$bounds[[v]][2])}
  for(v in names(p$levels)){z<-as.character(d[[v]]); z[is.na(z)|!z%in%p$levels[[v]]]<-p$mode[[v]]; d[[v]]<-factor(z,levels=p$levels[[v]])}
  d<-d[p$vars]; if(anyNA(d)) stop("预处理后仍有缺失"); d
}
prep<-fit_prep(train); z_train<-apply_prep(train,prep)
constant_vars<-names(z_train)[vapply(z_train,function(x) length(unique(x))<2,logical(1))]; prep$vars<-setdiff(candidates,constant_vars)
if(length(prep$vars)<2) stop("去除训练集常量后候选变量不足2项")
z_train<-apply_prep(train,prep); z_internal<-apply_prep(internal,prep); z_external<-apply_prep(external,prep)
encode<-function(d){tt<-terms(~.,data=d); mm<-model.matrix(tt,d); a<-attr(mm,"assign"); keep<-a!=0; parent<-attr(tt,"term.labels")[a[keep]]
  mm<-mm[,keep,drop=FALSE]; storage.mode(mm)<-"double"; attr(mm,"parent")<-setNames(parent,colnames(mm)); mm}
x_train<-encode(z_train); x_internal<-encode(z_internal); x_external<-encode(z_external); parent_map<-attr(x_train,"parent")
stopifnot(identical(colnames(x_train),colnames(x_internal)),identical(colnames(x_train),colnames(x_external)),all(is.finite(x_train)),all(is.finite(x_internal)),all(is.finite(x_external)))
y_train<-factor(train$Outcome,levels=c(0,1)); cat("完成插补；训练常量剔除：",paste(constant_vars,collapse=", "),"\n")
cat(sprintf("训练集HGI = HbA1c - (%.6f + %.6f × Glucose)\n",prep$hgi[1],prep$hgi[2]))

# ==================== 6. LASSO：10折CV与变量提取 ====================
# 保留模板family=binomial、alpha=1、standardize=TRUE、10折deviance；所有选择都在训练集内。
# 修正旧模板注释1se/实际min不一致：下列提取和表格统一使用第1节lasso_rule；不自动改规则。
if(min(table(y_train))<nfolds) stop("训练集较少结局类别不足以进行设定的交叉验证折数")
set.seed(seed); foldid<-integer(nrow(train)); for(g in levels(y_train)){ii<-sample(which(y_train==g)); foldid[ii]<-rep(seq_len(nfolds),length.out=length(ii))}
set.seed(seed); cvfit<-cv.glmnet(x_train,y_train,family="binomial",alpha=1,standardize=TRUE,nfolds=nfolds,foldid=foldid,type.measure="deviance")
b_min<-as.matrix(coef(cvfit,s="lambda.min"))[,1]; b_1se<-as.matrix(coef(cvfit,s="lambda.1se"))[,1]; b<-if(lasso_rule=="lambda.1se") b_1se else b_min
selected_dummy<-setdiff(names(b)[b!=0],"(Intercept)"); select_lasso_feature<-unique(unname(parent_map[selected_dummy]))
lasso_table<-data.frame(Feature=colnames(x_train),Variable=unname(parent_map),Coef_min=unname(b_min[colnames(x_train)]),Coef_1se=unname(b_1se[colnames(x_train)]),Selected=colnames(x_train)%in%selected_dummy)
cc<-data.frame(lambda=log(cvfit$lambda),cvm=cvfit$cvm,cvup=cvfit$cvup,cvlo=cvfit$cvlo,nzero=cvfit$nzero)
lambda_table<-do.call(rbind,lapply(c("lambda.min","lambda.1se"),function(rule){i<-which.min(abs(cvfit$lambda-cvfit[[rule]])); data.frame(Rule=rule,Lambda=cvfit[[rule]],CV_deviance=cc$cvm[i],SE=cvfit$cvsd[i],Nonzero_coefficients=cc$nzero[i],Used=rule==lasso_rule)}))
cat("LASSO（",lasso_rule,"）选中",length(select_lasso_feature),"个原始变量：\n"); print(select_lasso_feature)

# ==================== 8. Boruta：参数、筛选和重要性 ====================
# 保留pValue=.01、mcAdj=TRUE、maxRuns=200、doTrace=2、holdHistory=TRUE、getImpRfZ、15线程。
# getImpRfZ默认500棵树；将模板嵌套ranger.params改为直接传参，使partition选项真正生效。
set.seed(seed); boruta_selection<-Boruta(x=z_train,y=y_train,pValue=boruta_p,mcAdj=TRUE,maxRuns=boruta_runs,doTrace=2,holdHistory=TRUE,
  getImp=getImpRfZ,num.trees=500L,num.threads=threads,respect.unordered.factors="partition")
decisions<-boruta_selection$finalDecision; select_boruta_feature<-names(decisions)[decisions=="Confirmed"]
boruta_table<-data.frame(Variable=names(decisions),attStats(boruta_selection),row.names=NULL); boruta_table<-boruta_table[order(-boruta_table$medianImp),]
cat("Boruta Confirmed：",length(select_boruta_feature),"项；Tentative：",sum(decisions=="Tentative"),"项。\n"); print(boruta_table,row.names=FALSE)

# ==================== 10. 交集、汇总表、下一步机器学习数据 ====================
# 按原始临床变量取交集：分类变量有任一虚拟列入选LASSO，即该变量入选；Boruta只取Confirmed。
# 不强制加入GV/SHR/HGI，不把Tentative当Confirmed；真实未入选状态保留在汇总表。
selected_common<-prep$vars[prep$vars%in%intersect(select_lasso_feature,select_boruta_feature)]
selection_table<-data.frame(Variable=candidates,LASSO=candidates%in%select_lasso_feature,Boruta=as.character(decisions[candidates]),Intersection=candidates%in%selected_common)
selection_table$Boruta[selection_table$Variable%in%constant_vars]<-"Constant (not tested)"
selected_columns<-names(parent_map)[parent_map%in%selected_common]
model_data<-function(raw,z) data.frame(stay_id=raw$stay_id,Outcome=factor(raw$Outcome,levels=c(0,1)),z[,selected_common,drop=FALSE],check.names=FALSE)
train_dt<-model_data(train,z_train); internal_dt<-model_data(internal,z_internal); nwicu_dt<-model_data(external,z_external)
unseen_table<-do.call(rbind,lapply(c("Internal","NWICU"),function(s){d<-if(s=="Internal") internal else external
  data.frame(Cohort=s,Variable=names(prep$levels),Unseen_levels_n=vapply(names(prep$levels),function(v) sum(!is.na(d[[v]])&!d[[v]]%in%prep$levels[[v]]),integer(1)))}))
settings<-list(seed=seed,train_ratio=train_ratio,missing_limit=missing_limit,lasso_rule=lasso_rule,nfolds=nfolds,boruta_p=boruta_p,boruta_runs=boruta_runs,
  boruta_mcAdj=TRUE,boruta_trees=500L,boruta_threads=threads,boruta_factors="partition",winsor_probs=winsor_probs,winsor_vars=winsor_vars)
# 下一步调参CV需在每个训练折重做预处理和筛选；下面保留原始拆分数据供其调用。
# 本节LASSO的CV曲线是选择lambda的内部过程，不作为独立预测效能结果。
out<-list(settings=settings,paths=paths,source_md5=tools::md5sum(paths),rename_map=rename_map,candidates=candidates,selected=selected_common,selected_columns=selected_columns,
  raw=list(train=train,internal=internal,nwicu=external),data=list(train=train_dt,internal=internal_dt,nwicu=nwicu_dt),prep=prep,
  x_train=x_train[,selected_columns,drop=FALSE],x_internal=x_internal[,selected_columns,drop=FALSE],x_nwicu=x_external[,selected_columns,drop=FALSE],
  parent_map=parent_map,lasso=cvfit,boruta=boruta_selection,selection_table=selection_table,missing_table=missing_table,cohort_table=cohort_table,foldid=foldid,
  constant_vars=constant_vars,unseen_table=unseen_table,session=sessionInfo())
saveRDS(out,file.path(out_dir,"private_checkpoints/特征工程结果.rds"))
sheets<-list("队列人数"=cohort_table,"缺失率筛查"=missing_table,"筛选交集"=selection_table,"LASSO系数"=lasso_table,"LASSO参数"=lambda_table,"Boruta重要性"=boruta_table,"新分类水平"=unseen_table)
for(n in names(sheets)) write.csv(sheets[[n]],file.path(out_dir,"results",paste0("特征工程_",n,".csv")),row.names=FALSE,fileEncoding="UTF-8")
saveRDS(list(lasso_table=lasso_table,lambda_table=lambda_table,boruta_table=boruta_table,selection=selection_table),file.path(out_dir,"private_checkpoints/特征工程图表对象.rds"))
for(n in names(out$raw)) write.csv(out$raw[[n]],file.path(out_dir,"private_model_inputs",paste0(n,"_原NA拆分.csv")),row.names=FALSE,fileEncoding="UTF-8")
stopifnot(nrow(mimic)==5195L,sum(mimic$Outcome)==751L,nrow(nwicu)==1140L,sum(nwicu$Outcome)==260L,!any(train$subject_id%in%internal$subject_id))
cat("\n完成：LASSO",length(select_lasso_feature),"项；Boruta Confirmed",length(select_boruta_feature),"项；交集",length(selected_common),"项。\n")
print(selection_table,row.names=FALSE); cat("\n三个指标：\n"); print(selection_table[selection_table$Variable%in%glycemic,],row.names=FALSE)
cat("\n结果目录：",out_dir,"\n交集变量：",paste(selected_common,collapse=", "),"\n")
if(!length(selected_common)) warning("交集为空；结果已保存，下一步不能直接建模，请先审核筛选结果。")
