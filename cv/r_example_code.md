
# ABC

pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance, htmltools,kableExtra)
conflict_prefer_all("dplyr", quiet = TRUE)
options(scipen = 999)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))




samp_priors <- function(n=1) {
  prior_samples <- tibble(
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
    weight_exam = runif(n,0,1)
  )
  return(prior_samples)
}

reject_abc <- function(simulation_function, data, num_iterations = 5000, n_try=500, return_dat="test_data") {
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  data <- data |> as.data.table()
  target_data <- case_when(
    return_dat == "test_data" ~ list(target_data_test <- data[data$expMode2 == "Test", ]), 
    return_dat == "train_data" ~ list( target_data_train <- data[data$expMode2 == "Train", ]),
    return_dat == "train_data, test_data" ~ list( target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ])
    ) |> pluck(1)

  tol <- target_data[, .(m = mean(y), sd = sd(y)), by = .(x,Block)][, .(tol = mean(sd) * tolM)] %>% as.numeric()

  start_tol = round(tol,3); 
  abc <- list(); abc$dist_sd <- vector("list", num_iterations)
  try_count=0;
  inc_count=0; # number of tolerance increases
  cur_tol_success=0;
  success_rate=0;
  closest_mean_errors <- numeric(0)

  t1=system.time({
  for(j in 1:num_iterations) {

    found=0;
    iter_count=0; # number of unsuccessful samples taken this iteration

    while(found==0) {
    try_count=try_count+1;
    iter_count = iter_count+1;

    current_theta <- samp_priors()
    sim_data <- simulation_function(data, current_theta$c, current_theta$lr, current_theta$weight_exam,
     input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) 

    # dist_sd <- target_data |> mutate(pred=sim_data,error=abs(y-pred)) |> group_by(id,condit,x) |> 
    #   summarise(mean_error=mean(error),.groups="keep") |> 
    #   group_by(id,condit) |> 
    #   summarise(mean_error=mean(mean_error),.groups="keep") 

    target_data[, pred := sim_data]
    target_data[, error := abs(y - pred)]
    intermediate_result <- target_data[, .(mean_error = mean(error)), by = .(id, condit, x,Block)]

    if(return_dat == "train_data, test_data") {
          # Calculate mean errors separately for training and testing, then take their average
          train_error <- intermediate_result[Block < 4, .(mean_error = mean(mean_error)),by = .(id, condit)]
          test_error <- intermediate_result[Block >= 4, .(mean_error = mean(mean_error)),by = .(id, condit)]
          # Ensure equal weighting by taking the simple average of train and test mean errors
          dist_sd <- merge(train_error, test_error, 
            by = c("id", "condit"), all = TRUE, suffixes = c("_train", "_test")) 
          dist_sd[, mean_error := (mean_error_train + mean_error_test) / 2][, c("mean_error_train", "mean_error_test") := NULL]
        } else {
          # Original computation for cases other than combined train and test data
          dist_sd <- intermediate_result[, .(mean_error = mean(mean_error)), by = .(id, condit)]
        }

    closest_mean_errors <- sort(c(closest_mean_errors, dist_sd$mean_error))
    if (length(closest_mean_errors) > 5) {
      closest_mean_errors <- head(closest_mean_errors, 5)  
    }

    if(dist_sd$mean_error< tol) {
      abc$dist_sd[[j]] <- cbind(current_theta,dist_sd,tol,inc_count,iter_count)
      found=1
      #try_count=try_count+1;
      cur_tol_success = cur_tol_success+1;
      success_rate = cur_tol_success/try_count;

      tol <- (tol) + ((mean(closest_mean_errors)-tol)/2) 

      } else if (try_count > n_try && success_rate < min_accept_rate){
       
        average_closest_error <- mean(closest_mean_errors)
        bump_tol <- abs(rnorm(n=1,mean=(3*iter_count/(n_try)),sd=1))
       # print(bump_tol)
        #tol <- (tol*tolInc) + (abs(average_closest_error-tol)/2) #* tolInc  # Adjust tolerance
        tol <- (tol) + ((average_closest_error-tol)/2) + bump_tol   # (max(0,min_accept_rate-success_rate)*tol)   #(min_accept_rate/(success_rate+.0001))
        #tol=tol*tolInc
        inc_count=inc_count+1;
        try_count=0;
        cur_tol_success=0;

      #  if (data$id[1] %in% watch_ids){
      #   message(paste0("increase tol(",round(tol,3),") for subject", data$id[1]," current iteration: ",j,"\n",
      #   "cur accept rate: ",round(success_rate,4),"\n",
      #    "avg closest error: ",round(average_closest_error,2), " new tol: ",round(tol,2),"\n"))
      #   }

      }
    }
  }  
  })

  abc <- rbindlist(abc$dist_sd) |> arrange(mean_error) |> relocate(id,condit) |> mutate(t=round(t1[3],2))
  best <- abc |> head(1) |> round_tibble(7)
  message((paste0("\n", data$id[1], " completed in: ", round(t1[3],2))))
  message((paste0("inc count: ",inc_count," Start tol: ",start_tol ," end tol: ", round(tol,2)," success rate: ",round(success_rate,4))))
  message((paste0("Best c: ", best$c, " Best lr: ", round(best$lr,3)," Best_w:",round(best$weight_exam,2)," Best error: ", round(best$mean_error,2),"\n\n\n")))
  return(abc)
}

####################################
####################################

run_abc_tests <- function(simulation_function, data_list, return_dat, ids) {
  p_abc()
  cat("\nRunning ABC Test: ",as.character(substitute(simulation_function))," ",return_dat, "\n", "Parallel Execution\n")

  t1 <- system.time({
    results <- future_map(data_list, 
                          ~reject_abc(simulation_function = simulation_function, 
                                      #prior_samples = prior_samples, 
                                      data = .x, 
                                      num_iterations = num_iterations, 
                                      n_try = n_try,
                                      return_dat = return_dat), 
                          .options = furrr_options(seed = TRUE, scheduling = Inf),
                          future.stdout = NA) %>% 
              setNames(ids)
  })
  
  print(t1[3])
  return(results)
}

####################################
####################################

args <- commandArgs(trailingOnly = TRUE)
num_iterations = ifelse(length(args) > 0, as.numeric(args[1]), 70)
n_try = ifelse(length(args) > 1, as.numeric(args[2]), 500)
tolM <<- ifelse(length(args) > 2, as.numeric(args[3]), .85)
tolInc <<- ifelse(length(args) > 3, as.numeric(args[4]), 1.01)
min_accept_rate <<- ifelse(length(args) > 4, as.numeric(args[5]), .1)

cMean <<- -6; 
cSig <<- 4.0; 
lrSig <<- 4.0
prior_samples <- samp_priors(n=1000) 

p_abc <- function(){
  message(paste0("\n","cMean: ",cMean," cSig: ",cSig," lrSig: ",lrSig,"\n", "tolM: ",
  tolM, " tolInc: ",tolInc," accept_rate: ",min_accept_rate, "\n",
  "num_iterations: ",num_iterations," n_try: ",n_try,"\n",
  "median-c=",round(median(prior_samples$c),4)," median-lr=",round(median(prior_samples$lr),4),"\n",
  " median-w=",round(median(prior_samples$weight_exam),3),"\n"))
  }

parallel <<- 1
#parallel <<- runif(1) <.5

if (parallel==1) {
  print(nc <- future::availableCores())
  future::plan(multicore, workers = nc-2)
  message("Running in parallel\n")
} else if(parallel==2){
spec <- make_spec()
#pc <- parallel::makeCluster(type='PSOCK', master=macpro, spec)
pc <- parallel::makeCluster(type='PSOCK', master=tg_m1, spec)

plan(cluster, workers = pc)
} else {
  message("Running in serial\n")
}


# parallel::stopCluster(pc)

run_function <- ifelse(parallel>0,run_abc_tests, run_abc_tests_serial)
 
####################################
# ids1 <- 1
# ids1 <- c(1,33,66)
# ids1 <- as.numeric(levels(ds$id))[1:8]


run_e1 <- function(){
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
#ds <- readRDS(here::here("data/e3_md_02-23-24.rds"))  |> as.data.table()
nbins <- 3
ds <- ds |> group_by(id) |> 
  mutate(Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) 
ids1 <- as.numeric(levels(ds$id))
subjects_data <-  ds |> filter(id %in% ids1)  %>% with(split(.,f =c(id), drop=TRUE))
save_folder <- paste0("e1_hybrid_n_iter_",num_iterations,"_ntry_",n_try,"_",format(Sys.time(),"%M%OS"))
dir.create(paste0("data/abc_reject/",save_folder))

#save_folder <- 'n_iter_200_ntry_300_5623'

t1<- system.time( hybrid_test <- run_function(full_sim_hybrid, subjects_data, "test_data", ids1) )
save_abc_test_results(hybrid_test, "hybrid", "Test", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(hybrid_test); gc()

t1<- system.time( hybrid_test_train <- run_function(full_sim_hybrid, subjects_data, "train_data, test_data", ids1) )
save_abc_test_results(hybrid_test_train, "hybrid", "Test_Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(hybrid_test_train); gc()

t1<- system.time( hybrid_train <- run_function(full_sim_hybrid, subjects_data, "train_data", ids1) )
save_abc_test_results(hybrid_train, "hybrid", "Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(hybrid_train)
rm(ds,ids1,subjects_data); gc()
print("finish e1")

}




### Modelling Results

```{r}
#| cache: false

#ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
ds <- readRDS(here::here("data/e2_md_02-23-24.rds"))  |> as.data.table()
nbins <- 3

fd <- readRDS(here("data/e2_08-21-23.rds"))
test <- fd |> filter(expMode2 == "Test") 
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

trainAvg <- fd |> filter(expMode2 == "Train") |> group_by(id) |> 
  mutate(tr=trial,x=vb,Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) |> 
  group_by(id,condit,vb,x,Block) |> 
  summarise(dist=mean(dist),y=mean(vx))

input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)
ids2 <- c(1,66,36)

#file_name <- "e2_n_iter_50_ntry_200_2506"
#file_name <- "n_iter_400_ntry_100_2944"
#file_name <- "e2_n_iter_100_ntry_200_3436"
file_name <- "e2_n_iter_200_ntry_300_2344"

ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 
```

```{r}
#| cache: false

generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}

#future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))


#  post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#     generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 50, return_dat="test_data")
#     })) |> 
#    select(Fit_Method,pp,-data) |>  
#    unnest(pp) |>  filter(expMode2=="Test") |> as.data.table()

# saveRDS(post_dat, here("data/model_cache/post_dat_e2.rds"))

post_dat <- readRDS(here("data/model_cache/post_dat_e2.rds"))



post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  mutate(error2 = y - pred) |>
  summarise(y = mean(y), pred = mean(pred), error = y - pred, error2=mean(error2)) |> as.data.table()

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
#rm(post_dat_avg)

post_dat_l <- post_dat_l |> mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))


# organize training data predictions
 pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
   generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 20, return_dat="train_data")
    })) |>
   select(Fit_Method,pp,-data) |>
  unnest(pp) |> as.data.table() |> filter(expMode2=="Train")

#saveRDS(pd_train, here("data/model_cache/pd_train.rds"))

#pd_train <- readRDS(here("data/model_cache/pd_train.rds"))

nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l  |>
  mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

#plan(sequential)

```

### Group level aggregations

```{r}
#| eval: true
#| label: tbl-htw-modelError
#| tbl-cap: "Mean model errors predicting testing data, aggregated over all participants and velocity bands. Note that Fit Method refers to how model parameters were optimized, while error values reflect mean absolute error for the 6 testing bands"
post_tabs <- abc_tables(post_dat,post_dat_l)

post_tabs$agg_pred_full |> 
  mutate(Fit_Method=rename_fm(Fit_Method)) |>
  flextable::tabulator(rows=c("Fit_Method","Model"), columns=c("condit"),
                       `ME` = as_paragraph(mean_error)) |> as_flextable()

#post_tabs$agg_pred_full |> pander::pandoc.table()

```



```{r fig.height=12,fig.width=11}
#| label: fig-htw-resid-pred
#| column: page-inset-right
#| fig-cap: A) Model residuals for each combination of training condition, fit method, and model. Residuals reflect the difference between observed and predicted values. Lower values indicate better model fit. Note that y axes are scaled differently between facets. B) Full posterior predictive distributions vs. observed data from participants.Points represent median values, thicker intervals represent 66% credible intervals and thin intervals represent 95% credible intervals around the median. 
#| fig-height: 19
#| fig-width: 12
##| layout: [[45,-5, 45], [100]]
##| fig-subcap: ["Model Residuals - training data", "Model Residuals - testing data","Full posterior predictive distributions vs. observed data from participants."]
train_resid <- pd_train |> group_by(id,condit,Model,Fit_Method, Block) |> 
  summarise(y = mean(y), pred = mean(pred), error = y - pred) |>
  ggplot(aes(x = Block, y = abs(error), fill=Model)) + 
  stat_bar + 
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, scales="free",ncol=2) +
  scale_fill_manual(values=wes_palette("AsteroidCity2"))+
  labs(title="Model Residual Errors - Training Stage", y="RMSE", x= "Training Block") +
  theme(legend.title = element_blank(), legend.position="top")

test_resid <-  post_dat |> 
   group_by(id,condit,x,Model,Fit_Method,rank) |>
   summarize(error=mean(abs(y-pred)),n=n()) |>
   group_by(id,condit,x,Model,Fit_Method) |>
   summarize(error=mean(error)) |>
  mutate(vbLab = factor(paste0(x,"-",x+200))) |>
  ggplot(aes(x = vbLab, y = abs(error), fill=Model)) + 
  stat_bar + 
  scale_fill_manual(values=wes_palette("AsteroidCity2"))+
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, axes = "all",ncol=2,scale="free") +
  labs(title="Model Residual Errors - Testing Stage",y="RMSE", x="Velocity Band") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 

group_pred <- post_dat_l |> 
  mutate(vbLab = factor(paste0(x,"-",x+200),levels=levels(testAvg$vb))) |>
  ggplot(aes(x=val,y=vbLab,col=Resp)) + 
  stat_pointinterval(position=position_dodge(.5), alpha=.9) + 
  scale_color_manual(values=wes_palette("AsteroidCity2"))+
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, axes = "all",ncol=2,scale="free") +
  labs(title="Posterior Predictions - Testing Stage",y="Velocity Band (lower bound)", x="X Velocity") +
theme(legend.title=element_blank(),axis.text.y = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


((train_resid | test_resid) / group_pred) +
  plot_layout(heights=c(1,1.5)) & 
  plot_annotation(tag_levels = list(c('A1','A2','B')),tag_suffix = ') ') & 
  theme(plot.tag.position = c(0, 1))

```


## Deviation Predictions

```{r}
#| fig-height: 12
#| fig-width: 11
 post_dat_l |> 
  mutate(vbLab = factor(paste0(x,"-",x+200),levels=levels(testAvg$vb))) |>
  ggplot(aes(x=condit,y=dist,fill=vbLab)) + 
  stat_bar + 
  #facet_wrap(~Resp)
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~Resp, axes = "all",ncol=3,scale="free")

```



pacman::p_load(tidyverse,data.table,abc,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 


tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()


avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> rbind(dsc |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()



input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer



generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 7),
    lr = runif(n, 0.000001, 7),
  )
  return(prior_samples)
}


full_sim_exam <- function(data, c, lr,pred_fun=exam.response, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  if (train_data$condit[1] != "Varied") {
    trainVec <- c(0, trainVec)
  }
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  
  test_prediction <- map_dbl(test_data$x, ~ pred_fun(.x, c, input_layer, 
                                                     output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}


full_sim_alm <- function(data, c, lr,pred_fun=alm.responseOnly, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

full_sim_alt_exam <- function(data, c, lr, pred_fun=alt_exam, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  trainVecY=train_data$y
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alt_exam(.x, c, input_layer, output_layer, train_results$wm,  trainVecX=trainVec, trainVecY=trainVecY))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}



# General function for ABC fits with model flexibility
run_abc_fits <- function(data, input_layer, output_layer, simulation_function, n_prior_samples = 1000, tol = 0.01, return_dat = "train_data,test_data") {
  
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  
  # Generate prior samples
  prior_samples <- generate_prior_c_lr(n_prior_samples)
  
  # Convert the string to a function
  if (is.character(simulation_function)) {
    simulation_function <- get(simulation_function, mode = "function")
  } else {
    simulation_function <- match.fun(simulation_function)
  }
  
  plan(multisession)
  # Run simulations for each prior sample
  simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE))
  
  # Extract target data
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  # ABC for train and test data
  abc_train_test <- abc(
    target = target_data_train_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results),
    tol = tol,
    method = "rejection"
  )
  
  # ABC for test data only
  abc_test <- abc(
    target = target_data_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[85:90, ]),
    tol = tol,
    method = "rejection"
  )
  
  # ABC for train data only
  abc_train <- abc(
    target = target_data_train,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[1:84, ]),
    tol = tol,
    method = "rejection"
  )
  
  # Return results
  tibble::lst(abc_train_test = abc_train_test, abc_test = abc_test, abc_train,data=data,n_prior_samples,prior_samples,tol, simulation_function, targets=tibble::lst(target_data_train_test, target_data_test))
}



n_priors=500000
args_list <- tibble::lst(
  abc_v_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_v_alt_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors),
  abc_v_alm=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_alm=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_c_alt_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors)
)

plan(multisession)
abc_list <- future_map(args_list, ~do.call(run_abc_fits, .x))


saveRDS(abc_list,here::here(paste0("data/model_cache/abc_group500k_",format(Sys.time(), "%H_%M_%OS"),".rds")))


@fig-design-e2 illustrates the design of Experiment 2. The stages of the experiment (i.e. training, testing no-feedback, test with feedback), are identical to that of Experiment 1. The only change is that Experiment 2 participants train, and then test, on bands in the reverse order of Experiment 1 (i.e. training on the softer bands; and testing on the harder bands). 
```{dot}
//| label: fig-design-e2
//| fig-cap: "Experiment 2 Design. Constant and Varied participants complete different training conditions. The training and testing bands are the reverse of Experiment 1. "
//| fig-width: 8.0
//| fig-height: 2.5
//| fig-responsive: false
//| column: screen-inset-right
digraph {
  graph [layout = dot, rankdir = LR]

  // define the global styles of the nodes
  node [shape = rectangle, style = filled]

  data1 [label = " Varied Training \n100-300\n350-550\n600-800", fillcolor = "#FF0000"]
  data2 [label = " Constant Training \n600-800", fillcolor = "#00A08A"]
  Test3 [label = "    Final Test \n  Novel With Feedback  \n800-1000\n1000-1200\n1200-1400", fillcolor = "#ECCBAE"]

  // edge definitions with the node IDs
  data1 -> Test1
  data2 -> Test1
  subgraph cluster {
    label = "Test Phase \n(Counterbalanced Order)"
    Test1 [label = "Test  \nNovel Bands  \n800-1000\n1000-1200\n1200-1400", fillcolor = "#ECCBAE"]
    Test2 [label = "  Test \n  Varied Training Bands  \n100-300\n350-550\n600-800", fillcolor = "#ECCBAE"]
    Test1 -> Test2
  }

  Test2 -> Test3
}

```







### Results

#### Testing Phase - No feedback. 

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw. 


##### Deviation From Target Band

Descriptive summaries testing deviation data are provided in @tbl-e2-test-nf-deviation and @fig-e2-test-dev. 
To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id). 

\begin{equation}
dist_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot band_{ij} + \beta_3 \cdot condit_{ij} \cdot band_{ij} + b_{0i} + b_{1i} \cdot band_{ij} + \epsilon_{ij}
\end{equation}

```{r}
#| label: tbl-e2-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]

result <- test_summary_table(testE2, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
result$constant |> kable()
result$varied |> kable()
# make kable table with smaller font size
# result$constant |> kbl(caption="Constant Testing - Deviation",booktabs=T,escape=F) |> kable_styling(font_size = 7)

```

```{r}
#| label: fig-e2-test-dev
#| fig-cap: E2. Deviations from target band during testing without feedback stage. 
testE2 |>  ggplot(aes(x = vb, y = dist,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
  labs(x="Band", y="Deviation From Target")
```




```{r}
#| label: tbl-e2-bmm-dist
#| tbl-cap: "Experiment 2. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band"
#contrasts(test$condit) 
# contrasts(testE2$vb)

modelName <- "e2_testDistBand_RF_5K"
e2_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE2,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp2 <- GetModelStats(e2_distBMM) |> kable(escape=F,booktabs=T)
mp2

e2_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),
          epred = TRUE, re_formula = NA) |> 
  pairs() |> gather_emmeans_draws()  |> 
   summarize(median_qi(.value),pd=sum(.value>0)/n()) |>
   select(contrast,Band=bandInt,value=y,lower=ymin,upper=ymax,pd) |> 
   mutate(across(where(is.numeric), \(x) round(x, 2)),
          pd=ifelse(value<0,1-pd,pd)) |>
   kable(caption="Contrasts")

coef_details <- get_coef_details(e2_distBMM, "conditVaried")


```


The model predicting absolute deviation showed a modest tendency for the varied training group to have lower deviation compared to the constant training group (β = `r coef_details$estimate`, 95% CI \[`r coef_details$conf.low`, `r coef_details$conf.high`\]),with 94% of the posterior distribution being less than 0. This suggests a potential benefit of training with variation, though the evidence is not definitive.

(SHOULD PROBABLY DO ALTERNATE ANALYSIS THAT ONLY CONSIDERS THE NOVEL EXTRAPOLATION BANDS)


```{r}
#| label: fig-e2-bmm-dist
#| fig-cap: E2. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals. 



e2_distBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
 condEffects()+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(testE2$vb), 
                     limits = c(0, 1400)) 
```


##### Discrimination between Velocity Bands

In addition to accuracy/deviation. We also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). @tbl-e2-test-nf-vx shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants. 

\begin{equation}
vx_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot bandInt_{ij} + \beta_3 \cdot condit_{ij} \cdot bandInt_{ij} + b_{0i} + b_{1i} \cdot bandInt_{ij} + \epsilon_{ij}
\end{equation}




```{r}
#| label: fig-e2-test-vx
#| fig-cap: E2 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band. 
#| fig-width: 11
#| fig-height: 9
testE2 %>% group_by(id,vb,condit) |> plot_distByCondit()

```


```{r}
#| label: tbl-e2-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 1

result <- test_summary_table(testE2, "vx","X Velocity" ,mfun = list(mean = mean, median = median, sd = sd))
result$constant |> kable()
result$varied |> kable()

```





```{r}
#https://nrennie.rbind.io/blog/2022-06-06-creating-flowcharts-with-ggplot2/
inNodes <- seq(1,6,1) %>% as.integer()
outNodes <- seq(300,1000,50)%>% as.integer()

stim <- "Stim"
resp <- "Response"
inFlow <- tibble(expand.grid(from=stim,to=inNodes)) %>% mutate_all(as.character)
outFlow <- tibble(expand.grid(from=outNodes,to=resp)) %>% mutate_all(as.character)

gd <- tibble(expand.grid(from=inNodes,to=outNodes)) %>% mutate_all(as.character) %>%
  rbind(inFlow,.) %>% rbind(.,outFlow)

g = graph_from_data_frame(gd,directed=TRUE)
coords2=layout_as_tree(g)
colnames(coords2)=c("y","x")
odf <- as_tibble(coords2) %>% 
  mutate(label=vertex_attr(g,"name"),
         type=c("stim",rep("Input",length(inNodes)),rep("Output",length(outNodes)),"Resp"),
         x=x*-1) %>%
  mutate(y=ifelse(type=="Resp",0,y),xmin=x-.05,xmax=x+.05,ymin=y-.35,ymax=y+.35)

plot_edges = gd %>% mutate(id=row_number()) %>%
  pivot_longer(cols=c("from","to"),names_to="s_e",values_to=("label")) %>%
                 mutate(label=as.character(label)) %>% 
  group_by(id) %>%
  mutate(weight=sqrt(rnorm(1,mean=0,sd=10)^2)/10) %>%
  left_join(odf,by="label") %>%
  mutate(xmin=xmin+.02,xmax=xmax-.02)

ggplot() + geom_rect(data = odf,
            mapping = aes(xmin = xmin, ymin = ymin, 
                          xmax = xmax, ymax = ymax, 
                          fill = type, colour = type),alpha = 0.5) +
  geom_text(data=odf,aes(x=x,y=y,label=label,size=3)) +
  geom_path(data=plot_edges,mapping=aes(x=x,y=y,group=id,alpha=weight)) +
  theme_void() + theme(legend.position = "none") 
  
   

```



```{r}


 
linear_function <- function(x) 2.2 * x + 30
exponential_function <- function(x) 200 * (1 - exp(-x/25))
quadratic_function <- function(x) 210 - (x - 50)^2 / 12

extrapLines <- list(geom_vline(xintercept=30,color="black",alpha=.4,linetype="dashed"),
                    geom_vline(xintercept=70,color="black",alpha=.4,linetype="dashed"))
 
linear_plot <- ggplot(deLosh_data$human_data_linear, aes(x, y)) +
    geom_point(shape=1) + stat_function(fun = linear_function, color = "black") +
  labs(y="Response Magnitude", title="Linear Function",x="") + extrapLines

exponential_plot <- ggplot(deLosh_data$human_data_exp, aes(x, y)) +
  geom_point(aes(shape = "Observed", color = "Observed"),shape=1) + 
  stat_function(aes(color = "True Function"),fun = exponential_function, geom="line")+
  labs(x="Stimulus Magnitude", title="Exponential Function",y="")  +
  extrapLines +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values = c("Observed" = "black", "True Function" = "black")) +
  theme(legend.title = element_blank(), legend.position="top") +
  guides(color = guide_legend(override.aes = list(shape = c(1, NA), 
                                                  linetype = c(0, 1))))

quadratic_plot <- ggplot(deLosh_data$human_data_quad, aes(x = x, y = y)) +
  geom_point( shape = 1) +
  stat_function( fun = quadratic_function, geom = "line") +
  labs(title="Quadratic Function",x="",y="") + extrapLines

linear_plot + exponential_plot + quadratic_plot


```


## ALM Definition

::: column-margin
###### Input Activation

$$
a_i(X)=\exp \left|-\gamma \cdot\left[X-X_i\right]^2\right|
$$

###### Output activation

$$
o_j(X)=\Sigma_{i=1, M} w_{j i} \cdot a_i(X) 
$$

###### Output Probability

$$
P\left[Y_j \mid X\right]=o_j(X) / \Sigma_{k=1, L} o_k(X) 
$$

###### Mean Response

$$
m(X)=\Sigma_{j=1, L} Y_j \cdot P\left[Y_j \mid X\right] 
$$
:::


::: column-margin
###### Feedback Signal

$$
f_j(Z)=e^{-c\cdot(Z-Y_j)^2}
$$

###### Weight Updates

$$
w_{ji}(t+1)=w_{ji}(t)+\alpha \cdot {f_i(Z(t))-O_j(X(t))} \cdot a_i(X(t))
$$
:::





::: column-margin







###### Input node actvation

$$
P[X_i|X] = \frac{a_i(X)}{\\sum_{k=1}^Ma_k(X)}
$$

###### Slope Computation

$$
E[Y|X_i]=m(X_i) + \bigg[\frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1} - X_{i-1}} \bigg]\cdot[X-X_i]
$$
:::






### Generate Response

```{r}
#| code-fold: show
#| code-summary: "Toggle Code"
alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat) {
  input.activation <- exp(-c * (input.layer - input)^2) / sum(exp(-c * (input.layer - input)^2))
  output.activation <- weight.mat %*% input.activation
  output.probability <- output.activation / sum(output.activation)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}
```


### Update Weights Based on Feedback

```{r}
#| code-fold: show
#| code-summary: "Toggle Code"
#| code-overflow: wrap
alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  weight.mat[weight.mat < 0] = 0
  return(weight.mat)
}

alm.trial <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer,output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat))
}
```

### Exam Generalization

```{r}
#| code-fold: show
#| code-summary: "Toggle Code"
exam.response <- function(input, c, trainVec, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer = OUTPUT_LAYER_DEFAULT,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}
```




### Simulation Functions

```{r}
#| code-fold: show

# simulation function
alm.sim <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
  weight.mat <- matrix(0.00, nrow = length(output.layer), ncol = length(input.layer))
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
  dat <- dat %>% mutate(almResp = st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}


simOrganize <- function(simOut) {
  dat <- simOut$d
  weight.mat <- simOut$wm
  c <- simOut$c
  lr <- simOut$lr
  trainX <- unique(dat$x)
  
  almResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model = "ALM", resp = alm.response(x, c, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat = weight.mat)$mean.response)
  
  examResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model = "EXAM", resp = exam.response(x, c, trainVec = trainX, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat))
  
  organized_data <- bind_rows(almResp, examResp) %>% 
    mutate(type = first(dat$type),
           error = abs(resp - y),
           c = c,
           lr = lr,
           type = factor(type, levels = c("linear", "exponential", "quadratic")),
           test_region = ifelse(x %in% trainX, "train", 
                                ifelse(x > min(trainX) & x < max(trainX), "interpolate", "extrapolate")))
  organized_data
}


generateSimData <- function(density, envTypes, noise) {
  reps <- 200 / length(trainingBlocks[[density]])
  map_dfr(envTypes, ~ 
            generate.data(rep(trainingBlocks[[density]], reps), type = .x, noise)) |>
    group_by(type) |>
    mutate(block = rep(1:reps, each = length(trainingBlocks[[density]])),
           trial=seq(1,200))
}

simulateAll <- function(density,envTypes, noise, c = .2, lr = .2) {
  trainMat <- generateSimData(density, envTypes, noise)
  trainData <- map(envTypes, ~ alm.sim(trainMat %>% filter(type == .x), c = c, lr = lr))
  assign(paste(density),list(train=trainData, test=map_dfr(trainData, simOrganize) %>% mutate(density = density)))
}


```

### Simulate Training and Testing

```{r fig.width=10, fig.height=8}
#| column: page-inset-right

envTypes <- c("linear", "exponential", "quadratic")
densities <- c("low", "med", "high")
noise=0
INPUT_LAYER_DEFAULT <- seq(0, 100, 0.5)
OUTPUT_LAYER_DEFAULT <- seq(0, 250, 1)

c = 1.4
lr=.8

results <- map(densities, ~ simulateAll(.x, envTypes, noise, c, lr)) |>
  set_names(densities) 

trainAll <- results %>%
  map_df(~ map_df(.x$train, pluck, "d"), .id = "density") |>
  mutate(stage=as.numeric(cut(trial,breaks=20,labels=seq(1,20))),
         dev=sqrt((y-almResp)^2),
         density=factor(density,levels=c("low","med","high")),
         type=factor(type,levels=c("linear","exponential","quadratic"))) |>
  dplyr::relocate(density,type,stage)

simTestAll <- results |> map("test") |> bind_rows() |>
  group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)


```

### Training Data

```{r fig.width=10, fig.height=8}

trainAll %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
  stat_summary(geom="point",fun=mean,alpha=.4)+
  stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")

```






Predictions for Generalization

```{r fig.width=11,fig.height=10}
#| column: page-inset-right
#| 
simTestAll %>% ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
  facet_grid(density~type) + 
  theme_bw() + theme(legend.position="bottom")


```

Collpasing Across Density Levels gives us:

```{r}
#| eval: false
#| column: page-inset-right
simTestAll %>% group_by(type,model,x,y) %>% summarise(resp=mean(resp))  %>% ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  facet_grid(~type) + 
  theme_bw() + theme(legend.position="bottom")
```



### Table

| Model & Definition                    | R Code                                                                                                      |
|---------------------------------------|-------------------------------------------------------------------------------------------------------------|
| $a_i(X)=\exp \left|-\gamma \cdot\left[X-X_i\right]^2\right|$   | `exp(-c * (input.layer - input)^2)`        |
| $o_j(X)=\Sigma_{i=1, M} w_{j i} \cdot a_i(X)$                 | `weight.mat %*% input.activation`                                                     |
| $P\left[Y_j \mid X\right]=o_j(X) / \Sigma_{k=1, L} o_k(X)$   | `output.activation / sum(output.activation)`                                         |
| $m(X)=\Sigma_{j=1, L} Y_j \cdot P\left[Y_j \mid X\right]$    | `sum(output.layer * output.probability)`                                                  |
| $f_j(Z)=e^{-c\cdot(Z-Y_j)^2}$                                 | `exp(-c * (output.layer - corResp)^2)`                                                               |
| $w_{ji}(t+1)=w_{ji}(t)+\alpha \cdot {f_i(Z(t))-O_j(X(t))} \cdot a_i(X(t)$   | ` lr *(fz - output.activation) %*% t(input.activation)`  |
                                                      |
| $E[Y|X_i]=m(X_i) + \bigg[\frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1} - X_{i-1}} \bigg]\cdot[X-X_i]$ | ` trainVec[which.min(abs(input - trainVec))]; xUnder <- ...; xOver <- ...; mUnder <- ...; mOver <- ...; exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)`|




### Master Function for full simulation

```{r}
#| eval: false
#| include: false 
# Function that goes through every step of generating data, simulating training, and simulating generalization
full.sim <- function(c,lr,noise)
{
  
envTypes <- c("linear", "exponential", "quadratic")
lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
medDensityTrainBlock <- c(
    30.0, 31.5, 33.0, 34.5, 36.5, 38.5, 41.0, 43.5, 46.0,
    48.5, 51.5, 54.0, 56.5, 59.0, 61.5, 63.5, 65.5, 67.0, 68.5, 70.0
)
highDensityTrainBlock <- c(
    30.0, 30.5, 31.0, 32.0, 33.0, 33.5, 34.5, 35.5,
    36.5, 37.0, 38.0, 38.5, 39.5, 40.5, 41.5, 42.0, 43.0,
    43.5, 44.5, 45.5, 46.5, 47.0, 48.0, 48.5, 49.0, 51.0, 51.5, 52.0,
    53.0, 53.5, 54.5, 55.5, 56.5, 57.0, 58.0, 58.5, 59.5, 60.5, 61.5, 
    62.0, 63.0,63.5, 64.5, 65.5, 66.5, 67.0, 68.0, 69.0, 69.5, 70.0
)
# low density has 25 training blocks, medium has 10 blocks, high has 4 blocks. 
# generate training data, for each combination of environment type and density. Use purrr map functions. Rep each dataset by its number of blocks.
a
lowSim <- map(envTypes, ~ alm.sim(lowTrain %>% filter(type == .x), c = 1.4, lr = .4))
medSim <- map(envTypes, ~ alm.sim(medTrain %>% filter(type == .x), c = 1.4, lr = .4))
highSim <- map(envTypes, ~ alm.sim(highTrain %>% filter(type == .x), c = 1.4, lr = .4))

simAll <- rbind(bind_rows(lowSim %>% map("d")) %>% mutate(density = "low"), 
                bind_rows(medSim %>% map("d")) %>% mutate(density = "med"), 
                bind_rows(highSim %>% map("d")) %>% mutate(density = "high"))

simAll <- simAll %>% mutate(stage=as.numeric(cut(trial,breaks=20,labels=seq(1,20))),
                                      dev=sqrt((y-almResp)^2),
                                  #reorder density factor levels
                            density=factor(density,levels=c("low","med","high")),
                            type=factor(type,levels=c("linear","exponential","quadratic"))) %>%
                            dplyr::relocate(density,type,stage)

lowSimTest <- map_dfr(lowSim,simOrganize) %>% mutate(density = "low")
medSimTest <- map_dfr(medSim,simOrganize) %>% mutate(density = "med")
highSimTest <- map_dfr(highSim,simOrganize) %>% mutate(density = "high")

simTestAll <- rbind(lowSimTest,medSimTest,highSimTest) %>% group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)

return(list(simAll=list(simAll),simTestAll=list(simTestAll)))
  
}
```

### Simulations with noise

```{r fig.height=10,fig.width=11}
#| eval: false
#| include: false
k = full.sim(c=1.4,lr=.4,noise=2.0)
k4 = full.sim(c=1.4,lr=.4,noise=4.0)

# run simulation with noise=10, 3 times, average results together. 
k10 = map_dfr(1:3, ~ full.sim(c=1.4,lr=.4,noise=10.0)) %>% group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)



k %>% pluck("simAll") %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
  stat_summary(geom="point",fun=mean,alpha=.4)+
  stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")


k4 %>% pluck("simAll") %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
  stat_summary(geom="point",fun=mean,alpha=.4)+
  stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")


k %>% pluck("simTestAll") %>%ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
  facet_grid(density~type) + 
  theme_bw() + theme(legend.position="bottom")

k4 %>% pluck("simTestAll") %>%ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
  facet_grid(density~type) + 
  theme_bw() + theme(legend.position="bottom")

k10 %>% pluck("simTestAll") %>%ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
  facet_grid(density~type) + 
  theme_bw() + theme(legend.position="bottom")


```
