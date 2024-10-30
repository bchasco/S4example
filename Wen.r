library(tidyr)
library(dplyr)

rm(list=ls())


# #Build simple data
locs <- c('CHImark','CHI','WEN','MCJ','JDJ','BON')
states <- c("Sub", "Yr", "Dead")

sheetNames <- c("2021 ChiwCompHist"
                ,"2022 ChiwCompHist"
                ,"Updated2023 ChiwComHist" #had to change a couple of the header names, Event Site Code Value	Event Capture Method Name
                ,"2021 NasCompHist"
                ,"2022 NasCompHist"
                ,"2023 NasCompHist"
                ,"NasRST2021SubsCompHist"
                ,"NasRST2022YCompHist"
                ,"NasRST2022SubsCompHist"
                ,"NasRST2023YCompHist"
                ,"NasRST2023SubsCompHist"
                ,"NasRST2023YCompHist"
                ,"ChiwRST2021SubsCompHist"
                ,"ChiwRST2022YCompHist"
                ,"ChiwRST2022SubsCompHist"
                ,"ChiwRST2023YCompHist"
                ,"ChiwRST2023SubsCompHist"
                ,"ChiwRST2023YCompHist"
                ,"LWenRST2022YCompHist"
                ,"LWenRST2023YCompHist"
                ,"LWenRST2024YCompHist"
)

sheets <- list()
for(i in seq_along(sheetNames)){
  sheets[[i]] <- readxl::read_xlsx('data/Upper Wenatchee 11W Remote Tagging_AllTags.xlsx', sheet = sheetNames[i])
  names(sheets[[i]]) <- names(sheets[[1]])
}

dat <-
  do.call('rbind',sheets)

recode <- read.csv("data/location.csv")

raw <- dat %>%
  # slice(-grep("LWB",dat$`Event Site Code Value`)) %>% #remove barge data
  group_by(`Tag Code`) %>%
  mutate(init_year = min(`Event Year YYYY`)) %>%
  mutate(event_year = `Event Year YYYY`) %>%
  filter((event_year - init_year)<=1) %>% #remove adults
  mutate(event_day = lubridate::yday(`Event Date MMDDYYYY`)) %>%
  mutate(stage = ifelse(event_day<=181,"yr","sub")) %>% #<July 1 is a yearling
  group_by(`Tag Code`,init_year,`Event Site Code Value`) %>%
  mutate(max_event_day = max(event_day)) %>%  #find the max detection for a site if there are dups
  filter(event_day == max_event_day) %>%
  arrange(`Tag Code`,event_year,event_day) %>%
  distinct(`Tag Code`,`Event Year YYYY`,`Event Site Code Value`, event_day, .keep_all = TRUE) %>%
  group_by(`Tag Code`) %>%
  mutate(init_stage = first(stage), init_loc = first(`Event Site Code Value`), init_date = first(`Event Date MMDDYYYY`)) %>%
  mutate(stage = ifelse(init_stage == "yr", "yr", stage)) %>%
  mutate(cum_time = julian(as.Date(`Event Date MMDDYYYY`), origin = as.Date("1970-01-01")) - julian(as.Date(init_date), origin = as.Date("1970-01-01"))) %>%
  filter(cum_time <= 365)

raw$code <- recode$Code[match(raw$`Event Site Code Value`,recode$Location)]
raw$IDnum <- recode$Idnum[match(raw$`Event Site Code Value`,recode$Location)]

raw <- raw %>%
  group_by(`Tag Code`) %>%
  mutate(init_code = first(code), init_num = first(IDnum)) %>%
  filter(IDnum!="")


locs <- c(1,2,4,5,6,7,8,9)
dimnames <- data.frame("Caploc" = c("Trib","First_Trap","Tumwater","Wen","RIS_RIA","MCJ","JDJ","BON","TWX_EST"),
                       "dimID" = c(1,2,3,4,5,6,7,8,9))
pvt <- raw %>%
  mutate(stageID = ifelse(stage == "sub", 1, 2)) %>%
  filter(init_num%in%(locs)) %>% #be suer you don't take capture histories for fish tagged upstream
  select(`Tag Code`, stageID, IDnum, code, init_year) %>%
  left_join(dimnames, by = c("IDnum" = "dimID")) %>%
  mutate(loc = ifelse(!is.na(Caploc), Caploc, IDnum)) %>%
  filter(loc%in%dimnames$Caploc[dimnames$dimID%in%locs]) %>% #be suer you don't take capture histories for fish tagged upstream
  mutate(loc = factor(loc, levels = dimnames$Caploc)) %>%
  group_by(`Tag Code`, loc, init_year) %>%
  summarise(stageID = first(stageID)) %>%
  pivot_wider(names_from = loc, values_from = stageID, values_fill = 3) %>%
  select('Tag Code', init_year, dimnames$Caploc[dimnames$dimID%in%locs]) %>%
  ungroup() #%>%
  # filter(Trib == 1)


ch <- pvt%>%
  # select(-(dimnames$Caploc[!(dimnames$dimID%in%locs)]))
  mutate(ch = apply(select(., dimnames$Caploc[dimnames$dimID%in%locs]),1,paste, collapse = "")) %>%
  group_by(ch, init_year) %>%
  # group_by(ch) %>%
  summarise(n = n())

ch_i <- as.data.frame(t(sapply(ch$ch,function(x){as.integer((strsplit(x,split = ""))[[1]])})))
# ch_i[ch_i==2] <- 1
# ch_i[ch_i==3] <- 2

ch_i <- ch_i[t(apply(ch_i,1,function(x){sum(x<3)}))!=0,]
U <- apply(ch_i,1,function(x){min(which(x<3,arr.ind = TRUE))})
names(ch_i) <- dimnames$Caploc[dimnames$dimID%in%locs]
init_loc <- apply(ch_i,1,function(x){min(which(x<3,arr.ind = TRUE))})
init_state <- apply(ch_i,1,function(x){min(x[x<3])})
ch_i <- data.frame(ch_i,n = ch$n, Year = ch$init_year, init = init_loc, init_loc = init_loc, init_state = init_state)
ch_i$y_i <- as.integer(factor(ch_i$Year))


MR_settings <- list(
  data = na.omit(ch_i),
  locs = dimnames$Caploc[dimnames$dimID%in%locs],
  states = states,
  state_frm = ' ~ LH',
  model_type = "MS",
  frms = list(
    phi = 'int ~  1' #Survival
    ,p = 'int ~  1' #Detection probability
    ,lam = 'int ~ 1'  #joint last location
    ,eta =  'int ~ (1|loc)' #transition probability
  ),
  make_sd = TRUE
)


tmb_list <- new("tmb_list",
                MR_settings = MR_settings)

tmb_obj <- make_tmb_lists(tmb_list)

st <- Sys.time()
tmb_model <- build_tmb(tmb_obj)

tmb_model@TMB$rep$gamma[1,2,,]
tmb_model@TMB$rep$gamma[2,2,,]

tmb_model@TMB$rep$phi[,1,]
tmb_model@TMB$rep$phi[,2,]

# print(round(plogis(tmb_model@TMB$opt$par),2))
# print(tmb_model@TMB$sd)
# print(apply(tmb_model@TMB$rep$gamma,c(2,3,4),mean)[2,,])
# print(Sys.time()- st)
#
