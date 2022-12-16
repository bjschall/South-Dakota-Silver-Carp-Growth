###########################################
##########  Load in Libraries  ############
###########################################
library(Matrix)
library(tidyverse)
library(brms)
library(readxl)
library(tidybayes)

###########################################
#############  Load in data  ##############
###########################################
df <- read_excel("SVC_Final_Ages.xlsx", sheet = "Sheet1")
df <- df %>% filter(Age>1)

#create function showing previous growth equation from Hayer et al. 2014
myfun <- function(Age) {
  1223.709 * (1 - exp(-0.173*(Age-0)))
}
###########################################
########  Create the von B Equation  ######
###########################################
form_1 <- bf(TL ~ Linf *(1- exp(-K*(Age-t0))), #Von B Equation
             Linf ~ Waterbody,
             K ~ Waterbody,
             t0 ~ Waterbody,
             nl = TRUE)

###########################################
##########  Prior Simulations  ############
###########################################
#set up von B values and SDs
priors <- tibble(Linf = abs(rnorm(200, 900, 100)),
                K = abs(rnorm(200, .2, 0.1)),
                t0 = rnorm(200, 0, 1),
                iter = 1:200)

#create dataframe with potential values
prior_sims <- priors %>% 
  expand_grid(df %>% distinct(Age) %>% 
                bind_rows(tibble(Age = c(1:14)))) %>% 
  mutate(length_sims = Linf *(1- exp(-K*(Age-t0))))

#plot the curve
ggplot() + 
  geom_line(data=prior_sims, aes(x = Age, y = length_sims, group = iter))

###############################################
##########  Growth model by water  ############
###############################################
#get priors
get_prior(form_1, data=df)

#run the model
model <- brm(form_1, data = df,
             prior = c(prior(normal(900, 100), nlpar = "Linf"),
                       prior(normal(0, 50), nlpar = "Linf", coef="WaterbodyJames"),
                       prior(normal(.2, 0.1), nlpar = "K",lb=0),
                       prior(normal(0, 0.2), nlpar = "K", coef="WaterbodyJames"),
                       prior(normal(0, 1), nlpar = "t0"),
                       prior(normal(0, 1), nlpar = "t0", coef="WaterbodyJames"),
                       prior(exponential(0.25), class="sigma")),
             file = "SilverCarpVB.rds",
             family = gaussian(), 
             chains = 4, cores=4, iter=2000)

#check the output
print(model, digits=3)
pp_check(model, ndraws = 50)
plot(model)
bayes_R2(model)

###############################################
##########  Posterior predictions  ############
###############################################
#predict mean and CrI across observed ages
post_preds <- tibble(Age = seq(1, max(df$Age))) %>% 
  expand_grid(Waterbody = c("BigSioux", "James")) %>% 
  add_epred_draws(model, re_formula = NULL)

#summarize and see results
posts <- post_preds %>% 
  group_by(Age, Waterbody) %>% 
  mean_qi(.epred); posts

#predictive interval across observed ages
post_pred_int <- tibble(Age = seq(1, max(df$Age))) %>% 
  expand_grid(Waterbody = c("BigSioux", "James")) %>% 
  add_predicted_draws(model)

pred_int <- post_pred_int %>%  
  group_by(Age, Waterbody) %>% 
  mean_qi(.prediction); pred_int


###############################################
############# Plotting Posterior ##############
###############################################
posts %>% 
  ggplot(aes(x = Age, y = .epred)) + 
  geom_line(aes(color = Waterbody), size = 1) +
  geom_point(data=df, aes(y=TL, x=Age, fill=Waterbody),colour="black",pch=21, size = 2.5, position = position_jitter(width=.05))+
  ylab("Total length (mm)")+
  xlab("Age (years)")+
  scale_fill_manual(values=c("black","grey70"))+
  scale_color_manual(values=c("black","grey70"))+
  scale_x_continuous(breaks=seq(1:14))+
  scale_y_continuous(limits = c(0, 1200), breaks=seq(0, 1200, by=100))+
  #geom_ribbon(aes(fill = Waterbody, ymin = .lower, ymax = .upper, y=.epred), alpha = 0.3)+
  geom_ribbon(data=pred_int, aes(fill = Waterbody, ymin = .lower, ymax = .upper, y=.prediction), alpha = 0.3)+
  theme_classic()+
  theme(legend.position = c(.9, 0.3),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=18),
        axis.text = element_text(size = 14))+
  stat_function(fun = myfun, geom = "line", size=1, linetype="dotted")+
  guides(fill="legend")

#ggsave("Figure 1. Silver Carp Growth Curves.jpeg", plot=last_plot(),width=8, heigh=5, units="in", dpi=600)

###############################################
########## Probability Statements #############
###############################################
###Linf
Linf_vals <- model %>% 
  spread_draws(b_Linf_Intercept, b_Linf_WaterbodyJames) %>% 
  rename(Linf_BigSioux = "b_Linf_Intercept") %>% 
  mutate(Linf_James = Linf_BigSioux + b_Linf_WaterbodyJames); Linf_vals

#probability the Big Sioux Linf value is less than the value in Hayer et al. (2014)
sum(Linf_vals$Linf_BigSioux<1223)/4000

#probability the James Linf value is less than the value in Hayer et al. (2014)
sum(Linf_vals$Linf_James<1223)/4000

#probability the James Linf value is greater than the Big Sioux value
sum(Linf_vals$b_Linf_WaterbodyJames>0)/4000


###K
K_vals <- model %>% 
  spread_draws(b_K_Intercept, b_K_WaterbodyJames) %>% 
  rename(K_BigSioux = "b_K_Intercept") %>% 
  mutate(K_James = K_BigSioux + b_K_WaterbodyJames); K_vals

#probability the Big Sioux K value is greater than the value in Hayer et al. (2014)
sum(Linf_vals$Linf_BigSioux>0.173)/4000

#probability the James K value is greater than the value in Hayer et al. (2014)
sum(Linf_vals$Linf_James>0.173)/4000

#probability the James K value is greater than the Big Sioux value
sum(K_vals$b_K_WaterbodyJames>0)/4000


#####
###probability that ages 1-6 are larger from Hayer et al. 2014
#Big Sioux
sum(filter(post_pred_int, Age==1, Waterbody=="BigSioux")$.prediction<myfun(1))/4000
sum(filter(post_pred_int, Age==2, Waterbody=="BigSioux")$.prediction<myfun(2))/4000
sum(filter(post_pred_int, Age==3, Waterbody=="BigSioux")$.prediction<myfun(3))/4000
sum(filter(post_pred_int, Age==4, Waterbody=="BigSioux")$.prediction<myfun(4))/4000
sum(filter(post_pred_int, Age==5, Waterbody=="BigSioux")$.prediction<myfun(5))/4000
sum(filter(post_pred_int, Age==6, Waterbody=="BigSioux")$.prediction<myfun(6))/4000

#James
sum(filter(post_pred_int, Age==1, Waterbody=="James")$.prediction<myfun(1))/4000
sum(filter(post_pred_int, Age==2, Waterbody=="James")$.prediction<myfun(2))/4000
sum(filter(post_pred_int, Age==3, Waterbody=="James")$.prediction<myfun(3))/4000
sum(filter(post_pred_int, Age==4, Waterbody=="James")$.prediction<myfun(4))/4000
sum(filter(post_pred_int, Age==5, Waterbody=="James")$.prediction<myfun(5))/4000
sum(filter(post_pred_int, Age==6, Waterbody=="James")$.prediction<myfun(6))/4000

###############################################
##########  Sensitivity Analysis  #############
###############################################
#sensitivity model
sens_model <- brm(form_1, data = df,
             prior = c(prior(normal(900, 200), nlpar = "Linf"),
                       prior(normal(0, 100), nlpar = "Linf", coef="WaterbodyJames"),
                       prior(normal(.2, 0.2), nlpar = "K",lb=0),
                       prior(normal(0, 0.4), nlpar = "K", coef="WaterbodyJames"),
                       prior(normal(0, 2), nlpar = "t0"),
                       prior(normal(0, 2), nlpar = "t0", coef="WaterbodyJames"),
                       prior(exponential(0.125), class="sigma")),
             file = "SilverCarpVB_sens.rds",
             family = gaussian(), 
             chains = 4, cores=4, iter=2000)

#check the output
print(sens_model, digits=3)
pp_check(sens_model, ndraws = 50)
plot(sens_model)
bayes_R2(sens_model)

#compare models
loo(sens_model, model)
