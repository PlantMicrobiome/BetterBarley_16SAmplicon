#Title:"Barley domestication age affects growth and nutrient uptake at low nutrient availability"

#R-script for plant and soil chemical and physical data -


# Plant data: Scripts for analyses of nested linear models and figures for shoot dry weight (Figure 2 in Results),
# leaf nitrogen (Figure 3 in Results), and other leaf macro- and micronutrients (Figure 4 in Results).

# Soil data: linear regression of C/N, pH, NH4, and NO3 in organic fertilizer treatments (Figure 5 in Results) 

#########################################
##### SHOOT DRY WEIGHT LINEAR MODEL #####
#########################################

# 1) Import excel file: "Data_Shoot_Dw_N"

# 2) Test effect of fertilizer treatment (Treatment), domestication age (Age) and Cultivar (Cultivar) on shoot dry weight (Dw_mg)
# using nested linear mixed model (lme4 package) with nesting of Cultivar within Age:   

library(lme4)

Data_Shoot_Dw_N$Treatment <- factor(Data_Shoot_Dw_N$Treatment)
Data_Shoot_Dw_N$id <- factor(c(1:nrow(Data_Shoot_Dw_N)))

Data_Shoot_Dw_N$Cultivar <- factor(Data_Shoot_Dw_N$Cultivar, levels = c("Babushka","Langeland", "Salka", "Irina KWS", "Feedway", "Flair", "RGT Planet"))
Data_Shoot_Dw_N$Age <- factor(Data_Shoot_Dw_N$Age, levels = c("Old", "Modern"))
Data_Shoot_Dw_N$Treatment <- factor(Data_Shoot_Dw_N$Treatment, levels = c("No_fertilizer", "Organic_1", "Organic_2", "Organic_3", "Organic_4", "Organic_6", "Mineral"))

model_Dw <- lmer(Dw_mg ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar) +
             (0 + as.numeric(Treatment=="Mineral") +
                as.numeric(Treatment=="No_fertilizer") +
                as.numeric(Treatment=="Organic_1") +
                as.numeric(Treatment=="Organic_2") +
                as.numeric(Treatment=="Organic_3") +
                as.numeric(Treatment=="Organic_4") +
                as.numeric(Treatment=="Organic_6") || id),
           data = Data_Shoot_Dw_N,
           na.action = na.omit,
           control=lmerControl(check.nobs.vs.nlev = "ignore",
                               check.nobs.vs.rankZ = "ignore",
                               check.nobs.vs.nRE="ignore")
)# 3) Inspect normality of model residuals using Q-Q and residual plot: 

residuals <- residuals(model_Dw, type = "pearson")

plot(residuals)

qqnorm(residuals)
abline(0,sd(residuals))

# 4) Use joint_tests function in emmeans package to produce analysis-of variance-like table (Type III SS) based on the linear model (Dw_LM_0):

library(emmeans)

(jt_Dw <-joint_tests(model_Dw))

# 5) Since we see a significant interaction between fertilizer treatment (Treatment) and domestication age (Age) we use the emmeans function in the emmeans
# package to calculate estimated marginal means (EMMs) based on the nested linear mixed model and do Tukey HSD on the EMMs, comparing each age group
# within each fertilizer treatment:  

Emmeans_Dw_LM_0 <- emmeans(model_Dw,~Age|Treatment)

#Visualize the results of the Tukey HSD using the pwpm function in emmeans: 

(pwpm_Dw <-pwpm(Emmeans_Dw_LM_0))

# 6) However, we also see that the effect of Domestication age depends on specific cultivar (Age:Cultivar).
# Therefore, we compare cultivars within each fertilizer treatment (have to use linear model to do this): 

m1 <- lmer(Dw_mg ~ Treatment + Cultivar + Treatment:Cultivar +
             (0 + as.numeric(Treatment=="Mineral") +
                as.numeric(Treatment=="No_fertilizer") +
                as.numeric(Treatment=="Organic_1") +
                as.numeric(Treatment=="Organic_2") +
                as.numeric(Treatment=="Organic_3") +
                as.numeric(Treatment=="Organic_4") +
                as.numeric(Treatment=="Organic_6") || id),
           data = Data_Shoot_Dw_N,
           na.action = na.omit,
           control=lmerControl(check.nobs.vs.nlev = "ignore",
                               check.nobs.vs.rankZ = "ignore",
                               check.nobs.vs.nRE="ignore")
)


emm1_m1  <-emmeans(m1,~Cultivar|Treatment, type="response")

#Visualize the results of the Tukey HSD using the pwpm function in emmeans: 
(pwpm_Cultivar <-pwpm(emm1_m1))


######################################
##### LEAF NITROGEN LINEAR MODEL #####
######################################

# 1) Test effect of fertilizer treatment (Treatment), domestication age (Age) and Cultivar (Cultivar) on leaf nitrogen (Shoot_N)
# using linear model with nesting of Cultivar within Age (Shoot_N is log transformed as it gives better Q-Q and residual plot):   

N_LM_0 <- lm(log(Shoot_N) ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_Shoot_Dw_N)

# 2) Inspect normality of model residuals using Q-Q and residual plot: 

qqnorm(residuals(N_LM_0))
abline(0,sd(residuals(N_LM_0)))

res <- resid(N_LM_0)
plot(N_LM_0$fitted.values, res)

# 3) Use joint_tests function in emmeans package to produce analysis-of variance-like table (Type III SS) based on the linear model (N_LM_0):

library(emmeans)
(jt_N <-joint_tests(N_LM_0))

# 4) Since we see a significant effect of Age and not Cultivar (Age*Cultivar) we use the emmeans function in the emmeans
# package to calculate estimated marginal means (EMMs) based on the linear model and do Tukey HSD on the EMMs, comparing each
# dommestication age (old and modern) within each fertilizer treatment:  

Emmeans_N_LM_0 <- emmeans(N_LM_0,~Age|Treatment)

#Visualize the results of the Tukey HSD using the pwpm function in emmeans: 

(pwpm_Dw <-pwpm(Emmeans_N_LM_0))


###############################################
##### BOX PLOT (SHOOT DRY WEIGHT) Fig. 2A #####
###############################################

library(ggplot2)

(Dw_box <- ggplot(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Dw_mg)) +
    geom_boxplot(aes(fill = Age), outlier.shape = NA, lwd = 1) +
    scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB")) +
    geom_point(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Dw_mg, group = Age), color = "#000000", position = position_dodge(width = 0.75), size = 4) +
    geom_point(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Dw_mg, color = Cultivar, group = Age), position = position_dodge(width = 0.75), size = 2) +
    scale_color_manual(values = c("Babushka" = "#FF9800", "Langeland" = "#ED6167", "Salka" = "#878787", "Irina KWS" = "#E5E5E5", "Feedway" = "#90EAE4", "Flair" = "#AA3377", "RGT Planet" = "#74648c")) +
    theme_classic() +
    scale_x_discrete(labels = c("No_fertilizer" = "No fertilizer", "Organic_1" = "Organic 1", "Organic_2" = "Organic 2", "Organic_3" = "Organic 3", "Organic_4" = "Organic 4", "Organic_6" = "Organic 6", "Mineral" = "Mineral")) +
    scale_y_continuous(name = expression(bold("Shoot dry weight (mg plant"^{-1}*")"))) +
    
    theme(
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 5), face="bold"),
      axis.text.y = element_text(size = 14, color = "black", margin = margin(l = 5, r = 5), face="bold"),
      legend.text = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 10),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")  # Adjust top margin value as per your requirement
    ) +
    guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)))


ggsave('Dw_box.alt.color.pdf',Dw_box,width=10 ,height=6)
ggsave('Dw_box.alt.color.png',Dw_box,width=10,height=6, type="cairo")



####################################
##### BOX PLOT (leaf N) Fig. 2B ####
####################################

(N_box <- ggplot(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Shoot_N)) +
   geom_boxplot(aes(fill = Age), outlier.shape = NA, lwd = 1) +
   scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB")) +
   geom_point(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Shoot_N, group = Age), color = "#000000", position = position_dodge(width = 0.75), size = 4) +
   geom_point(data = Data_Shoot_Dw_N, aes(x = Treatment, y = Shoot_N, color = Cultivar, group = Age), position = position_dodge(width = 0.75), size = 2) +
   scale_color_manual(values = c("Babushka" = "#FF9800", "Langeland" = "#ED6167", "Salka" = "#878787", "Irina KWS" = "#E5E5E5", "Feedway" = "#90EAE4", "Flair" = "#AA3377", "RGT Planet" = "#74648c")) +
   geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
   theme_classic() +
   scale_x_discrete(labels = c("No_fertilizer" = "No fertilizer", "Organic_1" = "Organic 1", "Organic_2" = "Organic 2", "Organic_3" = "Organic 3", "Organic_4" = "Organic 4", "Organic_6" = "Organic 6", "Mineral" = "Mineral")) +
   scale_y_continuous(name = "Leaf N (% of dry weight)") +
   theme(
     axis.title = element_text(size = 16, face = "bold"),
     axis.title.x = element_blank(),
     axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 5), face="bold"),
     axis.text.y = element_text(size = 14, color = "black", margin = margin(l = 5, r = 5), face="bold"),
     legend.text = element_text(size = 13, face = "bold"),
     legend.title = element_text(size = 10),
     axis.line = element_line(linewidth = 1),
     axis.ticks = element_line(linewidth = 1),
     legend.position = "none",
     plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")  # Adjust top margin value as per your requirement
   ) +
   guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1)))


ggsave('N_box.alt.color.pdf',N_box,width=10 ,height=6)
ggsave('N_box.alt.color.png',N_box,width=10,height=6, type="cairo")


###################################################
##### MACRO- AND MICRONUTRIENTS LINEAR MODELS #####
###################################################

# 1) Import excel file: "Data_ICP"

# 2) Sort Data_ICP:

Data_ICP$Cultivar <- factor(Data_ICP$Cultivar, levels = c("Babushka","Langeland", "Salka", "Irina KWS", "Feedway", "Flair", "RGT Planet"))
Data_ICP$Age <- factor(Data_ICP$Age, levels = c("Old", "Modern"))

# 3) Test effects of fertilizer treatment (Treatment), domestication age (Age) and Cultivar (Cultivar) on
# macronutrients: P, K, S, Ca, Mg and micronutrients: Fe, Zn, Mn, and B, using linear models with nesting of
# Cultivar within Age. We use the same procedure as described above for Shoot dry weight and leaf nitrogen:

# Phosphorus (P) (log transformed for better model fit)

P_LM_0 <- lm(log(P) ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(P_LM_0))
abline(0,sd(residuals(P_LM_0)))

# Get list of residuals
res <- resid(P_LM_0)
# Produce residual vs. fitted plot
plot(P_LM_0$fitted.values, res)

(jt_P <- joint_tests(P_LM_0))

# Significant interaction between Treatment*Age, therefore we use emmeans to compare domestication age
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_P_LM_0 <- emmeans(P_LM_0,~Age|Treatment, type="response") 

(pwpm_P <-pwpm(Emmeans_P_LM_0))

# Potassium (K)

K_LM_0 <- lm(K ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(K_LM_0))
abline(0,sd(residuals(K_LM_0)))

# Get list of residuals
res <- resid(K_LM_0)
# Produce residual vs. fitted plot
plot(K_LM_0$fitted.values, res)

(jt_K <- joint_tests(K_LM_0))

# Significant effect of Age, therefore we use emmeans to compare domestication age 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_K_LM_0 <- emmeans(K_LM_0,~Age|Treatment)

(pwpm_K <-pwpm(Emmeans_K_LM_0))

# Sulfur (S)

S_LM_0 <- lm(S ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(S_LM_0))
abline(0,sd(residuals(S_LM_0)))

# Get list of residuals
res <- resid(S_LM_0)
# Produce residual vs. fitted plot
plot(K_LM_0$fitted.values, res)

(jt_S <-joint_tests(S_LM_0))

# Significant effect of Age, therefore we use emmeans to compare domestication age 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_S_LM_0 <- emmeans(S_LM_0,~Age|Treatment) 

pwpm_S<-pwpm(Emmeans_S_LM_0)

# Calcium (Ca) (log transformed for better model fit)

Ca_LM_0 <- lm(log(Ca) ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Ca_LM_0))
abline(0,sd(residuals(Ca_LM_0)))

# Get list of residuals
res <- resid(Ca_LM_0)
# Produce residual vs. fitted plot
plot(Ca_LM_0$fitted.values, res)

(jt_Ca <-joint_tests(Ca_LM_0))

# Significant effect of Age:Cultivar, therefore we use emmeans to compare cultivar 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral): 

Emmeans_Ca_LM_0 <- emmeans(Ca_LM_0,~Cultivar|Treatment, type="response") 

pwpm_Ca<-pwpm(Emmeans_Ca_LM_0)

# Magnesium (Mg)

Mg_LM_0 <- lm(Mg ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Mg_LM_0))
abline(0,sd(residuals(Mg_LM_0)))

# Get list of residuals
res <- resid(Mg_LM_0)
# Produce residual vs. fitted plot
plot(Mg_LM_0$fitted.values, res)

(jt_Mg <-joint_tests(Mg_LM_0))

# Only significant effect of fertilizer treatment, therefore use emmeans to compare fertilizer treatments: 

Emmeans_Mg_LM_0 <- emmeans(Mg_LM_0,~Treatment) 

(pwpm_Mg<-pwpm(Emmeans_Mg_LM_0))

# Iron (Fe)

Fe_LM_0 <- lm(Fe ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Fe_LM_0))
abline(0,sd(residuals(Fe_LM_0)))

# Get list of residuals
res <- resid(Fe_LM_0)
# Produce residual vs. fitted plot
plot(Fe_LM_0$fitted.values, res)

(jt_Fe <-joint_tests(Fe_LM_0))

# Significant effect of Age, therefore we use emmeans to compare domestication age 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):  

Emmeans_Fe_LM_0 <- emmeans(Fe_LM_0,~Age|Treatment) 

(pwpm_Fe <-pwpm(Emmeans_Fe_LM_0))

# Zinc (Zn) (log transformed for better model fit)

Zn_LM_0 <- lm(log(Zn) ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Zn_LM_0))
abline(0,sd(residuals(Zn_LM_0)))

# Get list of residuals
res <- resid(Zn_LM_0)
# Produce residual vs. fitted plot
plot(Zn_LM_0$fitted.values, res)

(jt_Zn <-joint_tests(Zn_LM_0))

# Significant effect of Age, therefore we use emmeans to compare domestication age 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_Zn_LM_0 <- emmeans(Zn_LM_0,~Age|Treatment, type="response") 

(pwpm_Zn <-pwpm(Emmeans_Zn_LM_0))

# Copper (Cu)

Cu_LM_0 <- lm(Cu ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Cu_LM_0))
abline(0,sd(residuals(Cu_LM_0)))

# Get list of residuals
res <- resid(Cu_LM_0)
# Produce residual vs. fitted plot
plot(Cu_LM_0$fitted.values, res)

(jt_Cu <-joint_tests(Cu_LM_0))

# Significant effect of Age, therefore we use emmeans to compare domestication age 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_Cu_LM_0 <- emmeans(Cu_LM_0,~Age|Treatment) 

pwpm(Emmeans_Cu_LM_0)

# Boron (B)(log transformed for better model fit)

B_LM_0 <- lm(log(B) ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(B_LM_0))
abline(0,sd(residuals(B_LM_0)))

# Get list of residuals
res <- resid(B_LM_0)
# Produce residual vs. fitted plot
plot(B_LM_0$fitted.values, res)


(jt_B <-joint_tests(B_LM_0))

# Significant effect of Cultivar, therefore we use emmeans to compare Cultivar 
# within fertilizer treatments (No fertilizer, Organic 4, Mineral):

Emmeans_B_LM_0 <- emmeans(B_LM_0,~Cultivar|Treatment, type="response") 

(pwpm_B <-pwpm(Emmeans_B_LM_0))

# Manganese (Mn)

Mn_LM_0 <- lm(Mn ~ Treatment + Age/Cultivar + Treatment:(Age/Cultivar), data = Data_ICP)

qqnorm(residuals(Mn_LM_0))
abline(0,sd(residuals(Mn_LM_0)))

# Get list of residuals
res <- resid(Mn_LM_0)
# Produce residual vs. fitted plot
plot(Mn_LM_0$fitted.values, res)

(jt_Mn <-joint_tests(Mn_LM_0))

# Only significant effect of fertilizer treatment, therefore we use emmeans to compare fertilizer treatments:

Emmeans_Mn_LM_0 <- emmeans(Mn_LM_0,~Treatment, type="response") 

pwpm(Emmeans_Mn_LM_0)

###########################################################################################
##### BOX PLOT OF NUTRIENTS AFFECTED BY DOMESTICATION AGE: K, P, S, Fe, Zn, Cu Fig. 3 #####
###########################################################################################

#K boxplot
(K_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=K/10000))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=K/10000, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=K/10000, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name="Leaf K (% of dry weight)", limits = c(0,7))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10), face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('K_box.alt.color.pdf',K_box,width=10 ,height=6)
ggsave('K_box.alt.color.png',K_box,width=10,height=6, type="cairo")

#P boxplot
(P_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=P/10000))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=P/10000, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=P/10000, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name="Leaf P (% of dry weight)", limits = c(0.15,0.45))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10),face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('P_box.alt.color.pdf',P_box,width=10 ,height=6)
ggsave('P_box.alt.color.png',P_box,width=10,height=6, type="cairo")


#S boxplot
(S_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=S/10000))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=S/10000, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=S/10000, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name="Leaf S (% of dry weight)", limits = c(0,0.7))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10), face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('S_box.alt.color.pdf',S_box,width=10 ,height=6)
ggsave('S_box.alt.color.png',S_box,width=10,height=6, type="cairo")


#Fe boxplot
(Fe_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=Fe))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Fe, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Fe, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 100, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name=expression(bold(Leaf~Fe~(ug~g^-1~dry~weight))), limits = c(25,125))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10), face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('Fe_box.alt.color.pdf',Fe_box,width=10 ,height=6)
ggsave('Fe_box.alt.color.png',Fe_box,width=10,height=6, type="cairo")

#Zn boxplot
(Zn_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=Zn))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Zn, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Zn, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 20, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name=expression(bold(Leaf~Zn~(ug~g^-1~dry~weight))), limits = c(10,80))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10), face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('Zn_box.alt.color.pdf',Zn_box,width=10 ,height=6)
ggsave('Zn_box.alt.color.png',Zn_box,width=10,height=6, type="cairo")

#Cu boxplot
(Cu_box <- ggplot(data = Data_ICP, aes(x=Treatment,y=Cu))+
  geom_boxplot(aes(fill=Age),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("Old" = "#FFFFFF", "Modern" = "#0BA1CB"))+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Cu, group=Age),color="#000000", position = position_dodge(width = 0.75), size=5)+
  geom_point(data=Data_ICP, aes(x=Treatment, y=Cu, color=Cultivar, group=Age), position = position_dodge(width = 0.75), size=3)+
  scale_color_manual(values = c("Babushka"="#FF9800","Langeland"="#ED6167","Salka"="#878787","Irina KWS"="#E5E5E5","Feedway"="#90EAE4", "Flair"="#AA3377", "RGT Planet"="#74648c"))+
  geom_hline(yintercept = 6, linetype = "dashed", color = "#CBC4AA", linewidth=1)+
  theme_classic()+
  scale_x_discrete(labels = c("a"="No fertilizer", "e"="Organic 4", "h"="Mineral Fertilizer"))+
  scale_y_continuous(name=expression(bold(Leaf~Cu~(ug~g^-1~dry~weight))), limits = c(4,10))+ 
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 16, color = "black", margin=margin(t=10), face="bold"))+
  theme(axis.text.y = element_text(size = 20, color = "black", margin=margin(r=10), face="bold"))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1)))

ggsave('Cu_box.alt.color.pdf',Cu_box,width=10 ,height=6)
ggsave('Cu_box.alt.color.png',Cu_box,width=10,height=6, type="cairo")

####################################################################################################################
##### LINEAR REGRESSION OF SOIL PARAMETERS: C:N, pH, and NH4-N + NO3-N IN ORGANIC FERTILIZER TREATMENTS Fig. 4 #####
####################################################################################################################

# 1) Import excel file: "Data_Soil_CN_NH4NO3_pH"

# Remove calculated C:N (Type = 2) and mineral fertilizer treatment as it does not fit linearly with the other treatments 
Data_1 <-subset(Data_Soil_CN_NH4NO3_pH, Type!="2" & Treatment!="h")
na_Data_1 <- na.omit(Data_1)

# Isolate calculated C:N (Based on measured C and N in soil and organic fertilizer at experiment start)
Data_2 <-subset(Data_Soil_CN_NH4NO3_pH, Type=="2" & Code!= "cal6") 

library(ggplot2)
library(ggpmisc)

# Soil NH4 and NO3 regressions

#linear model NH4
lm_NH4 <- lm(log1p(NH4) ~ Added_fert_g, data=na_Data_1)
summary(lm_NH4)

qqnorm(residuals(lm_NH4))
abline(0,sd(residuals(lm_NH4)))

res <- resid(lm_NH4)
plot(lm_NH4$fitted.values, res)

#linear model NO3
lm_NO3 <- lm(log(NO3+1) ~ Added_fert_g, data=na_Data_1)
summary(lm_NO3)

qqnorm(residuals(lm_NO3))
abline(0,sd(residuals(lm_NO3)))

res <- resid(lm_NO3)
plot(lm_NO3$fitted.values, res)

x.axis.label <- expression(bold("Organic fertilizer (g pot"^"-1"~")"))
y.axis.label <- expression(bold("Soil NH"[4]^"+"*"-N and NO"[3]^"-"*"-N (ug g"^-1* "dry weight)"))

NH4_NO3_plants_reg <- ggplot(na_Data_1, aes(x = Added_fert_g)) +
  geom_smooth(aes(y = NH4), method = 'lm', formula = y ~ x, fullrange = TRUE, color = "#CE9C69", fill = "#CE9C69") +
  geom_smooth(aes(y = NO3), method = 'lm', formula = y ~ x, fullrange = TRUE, color = "#666666", fill = "#666666") +
  scale_y_continuous(
    trans = 'log1p',
    breaks = c(0, 50, 100, 150),
    labels = c(0, 50, 100, 150),
    limits = c(0, 160),
    expand = c(0, 0)
  ) +
  geom_point(aes(y = NH4), size = 5, shape = 20, color = "#000000") +
  geom_point(aes(y = NH4), size = 3, shape = 20, color = "#CE9C69") +
  geom_point(aes(y = NO3), size = 5, shape = 20, color = "#000000") +
  geom_point(aes(y = NO3), size = 3, shape = 20, color = "#666666") +
  
  xlim(0, 30) +
  scale_fill_manual(values = c("#666666", "#CE9C69")) +
  theme(legend.position = "bottom") +
  theme_minimal() +
  xlab(x.axis.label) +
  ylab(y.axis.label) +
  theme(axis.text.x = element_text(size = 18, color = "black",face="bold")) +
  theme(axis.text.y = element_text(size = 18, color = "black",face="bold")) +
  theme(axis.title.x = element_text(size = 20, color = "black",face="bold")) +
  theme(axis.title.y = element_text(size = 20, color = "black", face="bold")) +
  theme(legend.text = element_text(size = 12, color = "black")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

ggsave('NH4_NO3_plants_reg.pdf',NH4_NO3_plants_reg,width=8,height=8)
ggsave('NH4_NO3_plants_reg.png',NH4_NO3_plants_reg,width=8,height=8,type="cairo")


# Soil C:N regression

# Isolate calculated C:N (calculations based on measured C and N in soil and organic fertilizer at experiment start)
Data_2 <-subset(Data_Soil_CN_NH4NO3_pH, Type=="2") 

# Fit the linear models
lm_CN <- lm(C_N_end ~ Added_fert_g, data = Data_1)
summary(lm_CN)

x.axis.label <- expression(bold("Organic fertilizer (g pot"^"-1"~")"))
y.axis.label="Soil C:N" #<-- define the label of x axis
(CN_plants_reg=ggplot(Data_1,aes(x=Added_fert_g,y=C_N_end)) +
    geom_smooth(method='lm',fullrange = TRUE, color="#666666", fill="#666666") +
    geom_point(colour = "#000000", size = 3) +
    geom_point(colour = "#666666",size = 2) +
    #stat_poly_eq(formula=y ~ x,aes(label=paste(..eq.label.., ..rr.label.., ..p.value.label..,sep = "~~~")),parse=TRUE,label.x.npc = "left",label.y.npc = "bottom",rr.digits=2, size=6, color="#666666", vjust = -5, hjust = -0.1)+
    geom_point(data=Data_2, aes(x=Added_fert_g, y=C_N_end), color="black", size = 3)+
    geom_point(data=Data_2, aes(x=Added_fert_g, y=C_N_end), color="red", size = 2)+
    xlim(0, 30)+
    theme_minimal()+
    xlab(x.axis.label)+
    ylab(y.axis.label)+
    theme(axis.text.x=element_text(size=18,color="black", face="bold"))+
    theme(axis.text.y=element_text(size=18,color="black", face="bold"))+
    theme(axis.title.x=element_text(face="bold",size=20,color="black"))+
    theme(axis.title.y=element_text(face="bold", size=20,color="black"))+
    theme(legend.text=element_text(size=12,color="black"))+
    theme(legend.title=element_blank()))

ggsave('CN_plants_reg.pdf',CN_plants_reg,width=8,height=8)
ggsave('CN_plants_reg.png',CN_plants_reg,width=8,height=8,type="cairo")

x.axis.label <- expression(bold("Organic fertilizer (g pot"^"-1"~")"))
y.axis.label=expression(bold("Soil pH (0.01 M CaCl"["2"]*")")) #<-- define the label of x axis
(pH_reg=ggplot(Data_1,aes(x=Added_fert_g,y=Soil_pH)) +
    geom_smooth(method='lm',fullrange = TRUE, color="#666666", fill="#666666") +
    geom_point(colour = "#000000", size = 3) +
    geom_point(colour = "#666666",size = 2) +
    #stat_poly_eq(formula=y ~ x,aes(label=paste(..eq.label.., ..rr.label.., ..p.value.label..,sep = "~~~")),parse=TRUE,label.x.npc = "right",label.y.npc = "top",rr.digits=3, size=6, color="#666666")+
    xlim(0, 30)+
    ylim(4,8)+
    theme_minimal()+
    xlab(x.axis.label)+
    ylab(y.axis.label)+
    theme(axis.text.x=element_text(size=18,color="black", face="bold"))+
    theme(axis.text.y=element_text(size=18,color="black", face="bold"))+
    theme(axis.title.x=element_text(face="bold",size=20,color="black"))+
    theme(axis.title.y=element_text(face="bold", size=20,color="black"))+
    theme(legend.text=element_text(size=12,color="black"))+
    theme(legend.title=element_blank()))

ggsave('pH_plants_reg.pdf',pH_reg,width=8,height=8)
ggsave('pH_plants_reg.png',pH_reg,width=8,height=8,type="cairo")



