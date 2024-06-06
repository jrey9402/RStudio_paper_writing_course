# load packages of general use across the project --------------
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(readr)
library(ggpubr)
library(rstudioapi)
library(magick)
library(janitor)
library(tinytable)
library(rio)
library(ggplot2)
library(multcompView)
library(graphics)
library(readxl)
library (forcats)
library(ggthemes)
library(wesanderson)
#Import datasets
data_Jose1 <- rio::import("Z:/José Rey CC/!!Projects/III.NCED3/V.NGA/Peels experiments/miRNANGA1_4/5_220124_ koPeels with Cys_/220124_data_cys.xlsx")
data_Jose2 <- rio::import("Z:/José Rey CC/!!Projects/III.NCED3/V.NGA/Peels experiments/miRNANGA1_4/7weeks 140224/140224_data_cys.xlsx")

#Export data to gitrepo
rio::export(data_Jose1, "analysis/data/Jose/nga4ko_5weeks_cys.csv")
rio::export(data_Jose2, "analysis/data/Jose/nga4ko_7weeks_cys.csv")

#ANOVA1-----------

model1 = aov(length ~ Treatment+genotype + Treatment*genotype, data=data_Jose1)

#table with factors, means and standard deviation 
data_summ1 = group_by(data_Jose1, genotype, Treatment) %>% 
  summarise(mean=mean(length), sd=sd(length)) %>% 
  arrange(desc(mean))
#View(data_sum)

#tukey test1

Tukey1= TukeyHSD(model1)
print(Tukey1)

Cld1=multcompLetters4(model1, Tukey1)
print(Cld1)

#add letter to table
Cld1= as.data.frame.list(Cld1$`Treatment:genotype`)
data_summ1$Cld1=Cld1$Letters
print(data_summ1)

# Plot data - Jose 
plot_Jose1= ggplot(data_Jose1, aes(x=forcats::fct_relevel(Treatment,"Control","Cys","ABA"), y=length, fill = Treatment,palette="Spectral"))+
  geom_boxplot(outlier.shape = NA) +facet_grid(.~forcats::fct_relevel(genotype,"WT","nga4ko","aba3")) +
  theme_base()+ 
  scale_y_continuous(limits = c(0, 15.5))+
  guides(fill = guide_legend(title = "Treatment"))+
geom_text(data = data_summ1, aes(label=Cld1,x=Treatment, y=mean +sd), position=position_dodge2(0), vjust =-5) +
  labs(x="Treatments", y = "Stomatal aperture (µm)") +
  theme(legend.title = element_text(colour="black", size=10, 
                                    face="bold"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))



#ANOVA2-------
  
  model2 = aov(length ~ Treatment+genotype + Treatment*genotype, data=data_Jose2)
  
  #table with factors, means and standard deviation 
  data_summ2 = group_by(data_Jose2, genotype, Treatment) %>% 
    summarise(mean=mean(length), sd=sd(length)) %>% 
    arrange(desc(mean))
  #View(data_sum)
  
  #tukey test2
  
  Tukey2= TukeyHSD(model2)
  print(Tukey2)
  
  Cld2=multcompLetters4(model2, Tukey2)
  print(Cld1)
  
  #add letter to table
  Cld2= as.data.frame.list(Cld2$`Treatment:genotype`)
  data_summ2$Cld2=Cld2$Letters
  print(data_summ2)
  
  # Plot data - Jose 
  plot_Jose2= ggplot(data_Jose2, aes(x=forcats::fct_relevel(Treatment,"Control","Cys","ABA"), y=length, fill = Treatment,palette="Spectral"))+
    geom_boxplot(outlier.shape = NA) +facet_grid(.~forcats::fct_relevel(genotype,"WT","nga4ko","aba3")) +
    theme_base()+ 
    scale_y_continuous(limits = c(0, 15.5))+
    guides(fill = guide_legend(title = "Treatment"))+
    geom_text(data = data_summ2, aes(label=Cld2,x=Treatment, y=mean +sd), position=position_dodge2(0), vjust =-5) +
    labs(x="Treatments", y = "Stomatal aperture (µm)") +
    theme(legend.title = element_text(colour="black", size=10, 
                                      face="bold"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))


  #style plots
  args(theme)
  
  theme_plots <- theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.key.size = unit(7, "mm"),
      legend.title.position = "top",
      legend.background = element_rect(color = "grey"),
      plot.title.position = "panel"
    )
  plot_Jose1 <- plot_Jose1 +
    theme_plots
  plot_Jose1
  
  plot_Jose2 <- plot_Jose2 +
    theme_plots
  plot_Jose2

#Saving figures----
  ggsave( "analysis/data/Jose/pictures/plot_Jose1.png",
          limitsize = FALSE,
          units = c("px"), plot_Jose1,
          width = 1600, height = 1400
  )

  ggsave( "analysis/data/Jose/pictures/plot_Jose2.png",
          limitsize = FALSE,
          units = c("px"), plot_Jose2,
          width = 1600, height = 1400
  )
  
#Panels creation----  
  img1 <- readPNG("analysis/data/Jose/pictures/plot_Jose1.png")
  img2 <- readPNG("analysis/data/Jose/pictures/plot_Jose2.png")
  
  #convert to panels
  panel_JoseA <- ggdraw() + draw_image(img1)
  panel_JoseB <- ggdraw() + draw_image(img2)  

#define layout with textual representation
  layout <- "
AB"
  
#assemble multipanel figure based on layout
  Figure_Jose <- panel_JoseA + panel_JoseB + plot_Jose1 + plot_Jose2 +
    plot_layout(design = layout, heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 12, face='plain'))

  ggsave(
    "analysis/data/Jose/figures/Figure_Jose.png", limitsize = FALSE, 
    units = c("px"), Figure_Jose, width = 4000, height = 1600,
    bg = "white"
  )
  