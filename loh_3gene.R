



#######FIGURE 1B#######

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)

setwd("~/Desktop/loh_paper/from_220301/num_of_loh")
loh <- read.table("~/Desktop/loh_paper/from_220301/identified_iloh/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid.txt")
loh_count <- loh %>% 
  group_by(V6) %>% 
  summarise(n=n())

mean(loh_count$n)
median(loh_count$n)

ggplot(loh_count, aes(x = reorder(V6, -n), y = n)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 40.66667, linetype = "dashed", color = "red") + #mean
  geom_hline(yintercept = 39, linetype = "dashed", color = "blue") + #median
  scale_y_continuous(breaks = seq(0,110,20))+
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.4, colour = "gray80"),
        panel.grid.minor.x = element_line(size = 0.15, colour = "gray80"))

