



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


#######FIGURE 2A and Figure 2-figure supplement 1#######
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)
setwd("~/Desktop/loh_paper/from_220301/corr_age/")
loh_age <- read.table("~/Desktop/loh_paper/from_220301/identified_iloh/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_header.txt"
                      , header = T)

ge12 <- read.csv("~/Desktop/loh_paper/from_220301/raw_data/age_info/second_gen.dnms.summary.csv", header = T)

ge12_edit <- ge12[, c(5,10,11,12,19,21)]
names(ge12_edit) <- c("d1_age", "fam", "m1_id", "m1_age", "d1_id", "ge2_id")

loh_aged <- loh_age %>% 
  left_join(ge12_edit,
            by = c('sam_id'='ge2_id'))

#paternal
ggscatter(loh_aged_count, x = "d1_age", y = "n",
          fullrange = TRUE,
          conf.int = FALSE, xlim = c(15, 50),
          ylab = "Number of I-LOHs", xlab = "Paternal age",
          font.x = 15, font.y = 15,
          font.xtickslab = 14, font.ytickslab = 14,
          ggtheme = theme_bw()) + 
  stat_cor(label.x = 41) +
  theme(aspect.ratio = 1)

#maternal
ggscatter(loh_aged_count, x = "m1_age", y = "n",
          fullrange = TRUE,
          conf.int = FALSE, xlim = c(15, 50),
          ylab = "Number of I-LOHs", xlab = "Maternal age",
          font.x = 15, font.y = 15,
          font.xtickslab = 14, font.ytickslab = 14,
          ggtheme = theme_bw()) + 
  stat_cor(label.x = 41) +
  theme(aspect.ratio = 1)

#paternal#Figure 2
ggscatter(loh_aged_count, x = "d1_age", y = "n",
          fullrange = TRUE,
          conf.int = FALSE, xlim = c(15, 50),
          ylab = "Number of I-LOHs", xlab = "Paternal age",
          font.x = 15, font.y = 15,
          font.xtickslab = 14, font.ytickslab = 14,
          ggtheme = theme_bw()) + 
  stat_cor(label.x = 41) +
  theme(aspect.ratio = 1)

#maternal#Figure 2-figure supplement 1
ggscatter(loh_aged_count, x = "m1_age", y = "n",
          fullrange = TRUE,
          conf.int = FALSE, xlim = c(15, 50),
          ylab = "Number of I-LOHs", xlab = "Maternal age",
          font.x = 15, font.y = 15,
          font.xtickslab = 14, font.ytickslab = 14,
          ggtheme = theme_bw()) + 
  stat_cor(label.x = 41) +
  theme(aspect.ratio = 1)


#######FIGURE 2C#######
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)

setwd("~/Desktop/loh_paper/from_220301/loh_origin/")
loh <- read.table("~/Desktop/loh_paper/from_220301/identified_iloh/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid.txt")
ori <- read.table("~/Desktop/loh_paper/from_220301/loh_origin/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_d2m2_v10_isnvori.txt")

names(loh) <- c("chrm", "start", "end", "fam", "sex", "ge2_id", "het", "isnv", "isnv_count", "min_length", "max_length")
names(ori) <- c("chrm","start", "end", "ref", "alt", "fam", "sex", "loh_info", "ori" )

loh_ori <- loh %>%
  left_join(ori, by = c("chrm", "start", "fam", "sex"))

write.table(loh_ori, "~/Desktop/loh_paper/from_220301/loh_origin/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_header_aged_ori.txt", 
            quote = F, sep = "\t", row.names = F)

loh_ori_donor <- read.table("~/Desktop/loh_paper/from_220301/loh_origin/allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_header_aged_ori_edit_uniq_donor.txt") 

loh_ori_edit_count_donor <- loh_ori_donor %>% 
  group_by(V6,V16) %>% 
  summarise(n=n())

ggboxplot(loh_ori_edit_count_donor, x = "V16", y = "n",
          color = "V16",
          palette = c("dodgerblue", "deeppink"),
          add = "jitter", ylab = "Number of donor", xlab = "donor",
          legend = "none",
          font.x = 15, font.y = 15,
          font.xtickslab = 14, font.ytickslab = 14,  
          ggtheme = theme_bw()) + stat_compare_means(label.x = 1.3)


#######FIGURE 3B#######
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)

setwd("~/Desktop/loh_paper/from_220301/germ_somatic/")
loh_case <- read.table("allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_sorted_intersect_extract_v2_filtering_trimmed_case.txt")
loh_case_edit <- loh_case[, c(1,2,3,6,12)]

loh_case_edit$V12 <- gsub("^case1$", "case13", loh_case_edit$V12)
loh_case_edit$V12 <- gsub("^case3$", "case13", loh_case_edit$V12)
loh_case_edit$V12 <- gsub("^case2$", "case24", loh_case_edit$V12)
loh_case_edit$V12 <- gsub("^case4$", "case24", loh_case_edit$V12)

loh_case_edit_count <- loh_case_edit %>% 
  group_by(V6,V12) %>% 
  summarise(n=n())

#back to back graph

library(magrittr)
library(ggpol)
library(dplyr)

loh_case_edit_count_reorder <- loh_case_edit_count %>% 
  mutate(V6 = reorder(V6, -n, mean)) %>% 
  mutate(n_mod = ifelse(V12 == "case13", -n, n)) 

loh_case_edit_count_reorder_xmax <- loh_case_edit_count_reorder %>% 
  mutate(x_max = ifelse(V12 == "case13", -45, 45))


ggplot(loh_case_edit_count_reorder_xmax, aes(x = n_mod, y = reorder(V6, -n), fill = V12)) + 
  geom_col() +
  facet_share(~ V12, dir = "h",scales = "free",  reverse_num = TRUE) +
  labs(x = "Number of LOH", y = "Sample ID") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c("red1", "goldenrod1")) +
  geom_blank(aes(x = x_max)) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray60", size = 0.2),
        panel.grid.minor.x = element_line(color = "gray70", size = 0.1),
        panel.grid.major.y = element_line(color = "gray70", size = 0.1),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10))


#######FIGURE 4A#######
