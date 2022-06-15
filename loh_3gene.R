



#######FIGURE 1B#######
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)

setwd("SET YOUR WORKING DIRECTORY")
loh <- read.table("raw_identified_all_loh.txt")
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
setwd("SET YOUR WORKING DIRECTORY")
loh_age <- read.table("raw_identified_all_loh.txt, header = T)

ge12 <- read.csv("second_gen.dnms.summary.csv", header = T) #YOU CAN GET THE "second_gen.dnms.summary.csv" from GitHub of Sasani, 2019

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

setwd("SET YOUR WORKING DIRECTORY")
loh <- read.table("raw_identified_all_loh.txt")
ori <- read.table("raw_identified_all_origin_of_loh.txt")

names(loh) <- c("chrm", "start", "end", "fam", "sex", "ge2_id", "het", "isnv", "isnv_count", "min_length", "max_length")
names(ori) <- c("chrm","start", "end", "ref", "alt", "fam", "sex", "loh_info", "ori" )

loh_ori <- loh %>%
  left_join(ori, by = c("chrm", "start", "fam", "sex"))


loh_ori_donor <- read.table("raw_identified_all_origin_of_loh.txt REARRANGED BY "donor"") 

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
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)

setwd("~/Desktop/loh_paper/from_220301/gc_content/")
loh_gc <- read.table("allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_sorted_iloh.txt")

loh_gc$V4 <- round(loh_gc$V4, 2)
loh_gc_more100 <- subset(loh_gc, V5 >= 100)

loh_gc_more100$grade <- cut(loh_gc_more100$V4, breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3,
                                                          0.35,0.4,0.45,0.5,0.55,0.6,
                                                          0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))


loh_gc_more100_count <- loh_gc_more100 %>% 
  group_by(grade) %>% 
  summarise(n=n())

loh_gc_more100_count_freq <- loh_gc_more100_count %>% 
  mutate(freq_bin = n/sum(n))

loh_gc_more100_count_loh <- aggregate(loh_gc_more100$V6,
                                      by=list(grade=loh_gc_more100$grade), FUN=sum)


loh_gc_more100_count_loh_freq <- loh_gc_more100_count_loh %>% 
  mutate(freq = x/sum(x))
  
loh_join <- loh_gc_more100_count_freq %>% 
  full_join(loh_gc_more100_count_loh_freq, by = c("grade"))

loh_join_freq <- mutate(loh_join, obs_exp = freq/freq_bin)

loh_join_freq_subset <- subset(loh_join_freq, obs_exp > 0)

loh_join_freq_subset_edit <- loh_join_freq_subset %>% 
  mutate(ratio = ifelse(obs_exp < 1 , -obs_exp, obs_exp))

loh_join_freq_subset_edit_edit <- loh_join_freq_subset_edit %>% 
  mutate(ratio_munusone = ratio - 1)
loh_join_freq_subset_edit_edit[, "grade_edit"] <- c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75)

 ggplot(loh_join_freq_subset_edit_edit, aes(x = grade, y = obs_exp - base )) +
     geom_bar(stat = "identity") +
     scale_y_continuous(breaks = c(-1, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,1  ), labels = function(obs_exp) obs_exp + base) +
     coord_cartesian(ylim = c(-1, 1)) +
     theme_minimal() +
     theme(axis.text.y = element_text(size= 14),
             axis.text.x = element_text(size = 14,angle = 40, hjust = 1))
   
       
#######FIGURE 4B#######
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggExtra)
library(reshape2)
library(ggrepel)

setwd("~/Desktop/loh_paper/from_220301/chr_territory/")
loh_terri <- read.table("allfam_phasing_loh_onlypm_nowhatdnv_v5_v16_range_reshaped_v10_isnv_updown_iloh_samid_sorted_iloh_len_by_chrom.bed")
loh_terri_melt <- melt(data=loh_terri, id.vars=1:4)
names(loh_terri) <- c("sam_id", "chr", "loh_count", "loh_len", "count_per", "len_per")
loh_terri_mean <- loh_terri %>% 
  group_by(chr) %>% 
  summarise(avg = mean(count_per), med = median(count_per))


gene_dens <- read.table("~/Desktop/loh_paper/from_211201/chr_territory/gencode.v19.annotation_protein_coding_gene_noLCR_norepeat_uniq_count_gene_per.txt")
names(gene_dens) <- c("chr", "dens")

terri_gene_dens <- loh_terri_mean %>% 
  full_join(gene_dens, by = c("chr"))

library(ggrepel)
ggplot(terri_gene_dens, aes(x = avg, y = dens)) +
  geom_point() +
  geom_text_repel(aes(label = chr),size = 3)
                        
gene_dens <- read.table("~/Desktop/loh_paper/from_211201/chr_territory/gencode.v19.annotation_protein_coding_gene_noLCR_norepeat_uniq_count_gene_per_1M.txt")
names(gene_dens) <- c("chr", "dens")

terri_gene_dens <- loh_terri_mean %>% 
  full_join(gene_dens, by = c("chr"))

ggplot(terri_gene_dens, aes(x = avg, y = dens)) +
  geom_point(color = "red", size = 2.5) +
  geom_text_repel(aes(label = chr),size = 5.5) +
  scale_y_continuous(breaks = c(0,10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_x_continuous(breaks = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.12)) +
  coord_cartesian(ylim = c(0,80), xlim = c(0.01, 0.12)) +
  theme_minimal() + 
  theme(axis.text.y = element_text(size= 16),
        axis.text.x = element_text(size = 16))
                        
                        
                      
                        
                      
