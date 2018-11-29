### data wrangling for Marginal Zone compartmnet
rm(list = ls())
gc()
setwd("~/Desktop/Git_repos/GDT_dynamics")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

library(tidyverse)

# source for FM 
parent_pop <- "Th.immature"
target_pop <- "naive"

# total counts and donor fractions for the source poppulation
source_donor <- readxl::read_excel(path = "data/gdT_numbers.xlsx", sheet = 5) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains(parent_pop))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_host <- readxl::read_excel(path = "data/gdT_numbers.xlsx", sheet = 6) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains(parent_pop))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
# calculating total counts, donor fractions
source_join <- full_join(source_host, source_donor, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor")) %>%
  mutate(total_counts = TH.immature.host + TH.immature.donor,
         fd = (TH.immature.donor)/total_counts)%>%
  select(-contains(".host"), -contains(".donor"))

# total counts and donor fractions for the source poppulation
gdt_donor <- readxl::read_excel(path = "data/gdT_numbers.xlsx", sheet = 5) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains(target_pop))%>% 
  mutate(Total_naive = SP.naive + LN.naive)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

gdt_host <- readxl::read_excel(path = "data/gdT_numbers.xlsx", sheet = 6) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains(target_pop))%>% 
  mutate(Total_naive = SP.naive + LN.naive)%>%  select(-contains("SP"), -contains("LN"))%>% 
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
# calculating total counts, donor fractions 
gdt_join <- full_join(gdt_host, gdt_donor, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_naive.host + Total_naive.donor,
         fd = Total_naive.donor/ total_counts)%>%
  select(-contains(".host"), -contains(".donor"))

# normalising donor fraction in FM by dividing with the donor fractions in the source compartment
gdt_fd <- gdt_join %>%
  full_join(source_join, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".gdt", ".source"))%>%
  mutate(Nfd = fd.gdt/ fd.source) %>%
  filter(Nfd <= 1.5)%>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("Nfd"))


gdt_counts <- gdt_join%>%
  right_join(gdt_fd, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"))%>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("total_counts"))


#plots
# for gdt niave
ggplot(gdt_counts) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(5e4, 1e6), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(60, 500), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = paste("Total counts: ", target_pop, sep = "")) +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(gdt_fd) +
  geom_point(aes(x = time.post.BMT, y = Nfd),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+
  ylim(0, 1.2) + theme_bw() + labs(x = "Days post t0", y = NULL, title = paste("Noirmalised donor fractions: ", target_pop)) +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

# for source i.e. immature cells in Thymus
ggplot(source_join) +
  geom_point(aes(x = age.at.S1K, y = total_counts), col=4, size = 4)+
  scale_y_log10(limits = c(5e3, 5e5), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 500), trans = "log10", breaks = c(300, 100, 30)) + 
  theme_bw() + labs(x = "Host age (days)", y = NULL, title = "Total counts: T2") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

ggplot(source_join) +
  geom_point(aes(x = time.post.BMT, y = fd),  size = 4)+
  geom_hline(yintercept = 1.00, linetype = 2, size =1.5, col="darkred")+
  ylim(0, 1.2) + theme_bw() + labs(x = "Days post t0", y = NULL, title = "Noirmalised donor fractions: SP FM") +
  theme(axis.text = element_text(size = 22),
        axis.title =  element_text(size = 20, face = "bold"),
        plot.title = element_text(size=20,  hjust = 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20, face = "bold"))

# writting files to the data directory 
write.csv(gdt_counts, "~/Desktop/Git_repos/GDT_dynamics/data/counts_gdt.csv")
write.csv(gdt_fd, "~/Desktop/Git_repos/GDT_dynamics/data/Nfd_gdt.csv")
write.csv(source_join, "~/Desktop/Git_repos/GDT_dynamics/data/source_gdt.csv")





