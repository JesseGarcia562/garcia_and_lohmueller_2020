library(tidyverse)
library(glue)
library(kableExtra)


compute_d_prime<-function(ld_df){

  ld_df<-ld_df %>% mutate(d_prime= case_when(
  d < 0 ~ d / pmin(pA*pB, (1-pA) *(1-pB)),
  d > 0 ~ d / pmin(pA*(1-pB), (1-pA) * pB  ),
  d == 0 ~ 0))
  return(ld_df)
  
}


low_recomb<-read_rds("high_frequency_ld_recombination_rate_1e-09_all.rds")


ld_df<-compute_d_prime(low_recomb)






set.seed(1)




low_recomb_summarised<-ld_df %>%
  mutate(AC.x=parse_number(AC.x)) %>%
  group_by(seed) %>%
  group_nest() %>%
  sample_n(3000) %>%
  mutate(group=1:3000) %>%
  unnest(cols=data)

low_recomb_summarised$AC_break<-cut(low_recomb_summarised$AC.x, c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))



df<-low_recomb_summarised %>% ungroup() %>%
  group_by( AC_break, Variation) %>%
  summarise(mean_d_prime=mean(d_prime), N=length(d_prime), standard_error=sd(d_prime)/sqrt(N))

df<-df %>% mutate(mean_d_prime=format(round(mean_d_prime,3),  nsmall=3))
df<-df %>% mutate(standard_error=formatC(standard_error, format="e", digits=2))


pretty_df<-bind_cols( df %>% filter(Variation == "Nonsynonymous"), df %>% filter(Variation == "Synonymous")) %>%
  select(-Variation...2, -Variation...7, -AC_break...6) %>%
  select(AC_break...1, N...4, mean_d_prime...3, standard_error...5, N...9, mean_d_prime...8, standard_error...10) %>%
  rename(`Allele count`=AC_break...1, 
         `Number of pairs` = N...4,
         Mean = mean_d_prime...3, 
         `Standard error` = standard_error...5, 
         `Number of pairs ` = N...9,
         'Mean ' = mean_d_prime...8,
         `Standard error ` = standard_error...10)

kable(pretty_df, format = "latex", booktabs=T) %>% 
  add_header_above(c(" "=1 , "Nonsynonymous Pairs" =3, "Synonymous Pairs" = 3)) %>% 
  kable_styling(latex_options = c("striped"), full_width = T, font_size = 6 ) %>%
  as_image(file="table_high_frequency_ld_low_recomb.png")







avg_recomb<-read_rds("high_frequency_ld_recombination_rate_1e-08_all.rds")

ld_df<-compute_d_prime(avg_recomb)


avg_recomb_summarised<-ld_df %>%
  mutate(AC.x=parse_number(AC.x)) %>%
  group_by(seed) %>%
  group_nest() %>%
  sample_n(3000)%>%
  unnest(cols=data)

avg_recomb_summarised$AC_break<-cut(avg_recomb_summarised$AC.x,c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))

  
df<-avg_recomb_summarised %>% ungroup() %>%
  group_by( AC_break, Variation) %>%
  summarise(mean_d_prime=mean(d_prime), N=length(d_prime), standard_error=sd(d_prime)/sqrt(N))


df<-df %>% mutate(mean_d_prime=format(round(mean_d_prime,3),  nsmall=3))
df<-df %>% mutate(standard_error=formatC(standard_error, format="e", digits=2))





pretty_df<-bind_cols( df %>% filter(Variation == "Nonsynonymous"), df %>% filter(Variation == "Synonymous")) %>%
  select(-Variation...2, -Variation...7, -AC_break...6) %>%
  select(AC_break...1, N...4, mean_d_prime...3, standard_error...5, N...9, mean_d_prime...8, standard_error...10) %>%
  rename(`Allele count`=AC_break...1, 
         `Number of pairs` = N...4,
         Mean = mean_d_prime...3, 
         `Standard error` = standard_error...5, 
         `Number of pairs ` = N...9,
         'Mean ' = mean_d_prime...8,
         `Standard error ` = standard_error...10)

kable(pretty_df, format = "latex", booktabs=T) %>% 
  add_header_above(c(" "=1 , "Nonsynonymous Pairs" =3, "Synonymous Pairs" = 3)) %>% 
  kable_styling(latex_options = c("striped"), full_width = T, font_size = 6 ) %>%
  as_image(file="table_high_frequency_ld_avg_recomb.png")

