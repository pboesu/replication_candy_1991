#budworm example from Dennis et al. 1986 / Candy 1991
library(dplyr)

#read raw data which has multiple columns side by side. columns are DDEG (degree days) TOT (total) LSF (life stage) NUM (count)
budworm_raw <- as.matrix(read.table("data/budworm_candy_1991_raw.txt"))
#reshape data
budworm_counts <- as.data.frame(rbind(budworm_raw[,1:4],budworm_raw[,5:8],budworm_raw[,9:12],budworm_raw[,13:16]))
names(budworm_counts) <- c('ddeg','total','stage','count')
budworm_counts %>% arrange(ddeg, stage) -> budworm_counts
readr::write_csv(budworm_counts, 'data/budworm_counts.csv')

budworm_table <- tidyr::pivot_wider(data = budworm_counts, names_from = stage, values_from = count, names_sort = TRUE, names_prefix = 'stage') %>% mutate(ddeg_cent = scale(ddeg, scale = FALSE))
readr::write_csv(budworm_table, 'data/budworm_table.csv')
# write out simplified latex version for manuscript
budworm_table_tex <- budworm_table %>%
  select(-total,-ddeg_cent) %>%
  knitr::kable(format = 'latex', align = 'r', row.names = FALSE,
               col.names = c('Degree Days',paste('Stage',1:7)),
               linesep = c('', '', '', '\\addlinespace'), escape = FALSE, booktabs = TRUE)
cat(budworm_table_tex, file = 'outputs/budworm_table.tex')

budworm_individuals <- budworm_counts %>% tidyr::uncount(count) %>% mutate(ddeg_cent = scale(ddeg, scale = FALSE))
readr::write_csv(budworm_individuals, 'data/budworm_individuals.csv')

