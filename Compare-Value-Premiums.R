# compares different value premiums 2022 05

# ====

path = 'D:/OneDrive/openap/Release_2022/CrossSection/Portfolios/Data/Portfolios/'

signalselect = 'BMdec'

csvlist = dir(path, pattern = 'PredictorAlt', full.names = T)
temp = dir(path, pattern = 'PredictorAlt', full.names = F) 

style = substr(temp, 19,nchar(temp)-4)


alldat = tibble()
for (i in 1:length(csvlist)){
  print(i)
  temp = fread(csvlist[[i]]) %>% 
    filter(signalname == signalselect, port == 'LS') %>% 
    select(signalname, date,ret) %>% 
    mutate(
      style = style[[i]]
    )
  
  alldat = rbind(alldat,temp)
}



# simple table
alldat %>% group_by(style) %>% 
  summarize(rbar = mean(ret, na.rm= T))

# break down by recent sample
alldat %>% 
  mutate(
    recent = year(date) >= 2005
  ) %>% 
  group_by(style, recent) %>% 
  summarize(rbar = mean(ret, na.rm= T)) %>% 
  pivot_wider(
    names_from = recent, values_from = rbar
  )