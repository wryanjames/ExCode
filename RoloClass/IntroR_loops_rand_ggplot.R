# """ Review of vectors and matrices, for loops, random sampling, and ggplot
#     author: 'Advanced Ecology: Populations and Communities'
#     date: 10/14/2021 """

library(tidyverse)
library(truncnorm)


## Vector and matrix operations in R

#Vectors are one the primary data structures in R and can be made in different ways depending on the type of data stored. Math can be done on numeric vectors. 

  

# make a numeric vector
a = c(1.1,5,3,4)
a

# make a integer vector
b = 1:15
b

# make a character vector 
c = c('a', 'b', 'c')
c 



# summary statistics of sequence of numbers
d = rnorm(n = 40, mean = 6, sd = 2)
mean(d) #mean
median(d) #median
sd(d) #standard deviation
quantile(d, c(0.025, 0.075)) # 95% CI

# math with vectors
e = c(3, 1, 5, 6)
e*4
e*e

 

# Another data structure that is useful in ecological uses of R are matrices. A matrix is made with the `matrix()` function with the basic syntax `matrix(data = , nrow = , ncol = , byrow = , dimnames = )`. 
# 
# `data =` the input vector which becomes the data elements of the matrix.
# 
# `nrow =` the number of rows to be created.
# 
# `ncol =` the number of columns to be created.
# 
# `byrow = FALSE` if `TRUE` then the input vector elements are arranged by row.
# 
# `dimname =` the names assigned to the rows and columns.

  
# make a matrix using nrow
m1 = matrix(c(1,2,3,4,5,6), nrow = 2)
m1 

# make a matrix using ncol
m2 = matrix(c(1,2,3,4,5,6), ncol = 2)
m2 

# make a matrix using nrow and by row
m3 = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = T)
m3

 

# Matrices can be indexed in 2 ways

  
# make a matrix using ncol
m2 = matrix(c(1,2,3,4,5,6), ncol = 2)
m2 

m2[1,2]
m2[4]
# make a matrix using nrow and by row
m3 = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = T)
m3

m3[1,2]
m3[4]

 

# Matrix algebra can be done on matrices. Also functions like `sum()`, `rowSums()`, and `colSums()`,`mean()`, `rowMeans()`, and `colMeans()`. 

  

m2 = matrix(c(1,2,3,4,5,6), ncol = 2)
m2 

m2 + 1 
m2/3

m3 = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = T)
m3

m2 + m3 
m3 - m2

# multply matrix with vector
v = c(1,2,3)

m2*v

# functions with matrices
m2
sum(m2)
mean(m2)
rowSums(m2)
rowMeans(m2)
colSums(m2)
colMeans(m2)

 

## `for` loops
# Another useful tool in programming is `for` loops. `for` loops repeat a process for a certain number of iterations. These can be useful iterate over a dataset or when using information in a time series. The `for` loop works over the number sequence indicated and does the code within the loop (inside of `{}`) for each number in the sequence. The iteration is typically indicated with `i`, but is just an object that is replaced at the begining of each loop and can be anything.

  
for(i in 1:10){
  print(i)
}

for(turtle in 5:10){
  print(turtle)
}

for(flower in 1:nrow(iris)){
  cat('The species for this iteration is ',
      as.character(iris$Species[flower]), '\n')
}


b = 1:10
for (i in 2:10){
  z = b[i] - b[i-1]
  
  cat('z =', z, 'b[i] =', b[i], 'b[i-1] =', b[i-1], '\n')
}


 

## Random sampling in R
# Generating random numbers or randomly sampling is commonly required in ecology. R has base functions that can be used to randomly generate numbers depending on the distribution wanted. Also packages like [truncnorm](https://www.rdocumentation.org/packages/truncnorm/versions/1.0-8/topics/truncnorm) and  `dplyr` in `tidyverse` are useful.  


#### Generate random numbers
  
# base R
rnorm(5, mean = 0, sd = 1)
runif(5, min = 0, max = 5)
rbinom(5, size = 2, prob = 0.5)
rbinom(5, size = 2, prob = 0.2)

 

  
#library(rtruncnorm)
# truncnorm to truncate normal distribution
rnorm(20, mean = 1, sd = 2)
truncnorm::rtruncnorm(20, mean = 1, sd = 2, a=-Inf, b=Inf)
truncnorm::rtruncnorm(20, mean = 1, sd = 2, a=0, b=2)
 

#### Random selection 
  
# base R 
b = 1:10
sample(b, size =2, replace = F)

c = 1:5 
sample(c, size = 6, replace = T)

# from list
sample(c('good', 'bad'), size = 8, replace = T)

# change probability
sample(c('good', 'bad'), size = 8, replace = T, prob = c(0.2, 0.8))

# sample matrices from list 
m1 = matrix(c(6,5,4,3,2,1), ncol = 2)
m2 = matrix(c(1,2,3,4,5,6), ncol = 2)
m3 = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = T)

m = list(m1, m2, m3)

sample(m, size = 4, replace = T)

sample(m, size = 4, replace = T, prob = c(0.8, 0.1, 0.1))

 


  
#library(tidyverse)
df = tibble(cond = c('good', 'bad', 'ok'), prob = c(0.5, 0.3, 0.2))

dplyr::sample_n(df, size = 2, replace = F)

dplyr::sample_n(df, size = 3, replace = T)

dplyr::sample_n(df, size = 10, replace = T, weight = prob)
 

# We can bring what we learned together to simulate change in a population over time


  

year_type = tibble(cond = c('good', 'average', 'bad'),
                   birth = c(500, 300, 100),
                   death = c(100, 300, 500),
                   prob = c(0.25, 0.5, 0.25))

pop_start = 50000 
time_length = 100
reps = 10


for (i in 1:reps){
  trial = paste0('trial_',i)
  pop = tibble(time = 0:time_length, n = NA,
               trial = trial, year = NA, change = NA)
  pop$n[pop$time == 0] = pop_start
  
  for (t in 1:(nrow(pop)-1)){
    year = sample_n(year_type, size = 1, weight = prob)
    birth = round(rnorm(n = 1, mean = year$birth[1], sd = 50))
    death = round(rnorm(n = 1, mean = year$death[1], sd = 50))
    pop$n[pop$time == t] = birth - death + pop$n[pop$time == (t-1)]
    pop$year[pop$time == t] = year$cond[1]
    pop$change[pop$time == t] = birth - death
  }
  
  if (i ==1){
    tpop = pop
  }else{
    tpop = bind_rows(tpop,pop)
  }
}

head(tpop)


## Figures with `ggplot2` 
# The `ggplot2` package is part of the packages that load with `tidyverse` and has become the standard in ecology. The syntax builds upon on a base function and is very customizable [see cheat sheet](https://www.rstudio.com/resources/cheatsheets/). 
# 
# The base of all `ggplot2` begins with `ggplot()` and `geom_...()` are built upon them 


ggplot(tpop, aes(x = time, y = n))+
  geom_point()
 

# Show color based on trial number and line connecting dots

  
ggplot(tpop, aes(x = time, y = n, color = trial))+
  geom_point()+
  geom_line()
 

# Change labels and style of plot


  
ggplot(tpop, aes(x = time, y = n, color = trial))+
  geom_point()+
  geom_line()+
  labs(x = 'Time', y = 'Population size')+
  theme_classic()
 

# Modify the size of axis label text and legend position  

  
ggplot(tpop, aes(x = time, y = n, color = trial))+
  geom_point()+
  geom_line()+
  labs(x = 'Time', y = 'Population size')+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_blank())
 

# Modify the range of time on x axis

  
ggplot(tpop, aes(x = time, y = n, color = trial))+
  geom_point()+
  geom_line()+
  scale_x_continuous(limits = c(0,20))+
  labs(x = 'Time', y = 'Population size')+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_blank())
 

# Split each trial into own grid

  
ggplot(tpop, aes(x = time, y = n))+
  geom_point()+
  geom_line()+
  labs(x = 'Time', y = 'Population size')+
  facet_wrap(~trial)+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_blank())
 

# Modify the label and size of strip text


  
labels = c('trial_1' = 'Trial 1',
           'trial_2' = 'Trial 2',
           'trial_3' = 'Trial 3',
           'trial_4' = 'Trial 4',
           'trial_5' = 'Trial 5',
           'trial_6' = 'Trial 6',
           'trial_7' = 'Trial 7',
           'trial_8' = 'Trial 8',
           'trial_9' = 'Trial 9',
           'trial_10' = 'Trial 10')

ggplot(tpop, aes(x = time, y = n))+
  geom_point()+
  geom_line()+
  labs(x = 'Time', y = 'Population size')+
  facet_wrap(~trial, nrow = 2, labeller = as_labeller(labels))+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_blank(),
        strip.text = element_text(size = 12))
 



# Remake figure with the mean and min and max values from the multiple trials


  
ggplot(tpop, aes(x = time, y = n))+
  geom_pointrange(stat = "summary",
                  fun.min = 'min',
                  fun.max = 'max',
                  fun = 'mean')+
  labs(x = 'Time', y = 'Population size')+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
 

# Make box plot of change in population size based on type of year 

  
ggplot(tpop, aes(x = year, y = change))+
  geom_boxplot()+
  labs(x = 'Type of year', y = expression(Delta*'Population'))+
  theme_bw()
 

# Remove `NA` and add color to plot

  
ggplot(tpop %>% drop_na, aes(x = year, y = change, fill = year))+
  geom_boxplot()+
  labs(x = 'Type of year', y = expression(Delta*'Population'))+
  theme_bw()
 

# Change order of x axis and color of plot. Colors can be both hex code or from names that R has. A help website for picking colors is [here](https://rstudio-pubs-static.s3.amazonaws.com/3486_79191ad32cf74955b4502b8530aad627.html).


tpop = tpop %>% 
  mutate(year = factor(year, levels = c('good', 'average', 'bad')))%>%
  drop_na

colors = c('good' = 'darkgreen',
           'average' = 'goldenrod',
           'bad' = 'firebrick')

ggplot(tpop, aes(x = year, y = change, fill = year))+
  geom_boxplot()+
  labs(x = 'Type of year', y = expression(Delta*'Population'))+
  scale_fill_manual(values = colors)+
  theme_bw()
 

# Modify the labels and remove the legend

  
ggplot(tpop, aes(x = year, y = change, fill = year))+
  geom_boxplot()+
  labs(x = 'Type of year', y = expression(Delta*'Population'))+
  scale_x_discrete(labels = c('Good', 'Average', 'Bad'))+
  scale_fill_manual(values = colors)+
  theme_bw()+
  theme(axis.title = element_text(size = 18), 
        axis.text.y = element_text(size = 18, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18),
        legend.text = element_text(size = 12))
 

 