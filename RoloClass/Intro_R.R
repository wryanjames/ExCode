
#R is useful for basic operations and follows math rules (i.e. PEMDAS). R will all code on a line unless there is a `#` to the left.

## Getting to know the basics 

#R is a programming language that has become the standard in Ecology due to its flexibility and open source nature. R can be used from simple math to complex models and is very useful for generating figures. R, like all computer languages, using a specific syntax to run commands that it is programmed to do. In other words, R will only do what it is commanded to do, and therefore, many common errors are due to errors in syntax (e.g. misspellings, missed commas, or unclosed brackets). 

# This example gives a basic intro into R syntax that can be useful for ecological modeling. This script gives examples of how to:
#   
# 1.  Basic operations in R
# 2.  Assigning objects
# 3.  Types of data structures in R
# 4.  Functions in R
# 5.  Using Packages in R
# + How to install and load packages
# + tidyverse
# + tidy data
# + piping 
# 6.  Indexing
# 7.  Vector operations 
# + Custom Functions
# + `plyr` and `dplyr`
# + `purr`
# 8.  For loops
# + Plotting 


## Basic operations in R

# R is useful for basic operations and follows math rules (i.e. PEMDAS). R will all code on a line unless there is a `#` to the left.


# addition 
1+1 

1+1 # + 2 (won't run anything to right of #)

# subtraction
5-2 

# multiplication
4*5

# division
33/5

# exponents can be done 2 ways
2^2
2**2

# follows PEMDAS
1+5*4
# different answer than above
(1+5)*4


# Note the `[1]` appears next to your result. R is just letting you know that this line begins with the first value in your result. Some commands return more than one value, and their results may fill up multiple lines. For example, the command 100:130 returns 31 values; it creates a sequence of integers from 100 to 130. Notice that new bracketed numbers appear at the start of the first and second lines of output. These numbers just mean that the second line begins with that value. You can mostly ignore the numbers that appear in brackets:
  

100:130

## Assigning objects
# When working in R it is useful to store data as an object. Assigning objects can be done in multiple ways, but the most common are `<-` and `=`. These objects are stored in the R environment and can be called. Objects can be assigned multiple times, but only the last assignment is what is stored. Also it is important to know that R is case sensative and capital and lower case numbers are different.


# assign an object
a = 4 
a

b <- 23

a+3 

b/2

a*b

c = 8
c = 14
c

d = 15 
D = 1 
d
D


## Types of data structures in R
# R has 6 basic data types. (In addition to the five listed below, there is also raw which will not be discussed in this workshop.)
# 
# + integer
# + numeric (real or decimal)
# + character
# + logical
# + complex
# 
# integers are whole numbers
# 
# numeric are numbers with decimals. Integers and numeric are different because of how the underlying data is stored. 
# 
# characters are strings of letters and numbers (e.g. `"abc"` and `"b1x"`) and are designated in R by `" "`. When using characters, `" "` are required because in R letters without quotations are objects and `c = 'd'` is different than `c = d`
# 
# logical is `TRUE` or `FALSE`. One thing to note is that `T` is the same as `TRUE` and `F` is the same as `FALSE`. Because `T` and `F` are special in R they cannot be used to name objects (but `t` and `f` are ok because R is case sensative). This is true for other cases as well like `NA` and `NULL`.  
# 
# complex numbers have both real and imaginary parts (`1+4i`)
# 
# Elements of these data types may be combined to form data structures, such as atomic vectors. When we call a vector atomic, we mean that the vector only holds data of a single data type. A vector is the most common and basic data structure in R and is pretty much the workhorse of R. Technically, vectors can be one of two types:
#   + atomic vectors
# + lists
# although the term "vector" most commonly refers to the atomic types not to lists.
# 
# There are different ways to make vectors



# make a numeric vector
a = c(1.1,5,3,4)
a

# make a integer vector
b = 1:15
b

# make a character vector 
c = c('a', 'b', 'c')
c 




# Because characters can be both letters and numbers, numbers in a vector with letters are stored as a character. These cannot be used for math operations, but integers and numeric data types can be used for math. 


a = 4.4
a / 1 


b = 6L # L can be used to keep a numeric as an integer, R typically defaults to numeric
b*3

# character
c = '1'
c*4


# Another common way to store data is in a dataframe or tibble (special type of dataframe from the `tidyverse` package we will see below). This is a collection of atomic vectors with the same length. 


b = data.frame(c1 = c(1,2,3), c2 = c('a','b','c'))
b



## Functions in R
# R comes with functions that are used to do tasks. Functions take arguments to complete a task. Functions have the general format `function(argument1 = , argument2,...)` The types of data used and output of the function is specific to that function. Below are just a few useful examples. 


# summary statistics of sequence of numbers
a = c(1.1,5,3,4)
mean(a) #mean
median(a) #median
sd(a) #standard deviation
quantile(a, 0.5) # quantile at 0.5 (median)

# make a sequence of numbers
b = 1:15
b
c = seq(1,15,1) #more flexibility than :
c
seq(4,20,2)

# information about objects
d = c('a', 'b', 'c')
typeof(d) 
typeof(c)
length(d)

# dataframe/tibble specific functions
e = data.frame(c1 = c(1,2,3), c2 = c('a','b','c'))
names(e) # column names
nrow(e) # number of rows
length(e) # for dataframe number of columns
str(e)# structure of data


## Using Packages in R
# R comes with a lot of base functions that are available for use when you open R, but this does not contain all of the functions useful to your tasks in R. Since R is open source, many R users have created Packages that contain functions that can be downloaded. Which includes the very common `tidyverse`.

### How to install and load packages
# Packages can be downloaded from CRAN or from Github. To download directly from Github other packages are needed. 

# install.packages('tidyverse') #from cran


# Once downloaded, packages can be loaded into the R environment with `library()` function. Packages have to be loaded each R session. In addition functions can be called directly from a package with `::` in the format of `packageName::function()`. 

library(tidyverse)


### tidyverse
# [tidyverse](https://www.tidyverse.org/) is a collection of packages that use similar syntax and are used for data science in R. Coding in tidyverse is typically easy to read and understand, and has useful functions that have been adopted into newer versions of base R (e.g. piping). Tibbles are the tidyverse version of a dataframe.

c = tibble(c1 = c(1,2,3), c2 = c('a','b','c'))
c

#### tidy data
# Data is collected and stored in many different ways, which can make it difficult to analyze. One of the goals of tidyverse is to easily turn messy data into tidy data which can easily be analyzed. In tidy data:
#   
#   1. Every column is a variable.
# 2. Every row is an observation.
# 3. Every cell is a single value.
# 
# Two functions `pivot_longer()` and `pivot_wider()` are useful in manipulating data stored in rows and columns. ***Note that `pivot_longer()` and `pivot_wider()` have replaced `gather()` and `spread()` in newer versions of `tidyverse`


#tidying data 
stock = tibble(name = c('GOOG', 'AMC', 'GME'),
               Jan = c(1000, 2, 4),
               Feb = c(1010, 15, 30),
               March = c(1005, 25, 180))

df = stock %>% 
  pivot_longer(cols = Jan:March, 
               names_to = 'Month',
               values_to = 'Price')

df

# wide format
fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))
fish 

fish %>%
  pivot_wider(id_cols = species, 
              names_from = year,
              values_from = catch)


#### piping  
# Tidyverse has an operator `%>%` known as a pipe that is useful for when you want to do multiple actions to the same data. It takes the output of the left of the `%>%` and makes it the first argument of what is on the right. Allowing to reduce code and make things tidier.


# this code
df = as_tibble(mtcars)
df = filter(df, mpg > 20)
df = mutate(df, color = 'red')
df = select(df, mpg, cyl, color)

head(df)

# can become

df = mtcars %>%
  as_tibble()%>%
  filter(mpg > 20)%>%
  mutate(color = 'red')%>%
  select(mpg, cyl, color)

head(df)


##  Indexing
# Once data is stored in an object, being able to retrieve those values is useful. Referred to as indexing, the syntax is specific to how the data is stored. With indexing specific values within your object can be modified. 


# vector 
b = 1:15
# 3rd object 
b[3]

# make a character vector 
c = c('a', 'b', 'c')
c
# 2nd object
c[2]
# change 
c[2] = 'new'
c

# dataframe and tibbles
mtcars
# first column
mtcars[1]
# first row
mtcars[1,]
# 2nd row of first column
mtcars[2,1]
# can call specific columns (called as a vector)
mtcars$mpg
mtcars$cyl
#same for tibble
d = mtcars %>% as_tibble
d[1]
d$mpg
d$cyl
# specific row in specific column
mtcars$cyl[1]
d$cyl[1]



##  Vector operations 
# As we have seen above, we can do operations over vectors. We sometimes want to do this to vectors stored in dataframes/tibbles, and the `mutate()` function makes this easy. 


iris %>% 
  mutate(petalArea = Petal.Length*Petal.Width)

iris %>%
  mutate(petalArea = Petal.Length*Petal.Width,
         PetalSize = if_else(condition = petalArea > 0.2, true ='big',
                             false = 'small'))

iris %>%
  mutate(petalArea = Petal.Length*Petal.Width,
         PetalSize = if_else(condition = petalArea > 0.2, true ='big',
                             false = 'small'))%>%
  group_by(PetalSize)%>%
  summarize(mean = mean(Petal.Width),
            n())


### Custom functions
# So far we have used functions that are apart of base R or from different packages, but we can also build custom functions. Custom functions are useful when existing functions will not do the task at hand, or when combining functions over multiple times. These can be used for vector operations.


# custom funciton 
my_fx = function(x){
  b = x * 23
  return(b)
}

my_fx(5)

t = seq(2,10,2)
my_fx(t)

# custom function with multiple arguments
my_fx2 = function(x,y){
  b = x/y
  return(b)
}

my_fx2(4,5)

# custom fucntion with true false

my_fx3 = function(x = T){
  if (x == T){
    cat('x is true \n')
  }else{
    cat('x is not true \n')
  }
}

my_fx3(x = T)
my_fx3()
my_fx3(x = F)
my_fx3(x = 5)


# custom functions with vector operations
mi_km = function(mi){
  km = mi * 1.60934
  return(km)
}

mtcars %>%
  mutate(kmpg = mi_km(mpg))




# summarizing data
summ <- d %>% 
  group_by(., year) %>% 
  summarise(mean.count = mean(count),
            sd.count = sd(count),
            cv = sd.count/mean.count) %>% 
  data.frame()
summ


#dplyr Example 2
summ.3 = d %>% 
  group_by(., year) %>% 
  mutate(std.value = (count - mean(count))/sd(count),
         cum.sum = cumsum(std.value))%>% 
  data.frame()

summ.3



# ##  For loops
# Another useful tool in programming is `for` loops. For loops repeat a process for a certain number of iterations. These can be useful iterate over a dataset or when using information in a time series. The `for` loop works over the number sequence indicated and does the code within the loop (inside of `{}`) for each number in the sequence. The iteration is typically indicated with `i`, but is just an object that is replaced at the begining of each loop and can be anything.


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

d = seq(1,15, 2)
d
for(i in 1:length(d)){
  b = d[i] + 1
  cat('d =',d[i], 'b = d + 1 =', b, '\n' )
}

b = 1:10
for (i in 2:10){
  z = b[i] - b[i-1]
  
  cat('z =', z, 'b[i] =', b[i], 'b[i-1] =', b[i-1], '\n')
}


start = 10 
pop = tibble(time = 0:10, n = NA)
pop$n[pop$time == 0] = start
pop
for (t in 1:10){
  growth = rnorm(n =1, mean = 3, sd = 1)
  pop$n[pop$time == t] = growth + pop$n[pop$time == (t-1)]
}
pop


start = 1000 
pop = tibble(time = 0:10, n = NA)
pop$n[pop$time == 0] = start
pop
for (t in 1:10){
  growth = rnorm(n =1, mean = 300, sd = 50)
  death = rnorm(n = 1, mean = 300, sd = 50)
  pop$n[pop$time == t] = growth - death + pop$n[pop$time == (t-1)]
}
pop



#### Plotting 
# The `ggplot2` package is part of the packages that load with `tidyverse` and has become the standard in ecology. The syntax builds upon on a base function and is very customizable [see cheat sheet](https://www.rstudio.com/resources/cheatsheets/). 

start = 1000 
pop = tibble(time = 0:100, n = NA)
pop$n[pop$time == 0] = start

for (t in 1:100){
  growth = rnorm(n =1, mean = 300, sd = 50)
  death = rnorm(n = 1, mean = 300, sd = 50)
  pop$n[pop$time == t] = growth - death + pop$n[pop$time == (t-1)]
}

ggplot(pop, aes(x = time, y = n))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = 'Time', y = 'Size of population (n)')



for (i in 1:5){
  trial = paste0('trial',i)
  start = 1000 
  pop = tibble(time = 0:100, n = NA, trial = trial)
  pop$n[pop$time == 0] = start
  
  for (t in 1:100){
    growth = rnorm(n =1, mean = 300, sd = 50)
    death = rnorm(n = 1, mean = 300, sd = 50)
    pop$n[pop$time == t] = growth - death + pop$n[pop$time == (t-1)]
  }
  
  if (i ==1){
    tpop = pop
  }else{
    tpop = bind_rows(tpop,pop)
  }
}

ggplot(tpop, aes(x = time, y = n, color = trial))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = 'Time', y = 'Size of population (n)')+
  theme(legend.position = 'bottom')

