> By Michael Jeziorski
>
> starting Rstudio session [here](http://192.168.100.45:8787) 

- R is a oriented objects lenguaje.
- everything in R is an object.
- An object is a data structure having some attributes and methods which act on its attributes.

###### Cheats 1

Let's check the **installed packages** where you're working on:

```R
ip <- installed.packages()[,c(1,3:4)]
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
print(ip, row.names=FALSE)
```

Although we select package, package's version and Priority of the Package, there're others colnames of interest in the function installed.package:

```R
 [1] "Package"              
 [2] "LibPath"              
 [3] "Version"              
 [4] "Priority"             
 [5] "Depends"              
 [6] "Imports"              
 [7] "LinkingTo"            
 [8] "Suggests"             
 [9] "Enhances"             
[10] "License"              
[11] "License_is_FOSS"      
[12] "License_restricts_use"
[13] "OS_type"              
[14] "MD5sum"               
[15] "NeedsCompilation"     
[16] "Built"
```



###### Cheats 2

Some package have error during load/use because version of dependence need to `unload` some package  than cause de dependency error. Then:

```R
if("package:vegan" %in% search()) detach("package:vegan", unload=TRUE) 
```

##### Cheats3

use exclamative expression `!` in a logical function (as `na.rm`)  in order to say not. Example:

```R
if (!require(x,character.only = TRUE, quietly=TRUE))
      install.packages(x,dep=TRUE, repos='http://cran.us.r-project.org')
```

# Control structure 

Computation in R consists of sequentially evaluating *statements*. Statements, such as `x<-1:10` or `mean(y)`, can be separated by either a semi-colon or a new line. 

> Control structure in R

Both, **semicolon** and new lines can be used to separate statements. A semicolon always **indicates the end of a statement** while a new line may indicate the end of a statement. If the current statement is not syntactically complete new lines are simply ignored by the evaluator. If the session is interactive the prompt changes from `>` to `+`.

```R
x <- 0; x + 5
[1] 5
y <- 1:10
1; 2
[1] 1
[1] 2
```

Statements can be grouped together using braces `{’ and ‘}`. A group of statements is sometimes called a *block*. Single statements are evaluated when a new line is typed at the end of the syntactically complete statement. Blocks are not evaluated until a new line is entered after the closing brace. In the remainder of this section, *statement* refers to either a single statement or a block.

```
> { x <- 0
+ x + 5
+ }
[1] 5
```

Braces are usefull in condition with `if`,  `while`,  looping with `for`,  `repeat` and `switch`.

>  `if`/`else` statements can be nested.

```
if ( statement1 ) {
    statement2
} else if ( statement3 ) {
    statement4
} else if ( statement5 ) {
    statement6
} else
    statement8
```

Reference: 
https://cran.r-project.org/doc/manuals/r-devel/R-lang.html



## Single object

> Set paths

```R
getwd()
setwd()

```

> Operation and variables in r

```R
3 + 4
3 * 4
3 / 4
x <- 3 ^4
```

> install, load and unload package

```R
install.packages("tidyverse")
library("tidyverse")
detach("package:tidyverse", unload=TRUE) 
```



### Type of objects

- Numeric or double (boolean numbers)
- Integer (Numeric intergers)
- Character (any ASCII character)
- Factor (3 levels of Control, experiment 1 and experiment 2)
- Logical (TRUE or FALSE)
- Complex

```R
class(x)
[1] "numeric"
```

Let's in to create vector, matrix, lists, data.frames which

> vectors

```R
y <- 1:10
z <- seq(1,10, by= 0.5)
z1 <- seq(11,20, by= 0.5)
```

> And manipulate it with stats

```R
y * 3
sum(y)
mean(y)
z+z1
```

Difference between data.frame and matrix?

> Data.frame is a object with more than one dimensions than save vectors from different classes
>
> **Matrix** is like a data.frame but, **every data-vector as the same class** within this object  (ie. Numeric, character, interger etc.).
>
> A **matrix** could be convert to data.frame but not necessary inversely

also can represententing values in a a variable and do more stuffs:

```R
ex1 <- c(2,4,8,12)
index <- 3

identical(ex1[3],
          ex1[index])
```

Also manipulate/combine more than one row from this vector using `c`

```R
ex1[c(2,4)] # select vector 2 and 4
ex1[c(1:3, 1)] # select from 1 to 3 and then back 1
ex1[-3] # remove vector 3
ex1[-c(3,4)] # remove more than one vector
```

```R
c(x,z) # combine two set of vectors
days <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
```



And manipulate logical condition

```R
ex1 > 5
ex1[ex1 > 5] # filter the vector by the first condition
```

other logical conditions:

`>` greater than

`<` minor than

`<=` minor or equal than

`>=` greater or equal than

`!=` Not equals to

`==` Equal to

`|` OR

`&` AND

## Two-dimensional objects

Let's use `data()` in order to list all two-dimensional data sets available in R. Then lets use `mtcars`. It is a data.frame (a lists of vector of the same length).

```R
str(mtcars)
```

> using str (structure) print  a summary of any object

Let's work in some examples of work with two-dimensional objects

```R
mtcars[2, 3] # select row 2 in the column 3
mtcars["Valiant", "hp"]  # select row of name "Valiant" in column "hp"
mtcars["Valiant",] # select all rows of name "Valiant"
mtcars[, "hp"] # select all columns of name "hp"

mtcars[, 1:4] # select all rows in columns 1 to 4
mtcars[, c(1, 11)] # select all raws in columns 1 and 11

mtcars[rownames(mtcars) =="Valiant",] # select all columns than contain rowname "Valiant"

mtcars$hp # due to it is atomic vector lets use $ to call the column hp
mtcars[, "hp"] # the same than above

mtcars[mtcars$hp <= 99, ] # in column hp condition to print all values less than 99, then back all the columns than adjust to this condition.
```

Lets modify rows in the object in order to modify the velocity of the fiat 128 from 66 (`mtcars["Fiat 128", "hp"]`) to 366 (`cars["Fiat 128", "hp"]`)

```R
cars <- mtcars
cars["Fiat 128", "hp"] <- 366
identical(cars["Fiat 128", "hp"], mtcars["Fiat 128", "hp"])
```

Figure out the faster car by `hp` column [ordering decreasing](https://www.statmethods.net/management/sorting.html). We're going to watch Fiat 128 in the top, following by Maserati Bora. 

```R
cars <- cars[order(cars$hp, decreasing = T), ] ; head(cars)
```

Use a **single vector** option in order to modify the row names. `rownames(cars)` is a single vector. lets show the first vector from this object.

```R
rownames(cars)[1]
rownames(cars)[1] <- "Mazda" # rename only this rowname
```

A complex version of the first:

```R
rownames(cars)[c(1,2)] <- c("mazda1", "mazda2")
```

## Condition examples

> Quiz: Lets figure out a condition than:
>
> Select values not equal to 4 in column cyl (ie. remove rows as value 4) or values are less than 30 in column mpg. Print all the columns  

```R
mtcars[mtcars$cyl != 4 | mtcars$mpg < 30, ]
```

> Good!

```R
colMeans(cars)
```

## Lists

> subsetting a list created as [dommies explain](http://www.cookbook-r.com/Numbers/Generating_random_numbers/) 

```R
ls <-  list(names=c("JK","RJ","EG","MZ","MS"),
    		age = sample(seq(20,30, by=1), 5),
            weight = runif(5, min=60, max=80),
            eight = runif(5, min=1, max=2))
```

input in the list object

```R
ls[[2]][1] # from column 2, select row 1
ls$age[1] # is the same than above !!
```

###  Approach to wragling data (Tidyverse)

```R
library(tidyverse)
```

`filter` = Choose rows based on condition

`select` = Choose columns

`arrange` = Sort rows based on column values

`mutate` = Convert data in one column into new data in a new column

`summarize` = Provide summary data for a column

`group by` = Group data based on a variable; often used with summarize

>  Examples 1:

a) by default tidyverse object remove colnames, lets convert it first

b) condition by column cyl and mpg and pipe to the next condition

c) arrange desceding the data.frame object descending by mpg

d) finally select first two columns by

```R
newcars <- rownames_to_column(mtcars) %>% # a
	filter( cyl == 4 & mpg < 30) %>% # b
	arrange(desc(mpg)) %>% # c
	select( -c(3:11)) # d

```

> Example 2

a) From mtcars object sort column `cyl`  then summarize the frequency of numbers by its mean

```R
mtcars %>% 
		group_by(cyl) %>% 
		summarize(mean = mean(mpg)) #a
```

# Wrangling real biological data



Using blat aligment dataset lets wrawling in R.

```R
rm(list = ls())
```

> When load data, remember:  not convert id names (usually rownames) as factor, ie. Turn on the option `stringAsFactors = FALSE` due to it are unique identifier instead of levels of factors in the object.

Load data

```R
blat <- read_tsv("file.psl", stringAsFactors = FALSE )
#or
blat <- read.psl("file.psl")
```

And start summary results as:

a) grouping the query aligments and summarise by count `n()`, then arrange descending

```R
blat %>%
	group_by(chr) %>%
	summarise(n()) %>%
	arrange(desc(`n()`)) # a
```

Determine number of unique genes per chromosome

```R
for i in seq()
```

...

# Infovis your data 

```R
library(ggplot2)
```

Using the mtcars dataset lets plot some stuffs.

> First load the plot in a object

```R
hp_vs_mpg <- ggplot(mtcars, aes(hp, mpg, color = as.factor(cyl))) + geom_point()
```

> Then, smooth the plot by lm Coefficients from `lm(mtcars$mpg ~ mtcars$hp)`

```R
hp_vs_mpg +
	geom_smooth( method = "lm", color = "orange", se = FALSE)
```

Good!

> Now, we can try wrap a different set of levels 

```R
hp_vs_mpg +
	facet_grid( .~ gear)
```



## The night

lets test some condition using `if` using the mtcars dataset

```R
a <- 5
if (a > 10) {
    	print("value is over 10")
} else {
    	print("value is equal to or under 10")
}
```

let's improve with a loop

```R
b <- sample(seq(1,100, by=1), 30)

for (i in 1:length(b)) {
    if (b[i] >10 ) {
        print ("value is over 10")
    } else {
        print("value is equal to or under 10")
    }
}

```

Also, using the family of apply functions ([paper here](https://www.datacamp.com/community/tutorials/r-tutorial-apply-family)) we can sustitute the loop

```R
lapply(b, fun)
# where fun:
fun <- function(x) {
    if (x > 10 ) {
     print ("value is over 10")
 } else {
     print("value is equal to or under 10")
 }
                   }
```

