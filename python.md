> Teacher: Dr. Augusto César Poot Hernández

"Good morning,

let's open python in the flavor `iphython`

```R
ssh -X ricardo@kneipe.lavis.unam.mx
qlogin
module load python/3.6.5 anaconda3/4.2.0 
source activate env_marisol
ipython
```

some syntax's:

```
install Python
install PIP
# IPYTHON is pretty good
3 + 5 
2 / 19
3 * 8
3**9 # exponencial

```

No need of declaring variables

case sensitive

Letters, digits, underscore

```python
x = 5
y = 4
z = x / y
print(z)
type(z)
float(z)
int(z)

a,b = 3,5 
print(a)
print(b)

```

## Strings

```python
s1 = "Monday"
s2 = "sleepy at UNAM early"
s3 = "I'm Ricardo"
# doueble " " is usefull when are typing a string as "I'm Ricardo"
```

usefull shortcuts

```python
s1 + s2 # print both strings together
s1 * 3 # print three times the string s1
s2.lower() # convert all the string to lowercases
s2.upper() # convert all the string to UPERCASES

s2.count("a") # count the times than 'a' are in the string
s2.replace("e", "k") # replace letter 'e' to 'k' from the string
s2.split() # split / cut the string
s2.split("n") # use a special character to separete the string
"y" in s2 # logical assigment. does 'y' in the string s2
len(s1)
```

Good! let's continue:

## Variables

- No need of declaring variables

- Case sensitive

- Letters, digits, underscore

- Don't start with digits

- Don't use keywords

```python
whos # print the number of variables created
Variable   Type     Data/Info
-----------------------------
a          int      3
b          int      5
s1         str      Monday 

```

using variables

```python
L1 = [2,3,'hello', 3.4]
L2 = [2,3,4]
L1[1] # Python start from 0 to n 
Out[70]: 3
```

Using complex objects calling

```python
L3 = [1,2, [10,20], [1, [5,6]]]
L3[2]
# Out[80]: [10, 20]
L3[2][1]
# Out[81]: 20
L1[-1] # call the last element of the vector
# Out[86]: 3.4
L1[2:5]
# Out[135]: ['HOLA', 'hello', 3.4]
L1[2:] # or print everything after vector 2
L1[:5] # or print everything before vector 5
```

Add, and append values into a vector

```python
L1.append(2500)
#
L1 + [3000]
#Out[95]: [2, 3, 'hello', 3.4, 2500, 3000]

L1.insert(2, "HOLA") # Insert an element in specific position in the vector
L1
# Out[97]: [2, 3, 'HOLA', 'hello', 3.4, 2500]
```

Good. Try this examples

| j    | u    | r    | i    | q    | u    | i    | l    | l    | a    |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9    |
| -10  | -9   | -8   | -7   | -6   | -5   | -4   | -3   | -2   | -1   |

```python
s = 'juriquilla'
s[1:4]
s
```

# Data types

## Dictionaries

> It type of objects are pretty usefull to manipulate big-data-sets

Unlike sequences, which are indexed by a range of numbers, dictionaries are indexed by *keys*, which can be any immutable type; strings and numbers can always be keys. Tuples can be used as keys if they contain only strings, numbers, or tuples; if a tuple contains any mutable object either directly or indirectly, it cannot be used as a key. Details [here](https://docs.python.org/2/tutorial/datastructures.html#dictionaries)

```python
biothing = {'blood':'red', 'legs':2, 'bionumber':[1,2,3] }

biothing['legs']
# 2
biothing['blood']
# red
biothing[(2,3)] =  3.1416
biothing
# {'blood': 'red', 'legs': 2, 'bionumber': [1, 2, 3], (2, 3): 3.1416}

```

and lets figure out _keys_ and _values_ with in the object `biothing`

```python
biothing.keys()
# Out[150]: dict_keys(['blood', 'legs', 'bionumber', (2, 3)])
biothing.values
# Out[151]: dict_values(['red', 2, [1, 2, 3], 3.1416])
 
biothing['inception'] = {1: 'a', 2: 'b'}
biothing['inception'][1]
# Out[154]: 'a'
```

Using the `dict` function

```python
anotherDict = [('mexico', 1), ('germany', 0), ('spain', 3)]
dict(anotherDict)
# Out[157]: {'mexico': 1, 'germany': 0, 'spain': 3}
```

Next example:

```python
A = [1,2,3,4,5]
B = [10,20,30,40,50,60]
list(zip(A,B)) # This is the rbinb or cbind version in python
```

## Sets

The [`sets`](https://docs.python.org/2/library/sets.html#module-sets) module provides classes for constructing and manipulating unordered collections of unique elements. Common uses include membership testing, removing duplicates from a sequence, and computing standard math operations on sets such as intersection, union, difference, and symmetric difference.

```python]
A = {2,3,4,5}
B = [2,3,4,5,5,6]
sets(B)
C = sets(B)

A = {2,2,2,3,4,4,5,7,7}

A.intersection(C)
A.union(C)
A.difference(C)

```

Good!, let's make new stuffs

```python
A = input('Your name ')
# Your name type_yourname
print(A)
#Out[9]: 'Ricardo'

input('Your age ')
#Your age type_your_age
print(B)
# Out[10]: '26'

print('hello ',  A, 'good to know your age', B)
# hello  Ricardo good to know your age 26


```

### Exercises

```python
print('Hello world')
A = input('Let me know your name')
B = input('... and your age, please ?')
print('Good to knwo you ', A, "I'm also ", B, " Years old bro!")

print(" Ok... bye, bye world and", A)
```

> Hardwork 2: A script than ask for two numbers and prints the sum of them

```python
print("Hello mathematica")
print("Let's sume some values")
A = input("Type your first and second number")
B = input("... second number")
C = int(A) + int(B)
print("Well, the sum of both values is: " + str(C))
```

> Hardwork 3:  Script than ask for your lunch's bill and print the final bill with 10 % of tips.

```python
print("Good, last question")
L = input(" What was you lunch bill?")
l = float(L)
p = l *  0.1
print(" Ok, you finall bill 'd be ", str(l+p))
```

### Boolean operators

**Spaces** are really important !!

- They define blocks of code
- Beware of mixing tabs and spaces
- Guideline: 4 spaces each level



```python
<, <=, >, >=. ==, !=, and, or
```

While

```
b = 10
a = 
```

## Functions

A difference from the `methods` functions have a general purpose to apply to any object. 

Python has libraries and modules with defined functions

```python
import math
math.sqrt(3)
# Out[37]: 1.7320508075688772
```

Also, you can rename the library name (where are the modules you want to use )

```python
import math as M

In [39]: M.sqrt(3)
Out[39]: 1.7320508075688772
```

Or import sepecific module from a package

> Note: this is a **dangerous** way due to code. Please try safety mode using the complete library name from module came from.

```python
from math import sqrt
sqrt(5)
```

```python
L2 =  ['Hello', 'How are you', "bye bye", "ok"]
random.choice(L2)
random.sample(L1, 5)
ramdom.sample(L2, 5)
```



## # Good examples

```python
In [4]: def sphereVolume(r):
   ...:     vol = (4/3.0) * math.pi * (r**3)
   ...:     return vol
   ...: 
   ...: 

# Then,   
sphereVolume(2)
# Out[5]: 33.510321638291124
```

> example 2

A module can be have many functions. Let's make a library of name `myFirstFunction.py ` with two functions as follow:

```python
def addTip(bill):
    return 1.1*bill

def equallyDivide(bill, n):
    return bill / n
```

Then in `ipython` use it:

```python
import myFirstFunction
myFirstFunction.addTip(100)
# Out[17]: 110.00000000000001
myFirstFunction.equallyDivide(100, 4)
# Out[18]: 25.0
```

Love it!

> Second example

If you're changing any into your functions, can use `%load_ext autoreload ` in order to upload the newer version of your library. (lets supose we change the `n` from our `myFirstFunction.py ` script); then

```python
%load_ext autoreload
%autoreload 2
myFirstFunction.equallyDivide(100)
```

## Reading files

Methods to read files are:

`read`: single **string** with the **whole** text

`readlines`: Single **string** with single lines

`readlines()`: **List of strings**, each string one line of the file

```python
f = open("/home/ricardo/python/files/totalEclipse.txttalEclipse.txt", "r")
o1= f.read() # read all elements 
f.seek(0) # back to the fist line from the object f
o2 = f.readlines() # read all elements as a list

f.close()

```

| Options    |                            |
| ---------- | -------------------------- |
| `r`        | reading                    |
| `r+`       | reading and writing        |
| `w`        | writing                    |
| `a`        | Append list                |
| `rb`, `wb` | Binary reading and writing |

Create files is similar to read it:

```python
h = open("/home/ricardo/python/files/somethingNew.txt", "w")
h.write("hello bioinformatics people")
h.close()
```

And append lines

```python
f1 = open("/home/ricardo/python/files/somethingNew.txt", "a")
cheesyLines = [" turn around", "ever now and then", " i get a little bit"]
f1.writelines(cheesyLines)
f1.close()
```

#### Exercises

1. Using python, create a new file **SUPERtotalEclipse.txt** with the same lines that **totalEclipse.txt**, but written in uppercase letters

> Bonus: replace the word `heart` with the word `moon`

```python
def UPERCONVERT(filename):
    i = open(filename, "r")
    content = i.read()
    o = open("SUPERtotalEclipse.txt", "w") # create out file
    # write the file by pipe the replace and upper function through your data content
    o.write( content.replace("heart", "moon").upper() )   
    o.close() # close it :)
    i.close()     
    return print("bye")
```

Then, import your functions and run into your data:

```python
import exercise1
exercise1.UPERCONVERT("/home/ricardo/python/files/totalEclipse.txt")
bye

```



2. Using python, create a function that concatenates two given files with a specific name (your function, thus , requires three strings as input). Then use your function to create the file **kareoke.txt** with the contents of the **SUPERtotalEclipse.txt** and **lifenMars.txt**

```python
def concatenate(file1, file2):
    f1 = open(file1, "r" )
    f2 = open(file2, "r")
    o = open("kareoke.txt", "w")
    content = f1.read() + f2.read()
    o.write(content)
    o.close
    return print("good")
```

Then, run

```python
import exercise2
exercise2.concatenate("/home/ricardo/python/files/totalEclipse.txt", "/home/ricardo/python/files/lifeInMars.txt")

```

## Comprehension lists

```python
nameList = ["octavio", "ali", "marcos", "pedro", "juan", "lucas", "pablo"] #a
nameL = [len(word) for word in nameList] # b
nameL2 = [word for word in nameList if len(word) < 5] # c
[len(word) for word in nameList if len(word) < 5] # d
```

a. Make a vector of names

b. Count the length of each element in the vector and save in new object

c.  Also apply logical filter for the string lenght (print elements)

d. or print the lenght of the elements resulted from the logical filter

```python

```





### Final words

#### Sorting files

for intergers of boolean

```python
a = [5,1,4,3]
b = [1.5, 4.3, 3.9, 9.6]
sorted(a)
# [1, 3, 4, 5]
sorted(b)
# [1.5, 3.9, 4.3, 9.6]

sorted(a, reverse=True)
# [5, 4, 3, 1]
```

Also, for intergers or boolean

```python
sorted(nameList)
sorted(nameList, key = len)
```



Good, lets do some def

in order to print the help from modules and commands try `help()`

```python
help(sorted)
#
import math
help(math)

```

```python
def lastLetter(word):
    return word[-1]
```



#### Visualizing

Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms.

> We want ot plot x**2

```python
x = [-1, 0, 1, 2]
y = [1, 0, 2, 4]
```

then, 

```python
import matplotlib.pyplot  as plt
plt.plot(x,y)
plt.show()
```

- matplotlib
- Pandas [ref](https://pandas.pydata.org/) 