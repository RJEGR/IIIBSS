`igraph` is a hug R library. Igraph is available as R, python and C libs. It can be installed in R by `install.packages(“igraph”)`

```R
library(igraph)
```
An empty graph/network with five spherical yellow vertex (edges) is created by:

```R
g <- make_empty_graph(n=5, directed=TRUE)
#V(g) means Vertex(g)
V(g)$color = "yellow"
V(g)$shape = "sphere"
```
```R
plot(g)
```

Add the next edges to the nework:
```
g <- add.edges(g, c(1,2, 1,3, 2,4, 3,4, 4,5))
plot(g)
```

Add a red spherical **vertices** (vertex) and add up to the network the next edges : 3->6, 6->5.

```R
g <- add.vertices(g, 1, color="red", shape="sphere")
plot(g)
```

Then add, **edges**

```R
g <- add.edges(g, c(3,6, 6,5))
plot(g)
```

The class (vector, matrix, data.frame, network, etc) of an object in R is obtained using the command `class(object)`. _Igraph_ provides its own class of object.

```R
class(g)
## [1] "igraph"
str(g)
```

Replace connection 1 to 3 with connection 3 to 1.

```R
g <- delete.edges(g, c(2))
g <- add.edges(g, c(3,1))
```

Name the nodes consecutively with the letters A-F.

```R
V(g)$name <- LETTERS[1:6]
V(g)
## + 6/6 vertices, named, from 8732b2a:
## [1] A B C D E F
E(g)
## + 7/7 edges from 8732b2a (vertex names):
## [1] A->B B->D C->D D->E C->F F->E C->A
```

## Degree

The command for the degree is the function **degree(igraph_object)**

```R
degree(g)
## A B C D E F 
## 2 2 3 3 2 2
```

Plot the network in such way that that size of the nodes is propotional to the number of input connections

```R
plot(g, layout=layout_nicely, vertex.size=degree(g, V(g), "in")*15+15,
     vertex.label.dist=0.5, edge.arrow.size=0.5)
```

Find and plot the degree distribution

```R
plot(degree_distribution(g), main="Degree distribution", xlab="Degree", ylab="Frequency")
hist(degree(g),col="salmon")
```

 The adjacency matrix is obtained by **get.adjacency(igraph_object)**

```R
adj_mat<-as.matrix(get.adjacency(g))
adj_mat
```

Plot a heat map from this matrix

```R
heatmap(adj_mat, Rowv=NA, Colv="Rowv")
```

Some special networks are available in igraph

- star (`make_star()`)
- ring (`make_ring()`)
- tree (`make_tree()`)
- lattice (`make_lattice()`)
- kautz (`make_kautz_graph()`)

Try to figure out what kind networks generate the previous commands before running them

# Real-data example

```R
complete_address <- "/home/ricardo/Rstudio/netwrks/"
nodes <- read.csv(paste0(complete_address,"Dataset1-Media-Example-NODES.csv"), header=T, as.is=T)
links <- read.csv(paste0(complete_address,"Dataset1-Media-Example-EDGES.csv"), header=T, as.is=T)
```

Let's examine data

```R
head(nodes)
```



|      |                     |            |            |               |
| ---- | ------------------- | ---------- | ---------- | ------------- |
| id   | media               | media.type | type.label | audience.size |
| s01  | NY Times            | 1          | Newspaper  | 20            |
| s02  | Washington Post     | 1          | Newspaper  | 25            |
| s03  | Wall Street Journal | 1          | Newspaper  | 30            |
| s04  | USA Today           | 1          | Newspaper  | 32            |
| s05  | LA Times            | 1          | Newspaper  | 20            |
| s06  | New York Post       | 1          | Newspaper  |               |

 ```R
head(links)
 ```

| rom  | to   | type      | weight |
| ---- | ---- | --------- | ------ |
| s01  | s02  | hyperlink | 22     |
| s01  | s03  | hyperlink | 22     |
| s01  | s04  | hyperlink | 21     |
| s01  | s15  | mention   | 20     |
| s02  | s01  | hyperlink | 23     |
| s02  | s03  | hyperlink | 21     |

## Creating an igraph object

....

## Concepts:

> Eccentricity
>
> 



Calculate the degree, average degree and degree distribution

```

# of the following networks
# Degree

k <- seq(1,10)
k <- matrix(k, 10,2)
N <-length(k) /2 

L <- apply(k, 1, function(x) sum(x)*0.5)
Lk = apply(k, 1, function(x) sum(x) / N )

Pk <- apply(k, 1, function(x) k / N)

head(Pk, n = 10)

plot(head(Pk, n = 10))
```

