# Class 6: R functions
Joseph Lo (PID: A18121493)

- [Background](#background)
- [Our first function](#our-first-function)
- [A second function](#a-second-function)
- [A Protein generating function](#a-protein-generating-function)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function.
- One or more input **arguments**
- The **body** the lines of R code that do the work

## Our first function

Let’s write a silly wee function called `add()` to add some numbers (the
input arguments)

``` r
add <- function(x, y) {
  x + y
}
```

Now we can use this function

``` r
add(100, 1)
```

    [1] 101

``` r
add(x=10, y=10)
```

    [1] 20

``` r
add(x=c(100, 1, 100), y=1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector `x` and `y`

``` r
add(x=c(100,1), y=c(100, 1))
```

    [1] 200   2

> Q. What if I give three inputs to the function?

``` r
#add(x=c(100,1), y=1, z=1)
```

> Q. What if I give only one input to the add function?

``` r
addnew <- function(x, y=1) {
  x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write our function with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input arguments “default”values by setting them equal to
some sensible value - e.g. y=1 in the `addnew()` function.

## A second function

Let’s try something more interesting: Make a sequence generating tool..

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size=4)
```

    [1]  9  8  5 10

> Q. Generate 9 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size=9)
```

    [1]  9  6  8  3  2  5 10  1  4

> Q. Generate 12 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size=12, replace=T)
```

     [1]  2  6  3  5  1  9 10  5 10  9  6  2

> Q. Write code for the `sample()` function that generates nucleotide
> sequences of length 6?

``` r
sample(c("A","C","G","T"), size = 6, replace=T)
```

    [1] "G" "C" "A" "C" "G" "C"

> Q. Write a first function `generated_dna()` that returns a **user
> specified length** DNA sequence:

``` r
generate_dna <- function(len=6) {
  sample(c("A","C","G","T"), size = len, replace=T)
}
```

``` r
generate_dna()
```

    [1] "A" "A" "C" "A" "G" "G"

> **Key-Points** Every function in R looks fundamentally the same in
> terms of its structure. Basically 3 things: name, input, and body

    name <- function(input) {
    body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments. With optional arguments having a
> set default value.

> Q. Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standard format like “AGTAGTA” rather
> than the vector “A”, “C”, “G”, “A”

``` r
generate_dna <- function(len=6, fasta=TRUE) {
  
  ans <- sample(c("A","C","G","T"), 
         size = len, replace=T)
 if(fasta) {  
   cat("Single-element vector output")
 ans <-  paste(ans, collapse = "")
 } else {
  cat("Multi-element vector output")
}
 return(ans)
}
```

``` r
generate_dna(fasta=T)
```

    Single-element vector output

    [1] "TCGTGC"

The `paste()` function - it’s job is to join up or stick together
(a.k.a. paste)input strings together

``` r
paste("alice", "loves R", sep=" ")
```

    [1] "alice loves R"

Flow control means where the R brain goes in your code

``` r
good_mood <- TRUE

if(good_mood) {
  cat("Great!")
} else{
  cat("Bummer")
}
```

    Great!

## A Protein generating function

> Q. Write a function, called `generate_protein()`, that generates a
> user specified length protein sequence.

> Q. Use that function to generate random protein sequences between
> length 6 and 12

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

There are 20 natural amino acids

``` r
aa <- c("A","R","N","D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
```

``` r
generate_protein <- function(len) {
  
  # The amino-acids to sample from
  aa <- c("A","R","N","D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  #Draw n=len amino-acids to make our sequence
  ans <- sample(aa, size = len, replace=T)
  ans <- paste(ans, collapse = "")
  return(ans)
}
```

``` r
myseq <- generate_protein(42)
myseq
```

    [1] "SYNFTEVMNREGLPADYDCMDMFKWRRQVCIFWMTSFNERHG"

> Q. Use that function to generate random protein sequences between
> length 6 and 12

``` r
generate_protein(6)
```

    [1] "CNLPVM"

``` r
generate_protein(7)
```

    [1] "KRPLCLH"

``` r
generate_protein(8)
```

    [1] "VAMDYGAV"

``` r
generate_protein(9)
```

    [1] "NSWPYIIGW"

``` r
generate_protein(10)
```

    [1] "WRPLRPCIVH"

``` r
generate_protein(11)
```

    [1] "DRAHIEPFCMR"

``` r
generate_protein(12)
```

    [1] "EHIKTIDMQMTE"

``` r
for(i in 6:12) {
  # FASTA ID line ">id"
  cat(">", i, sep="", "\n")
  # Protein sequence line
  cat(generate_protein(i), "\n")
}
```

    >6
    HEENHC 
    >7
    CEARGHA 
    >8
    PQNIWNEM 
    >9
    MHVFEYGSI 
    >10
    SENYGWGRPM 
    >11
    KYLWAMPWAHK 
    >12
    AIKCERLDPRPA 

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

Yes, after 8 sequences and above they are unique.
