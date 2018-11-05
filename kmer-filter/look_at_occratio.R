
my.dat <-
  read.table("Data/161010_Chinese_Spring_v1.0_pseudomolecules.occratio.10.45.fix")

my.dat <-
  read.table("Data/Hordeum_vulgare.Hv_IBSC_PGSB_v2.dna.toplevel.occratio.10.45.fix")

head(my.dat, 6)
nrow(my.dat)

names(my.dat) <-
  c("k",
    
    ## unique ....... # distribution of unique mers
    "uniq.n",
    "u.f",
    
    ## nonunique .... # distribution of non unique mers (counting each
    ##                  non unique mer only once)
    "rept.n",
    "r.f",
    
    ## nonuniquemulti # distribution of non unique mers (counting each
    ##                  non unique mer more than once)
    "rept.mn",
    "r.mf",
    
    ## total ........ # distribution of all mers (counting each non
    ##                  unique mer only once)
    "all.n",

    ## relative ..... # distribution of all mers (counting each non
    ##                  unique mer more than once)
    "all.mn"
    )

head(my.dat, 10)



## Plot total number of distinct kmers with k

## This plot includes kmers that occur once or multiple times, but
## counts each unique mer only once regardless
plot(all.n ~ k, data = my.dat,
     main="Number of distinct kmers, unique or not, against k")

## General observations, as k grows, there is a rapid growth in number
## of kmers. There is a fast growth in the number of kmers at a
## certain critical size.

## Is this the kmer size after which growth is 'stable'?
abline(v=21)



## Plot total number of uniqe kmers with k

## This plot only includes kmers that occur exactly once
plot(uniq.n ~ k, data = my.dat,
     main="Number of distinct unique kmers against k")

abline(v=21)

## General observations, at low values of k, all kmers are repetative
## (occur more than once).



## Both

matplot(x=my.dat$k, y=my.dat[,c('all.n', 'uniq.n')],
        type='b', pch=c(1,2))

abline(v=21)


## Ratio 
plot(uniq.n / all.n ~ k, data = my.dat)

abline(v=21)
abline(h=.8)
abline(h=.75)
abline(v=31)

abline(v=16, col=2)
abline(v=17, col=3)
abline(v=18, col=4)
abline(v=19, col=5)
abline(v=20, col=6)
abline(v=21, col=7)

plot(rept.n / all.n ~ k, data = my.dat,
     main="Fraction of distinct unique kmers against distinct unique kmers against k")










## Plot fraction of total number of uniqe kmers with k
plot(u.f ~ k, data = my.dat,
     ylim=c(0,1), main="Fraction of kmers unique with k")
##
abline(v=19, col=2)
abline(v=41, col=2)
abline(h=0.80, col=5)
abline(h=0.86, col=2, lwd=3)

plot(r.f ~ k, data = my.dat, type="b",
     ylim=c(0,1))
points(u.f ~ k, data = my.dat, type="b")
##
abline(v=19, col=2)
abline(v=41, col=2)
abline(h=0.80, col=5)
abline(h=0.86, col=2, lwd=3)




plot(r.mf ~ k, data = my.dat, type="b",
     ylim=c(0,1)#, main="Fraction of kmers unique with k"
     )


plot(rept.mn ~ k, data = my.dat, type="b"
     )

