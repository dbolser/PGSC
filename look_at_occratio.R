

my.dat <-
  read.table("Data/PGSC0003DMS.occratio.10.45.fix")

my.dat <-
  read.table("Data/PGSC0003DMB.occratio.10.45.fix")

head(my.dat)
nrow(my.dat)

names(my.dat) <-
  c("k",
    "uniq.mer.n",    "uniq.mer.f",
    "rept.mer.uc.n", "rept.mer.uc.f",
    "rept.mer.dc.n", "rept.mer.dc.f",
    "n.mer.uc",
    "n.mer.al")


plot(uniq.mer.f ~ k, data = my.dat,
     ylim=c(0,1), main="Fraction of kmers unique with k")
##
abline(v=19, col=2)
abline(v=41, col=2)
abline(h=0.86, col=2)
abline(h=0.95, col=2)


plot(rept.mer.uc.f ~ k, data = my.dat, type="b",
     ylim=c(0,1))
points(uniq.mer.f ~ k, data = my.dat, type="b")
##
abline(v=19, col=2)
abline(v=41, col=2)
abline(h=0.86, col=2)
abline(h=0.95, col=2)



plot(rept.mer.dc.f ~ k, data = my.dat, type="b",
     ylim=c(0,1)#, main="Fraction of kmers unique with k"
     )

