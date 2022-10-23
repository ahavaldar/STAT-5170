##########################
#### STAT 5170 Lab 4: ####
##########################


#### Question 3C: ####

n = 10
x.ts = c(5.8, 4.4, 6.2, 6.0, 6.9, 6.6, 7.4, 13.5, 16.3, 11.8)

## Step 0: Make a guess for mu! ##

mu.guess = mean(x.ts)


## Step 1:  Substitute mu.guess into the formula for phi1 derived in Part B i)! ##

phi.1 = (sum((x.ts[2:10]-mu.guess) * (x.ts[1:9]-mu.guess))) / (sum((x.ts[1:9]-mu.guess)^2))


## Step 2:  Substitute phi.1 into the formula for mu derived in Part B ii)! ##

mu = (1/(n-1)) * (((phi.1*sum(x.ts[1:9])) - sum(x.ts[2:10])) / (phi.1-1))

## Step 3:  Repeat Steps 1 and 2 until the numbers do not noticeably change ##

phi1.vect = numeric()
mu.vect = numeric()

for(i in 1:15){
  phi.1 = (sum((x.ts[2:10]-mu) * (x.ts[1:9]-mu))) / (sum((x.ts[1:9]-mu)^2))
  mu = (1/(n-1)) * (((phi.1*sum(x.ts[1:9])) - sum(x.ts[2:10])) / (phi.1-1))
  phi1.vect[i] = phi.1
  mu.vect[i] = mu
}

cbind(phi1.vect, mu.vect)

#### Estimator for sigw2 term: ####

sigw2 = (1/n) * sum(((x.ts[2:10]-mu.vect[15])-(phi1.vect[15]*(x.ts[1:9]-mu.vect[15])))^2)
sigw2
