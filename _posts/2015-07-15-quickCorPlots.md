---
layout: post
title:  Quick correlation scatter plots with p values
date:   2015-07-15
---

This function is quickly generating scatter plots with a correlation lines and p-values. 

#### The annotated scatter function:

```r
plotAnnotatedScatter <- function(x, y, pointCol=rgb(0,0,0,0.7), 
                                 legendPos="topleft", legendCex=1, ... ){
        # Generate a linear model summary
        fit <- lm(y ~ x)
        fitSum <- summary(fit)
        r2 <- fitSum$r.squared
        pVal <- fitSum$coefficients[2,4]
        
        # Format the legend for r and p values
        rp <- vector('expression',2)
        rp[1] <- substitute(expression(italic(R)^2 == valueA), 
                            list(valueA = format(r2,dig=3)))[2]
        rp[2] <- substitute(expression(italic(p) == valueB), 
                            list(valueB = format(pVal, digits = 2)))[2]
        
        # Plot the data
        plot(x, y, pch=19, cex=1.2, col=pointCol,...)
        
        # Add line for linear model fit
        abline(fit)
        
        # Add the legend
        legend(legendPos, legend = rp, bty = 'n', cex=legendCex)
}
```

#### First, generate some test data

```r
x <- rnorm(n=20, mean=7, sd=3)
y <- x+rnorm(20)
```

#### View the plot

```r
plotAnnotatedScatter(x, y)
```

![](/images/corPlots/corPlot1-1.png) 

#### Add axis labels

```r
plotAnnotatedScatter(x, y, xlab="X variable", ylab="Y variable")
```

![](/images/corPlots/corPlot2-1.png) 

#### Change point colour

```r
plotAnnotatedScatter(x, y, pointCol="red")
```

![](/images/corPlots/corPlot3-1.png) 

#### Change the legend font size

```r
plotAnnotatedScatter(x, y, legendCex=0.7)
```

![](/images/corPlots/corPlot4-1.png) 

#### Change legend position

```r
plotAnnotatedScatter(x, y, legendPos="bottomright")
```

![](/images/corPlots/corPlot5-1.png) 


