df <- data.frame(
    x=1:10,
    y=0
)

fit <- quantreg::rq(y~I(x)+I(x^2), data=df, tau=0.9)
set.seed(395)
ss <- summary(fit, se="boot")

df2 <- structure(list(infected = c(0.016, 0.376, 0.424, 0.32, 0.272, 
    0.464, 0, 0.016, 0, 0.024), sexual = c(0, 0, 0, 0, 0, 0, 0, 0.064, 
    0, 0)), .Names = c("x", "y"), row.names = c(NA, -10L
    ), class = "data.frame")

fit2 <- quantreg::rq(y~I(x)+I(x^2), data=df2, tau=0.9)
set.seed(320)
ss2 <- summary(fit2, se="boot")
